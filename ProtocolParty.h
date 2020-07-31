// Hanjun Li <lihanjun1212@gmail.com>

#ifndef PROTOCOLPARTY_H_
#define PROTOCOLPARTY_H_

#include <stdlib.h>
// #include "ProtocolTimer.h"
// #include <chrono>
// #include <libscapi/include/infra/Measurement.hpp>
// #include <libscapi/include/primitives/Matrix.hpp>
#include <bitset>
#include <emmintrin.h>
#include <libscapi/include/infra/Common.hpp>

#include <cmath>
#include <algorithm>
#include <thread>
#include <vector>
#include <iostream>
#include <fstream>
#include "Interpolate.h"
#include "Dispute.h"
#include "Communication.h"
#include "Tape.h"
#include "Share.h"
#include "HIM.h"
#include <libscapi/include/primitives/Mersenne.hpp>
#include <libscapi/include/cryptoInfra/Protocol.hpp>
#include <libscapi/include/comm/MPCCommunication.hpp>

#include <libscapi/include/circuits/ArithmeticCircuit.hpp>

#define flag_print_output false
#define msgPrefixLen 3                    // 3 byte + 2 _fieldByteSize

using namespace std;
// using namespace std::chrono;

// TODO: baseSharing

// TODO: minor: add output msg during dispute
// TODO: minor: make compatible w/ EMP(?)
// TODO: minor: some party does not follow protocol? (currently hangs)
// improve: use consistent full shares name and format
// improve: refactor matrix mult class
// improve: use matrix mult library
// improve: use move semantics(?)
// improve: time each part
template <class FieldType>
class ProtocolParty : public Protocol{

private:
  // -- global data
  int _numThreads = 4;     // at least 1
  int _K = 7;              // interpolation degree <=> 'shrink' factor
  int _iterations;
  string _inputFile;
  ArithmeticCircuit _circuit;
  TemplateField<FieldType> *_field;
  FieldType _zero;
  int _N, _T, _myId;            // N parties; T malicious
  int _nGates, _nInputs, _nOutputs, _nMults;
  int _fieldByteSize, _segSize, _nSegs, _authKeySize, _dispMsgSize;
  
  // -- communication channels, dispute sets, tapes
  Communication _comm;
  Dispute _disp;
  PartyTape<FieldType> _protoTape;
    
  // -- various matrices set in constructor
  vector<FieldType> _alpha_k;    // K distinct non-zero field elements
  vector<FieldType> _alpha_2km1; // 2K-1 distinct non-zero field element
  HIMp<FieldType> _mat_use_n;    // all _N points to 0
  vector<HIM<FieldType>> _mat_use_n_twist_vec; // all _N points (and 0) to i
  HIMp<FieldType> _mat_check_tp1;  // first t+1 points to rest N-t-1
  HIMp<FieldType> _mat_check_km1;  // first k points to rest 2k-1 - k
  VDMt<FieldType> _mat_vand_trans; // for extracting random shares
  vector<FieldType> _buffer;       // for interp to a sinlge element

  // -- local storage: global
  vector<Share<FieldType>> _gateShares; // my share of each gate
  vector<vector<FieldType>> _keys; // [i]-auth keys for nonDips (party j)
  vector<vector<vector<FieldType>>> _reckeyShares; // [i]-auth key shares received from j of k
  vector<vector<vector<FieldType>>> _dealtKeyShares; // [i]-auth key shares sent to j of k

  // -- local storage: per seg
  vector<Share<FieldType>> _coinShares;    // for public random coins
  int _coinOffset = 0;                     // next share for coin
  vector<Share<FieldType>> _padShares;     // for padding transcripts
  int _padOffset = 0;                      // next share for padding
  vector<Share<FieldType>> _refreshShares; // used for refresh only
  int _refreshOffset = 0;                  // next share for refresh
  vector<Share<FieldType>> _doubleShares;  // used for DN mult only
  int _doubleOffset = 0;        // next random double share to use

  HIMp<FieldType> _mat_to_T;    // for calculating share of T set
  HIMp<FieldType> _mat_from_T;  // for calculating share of nonT set
  vector<HIM<FieldType>> _mat_twist_vec;
  vector<bool> _TMask;

  // -- simple local helper functions
  __inline__
  void segMsg(int segStart, int segEnd, string msg) {
    cout << "[" << _myId << "] Seg: "
         << segStart << " - " << segEnd << ": "
         << msg  << endl;
  }

  __inline__
  FieldType getSecret(vector<FieldType>& fullShare) {
    assert(fullShare.size() == _N);
    _mat_use_n.MatrixMult(fullShare, _buffer);
    return _buffer[0];
  }

  __inline__
  FieldType getTwistSecret(vector<FieldType>& fullShare, int target) {
    assert(fullShare.size() == _N);
    vector<FieldType> inShares(_N, _zero);
    vector<int> inIds = allNButOne(target);
    for (int i=0; i<_N-1; i++) {
      int pid = inIds[i];
      inShares[i] = fullShare[pid];
    }
    _mat_use_n_twist_vec[target].MatrixMult(inShares, _buffer);
    return _buffer[0];
  }

  __inline__
  void getPadShares(int nRands, vector<Share<FieldType>> &result) {
      result.assign(_padShares.begin() + _padOffset,
                    _padShares.begin() + _padOffset + nRands);
      _padOffset += nRands;
  }

  __inline__
  void bytesToVector(vector<byte>& byteVec,
                     vector<FieldType>& elemVec) {
    int nElem = byteVec.size() / _fieldByteSize;
    elemVec.resize(nElem);
    for (int i=0; i<nElem; i++) {
      elemVec[i] = _field->bytesToElement(byteVec.data() + i * _fieldByteSize);
    }
  }

  __inline__
  bool checkTShares(vector<FieldType> &fullShares) {
    vector<FieldType> restTShares(_N-1-_T); // = T
    vector<FieldType> TShares(_T+1);
    if (_disp.hasCorrupt()) {
        vector<int> NonTSet, TSet;
        _disp.tAndNonTSetP(_myId, TSet, NonTSet);
        for (int i=0; i < _T+1; i++) {
          int pt = TSet[i];
          TShares[i] = fullShares[pt];
        }
        _mat_from_T.MatrixMult(TShares, restTShares);
        for (int i=0; i<_T; i++) {
          int pnt = NonTSet[i];
          // NOTE: also checking corr's share being zero
          // becuase fullShares[corr] is zero
          if (restTShares[i] != fullShares[pnt]){
            return false;
          }
        }
    } else {
      TShares.assign(fullShares.begin(), fullShares.begin() + _T+1);
      _mat_check_tp1.MatrixMult(TShares, restTShares);
      for (int i = 0; i < _N-_T-1; i++) {
        if (restTShares[i] != fullShares[_T + 1 + i]) {
          return false;
        }
      }
    }
    return true;
  }

  __inline__
  FieldType getSecretCheck(vector<FieldType> &fullShares) {
    if (!checkTShares(fullShares)) {
      cout << "inconsistent shares" << endl;
      abort();
    }
    return getSecret(fullShares);
  }

  __inline__
  void readMyInputs(int nInputs, vector<FieldType>& myInputs) {
    ifstream inFile(_inputFile);
    myInputs.resize(nInputs, _zero);
    long curInput;
    for (int i=0; i<nInputs; i++) {
      if (inFile.eof()) {
        cout << "input file too short: gets "
             << i+1 << "; needs " << nInputs;
        break;
      }
      inFile >> curInput;
      myInputs[i] = _field->GetElement(curInput);
    }
    inFile.close();
  }

  __inline__
  void makeTwist2TOShare(int target, vector<FieldType>& fullShare) {
    fullShare.assign(_N, _zero);
    vector<FieldType> inShare(_N, _zero);
    vector<int> inIds = allNButOne(_myId); // use my pos as the free point
    for (int i=0; i<_N-1; i++) {
      int pid = inIds[i];
      if (pid == target || _disp.isDisp(_myId, pid)) {
        continue;               // target pos fixed to zero
      }
      inShare[i] = _field->Random();
      fullShare[pid] = inShare[i];
    }
    _mat_use_n_twist_vec[_myId].MatrixMult(inShare, _buffer);
    fullShare[_myId] = _buffer[0];
  }

  __inline__
  FieldType makeTwist2TShare(int target, vector<FieldType>& fullShare) {
    fullShare.assign(_N, _zero);
    vector<FieldType> inShare(_N, _zero);
    vector<int> inIds = allNButOne(target);
    for (int i=0; i<_N-1; i++) {
      int pid = inIds[i];
      if (_disp.isDisp(_myId, pid)) {
        continue;
      }
      inShare[i] = _field->Random();
      fullShare[pid] = inShare[i];
    }
    _mat_use_n_twist_vec[target].MatrixMult(inShare, _buffer);
    return _buffer[0];
  }

  __inline__
  void makeTwistTShares(vector<FieldType>& secrets, vector<int>& targets,
                        vector<vector<FieldType>>& fullSharesVec) {
    int nShares = secrets.size();
    assert(nShares == targets.size());
    fullSharesVec.clear();
    fullSharesVec.resize(_N, vector<FieldType>(nShares));
    vector<FieldType> dispShares(_T+1, _zero), nonDispShares(_T+1);
    vector<int> TwistSet;
    for (int i=0; i<nShares; i++) {
      int target = targets[i];
      assert(!_disp.isDisp(_myId, target));
      int nDisp = _disp.twistSet(_myId, target, TwistSet);
      dispShares[_T] = secrets[i];
      for (int j=nDisp; j<_T-1; j++) {
        int pid = TwistSet[j];
        dispShares[j] = _field->Random();
        fullSharesVec[pid][i] = dispShares[j];
      }
      _mat_twist_vec[target].MatrixMult(dispShares, nonDispShares);
      for (int j=_T-1; j<2*_T; j++) {
        int pid = TwistSet[j];
        fullSharesVec[pid][i] = nonDispShares[j - (_T-1)];
      }
    }
  }

  __inline__
  void makeTShares(vector<FieldType>& secrets, vector<vector<FieldType>>& fullSharesVec,
                   bool nonTZero = false) {
    int nShares = secrets.size();
    fullSharesVec.clear();
    fullSharesVec.resize(_N, vector<FieldType>(nShares));    
    vector<int> NonTSet, TSet;
    _disp.tAndNonTSetP(_myId, TSet, NonTSet);
    vector<FieldType> nonTShares(_N-_T), TShares(_T+1);
    for (int i=0; i<nShares; i++) {
      nonTShares[_N-_T-1] = secrets[i]; // NonTSet's end is position 0
      if (!nonTZero) {
        for (int j = 0; j < _N - _T - 1; j++) {
          int pnt = NonTSet[j]; // vvv disp shares to 0
          nonTShares[j] = _disp.isDisp(_myId, pnt) ? _zero : _field->Random();
          fullSharesVec[pnt][i] = nonTShares[j];
        }
      }
      _mat_to_T.MatrixMult(nonTShares, TShares);
      for (int j=0; j<_T+1; j++) {
        int pt = TSet[j];
        fullSharesVec[pt][i] = TShares[j];
      }
    }
  }

  __inline__
  void make2TShares(int nShares, vector<FieldType>& secrets,
                    vector<vector<FieldType>>& fullSharesVec) {
    secrets.resize(nShares);
    fullSharesVec.clear();
    fullSharesVec.resize(_N ,vector<FieldType>(nShares));
    vector<FieldType> doubleShares(_N);
    for (int i=0; i<nShares; i++) {
      for (int j=0; j<_N; j++) {  // vvv disp shares = 0
        doubleShares[j] = _disp.isDisp(_myId, j) ? _zero : _field->Random();
        fullSharesVec[j][i] = doubleShares[j]; 
      }
      secrets[i] = getSecret(doubleShares);
    }
  }

  __inline__
  void intToElements(vector<int>& intVec, vector<FieldType>& elemVec) {
    int nInts = intVec.size();
    elemVec.resize(nInts);
    for (int i=0; i<nInts; i++) {
      elemVec[i] = _field->GetElement(intVec[i]+1);
    }
  }

  __inline__
  void matForTwist(int pid, vector<HIM<FieldType>>& mat_twist_vec, bool inv = false) {
    mat_twist_vec.resize(_N);
    vector<int> DispSet, NonDispSet;
    _disp.dispAndNonDispSet(pid, DispSet, NonDispSet);
    int nDisp = DispSet.size();
    vector<FieldType> DispElems(nDisp), NonDispElems(_N-nDisp);
    intToElements(DispSet, DispElems);
    intToElements(NonDispSet, NonDispElems);
    for (int i=0; i<_N - nDisp; i++) {
      int pid = NonDispSet[i];  // target id
      auto ptrI = NonDispElems.begin()+i;
      vector<FieldType> rest = DispElems;
      rest.insert(rest.end(), NonDispElems.begin(), ptrI);
      rest.insert(rest.end(), ptrI+1, NonDispElems.end());
      vector<FieldType> alpha_disp(rest.begin(), rest.begin()+_T-1);
      vector<FieldType> beta_nonDisp(rest.begin()+_T-1, rest.end());
      alpha_disp.push_back(_zero); // zero position last before end
      alpha_disp.push_back(*ptrI); // secret position in the end
      mat_twist_vec[pid].allocate(_T+1, _T+1, _field);
      if (inv) {
        mat_twist_vec[pid].InitHIMByVectors(beta_nonDisp, alpha_disp);
      } else {
        mat_twist_vec[pid].InitHIMByVectors(alpha_disp, beta_nonDisp);
      }
    }
  }

  __inline__
  void matForRefresh(int pid, HIMp<FieldType>& mat_from_T, HIMp<FieldType>& mat_to_T) {
    vector<int> TSet, NonTSet;
    _disp.tAndNonTSetP(pid, TSet, NonTSet);
    vector<FieldType> TElems, NonTElems;
    intToElements(TSet, TElems);
    intToElements(NonTSet, NonTElems);
    vector<FieldType> alpha_T = TElems;
    alpha_T.insert(alpha_T.end(), NonTElems.begin(), NonTElems.end());
    vector<FieldType> alpha_nonT = NonTElems;
    alpha_nonT.push_back(_zero);
    alpha_nonT.insert(alpha_nonT.end(), TElems.begin(), TElems.end());
    mat_from_T.allocate(_N - _T - 1, _T+1, _field);
    mat_from_T.InitHIMVectorAndsizes(alpha_T, _T+1, _N - _T - 1);
    mat_to_T.allocate(_T+1, _N - _T, _field);
    mat_to_T.InitHIMVectorAndsizes(alpha_nonT, _N - _T, _T+1);
  }

  __inline__
  vector<int> allNButOne(int p) {
    vector<int> result(_N-1);
    int idx = 0;
    for (int i=0; i<_N; i++) {
      if (i != p) {
        result[idx++] = i;
      }
    }
    return result;
  }

  __inline__
  void matFor2TTwist(vector<HIM<FieldType>>& mat_twist_2t_vec) {
    mat_twist_2t_vec.resize(_N);
    for (int i=0; i<_N; i++) {
      vector<int> inIDs = allNButOne(i);
      vector<FieldType> alpha;
      intToElements(inIDs, alpha);
      alpha.push_back(_zero);   // last pos is for zero
      vector<FieldType> beta(1, _field->GetElement(i));
      mat_twist_2t_vec[i].allocate(1, _N, _field);
      mat_twist_2t_vec[i].InitHIMByVectors(alpha, beta);
    }
  }

  __inline__
  void interpTranscripts(HIMp<FieldType>& mat_from_2k_to_lambda,
                         vector<mTranscript<FieldType>>& transcripts,
                         mTranscript<FieldType>& transcript) {
    int length = 2*_K-1;
    assert(transcripts.size() == length);
    // output transcripts values in transposed form
    vector<Share<FieldType>> rVec(length), r2Vec(length), tmpShare(1);
    vector<FieldType> eVec(length), e2Vec(length), zVec(length), tmpVal(1);
    for (int i=0; i<2*_K-1; i++) {
      rVec[i] = transcripts[i]._r;
      r2Vec[i] = transcripts[i]._r2;
      eVec[i] = transcripts[i]._e;
      e2Vec[i] = transcripts[i]._e2;
      zVec[i] = transcripts[i]._z;
    }
    mat_from_2k_to_lambda.MatrixMultShares(rVec, tmpShare);
    transcript._r = tmpShare[0];
    mat_from_2k_to_lambda.MatrixMultShares(r2Vec, tmpShare);
    transcript._r2 = tmpShare[0];
    mat_from_2k_to_lambda.MatrixMult(eVec, tmpVal);
    transcript._e = tmpVal[0];
    mat_from_2k_to_lambda.MatrixMult(e2Vec, tmpVal);
    transcript._e2 = tmpVal[0];
    mat_from_2k_to_lambda.MatrixMult(zVec, tmpVal);
    transcript._z = tmpVal[0];
  }

  __inline__
  void interpRelayTranscripts(HIMp<FieldType>& mat_from_2k_to_lambda,
                              vector<vector<FieldType>>& relayTranscripts,
                              vector<FieldType>& relayTranscript) {
    assert(relayTranscripts.size() == 2*_K-1);
    for (int i=0; i<2*_K-1; i++) {
      assert(relayTranscripts[i].size() == _N);
    }
    relayTranscript.resize(_N);
    vector<FieldType> tmp(2*_K-1);
    for (int i=0; i<_N; i++) {
      for (int j=0; j<2*_K-1; j++) {
        tmp[j] = relayTranscripts[j][i];
      }
      mat_from_2k_to_lambda.MatrixMult(tmp, _buffer);
      relayTranscript[i] = _buffer[0];
    }
  }

  __inline__
  void interpKingTranscripts(HIMp<FieldType>& mat_from_2k_to_lambda,
                             vector<vector<vector<FieldType>>>& kingTranscripts,
                             vector<vector<FieldType>>& kingTranscript) {
    assert(2 == kingTranscripts.size());
    for (int i=0; i<2; i++) {
      assert(kingTranscripts[i].size() == 2*_K-1);
      for (int j=0; j<2*_K-1; j++) {
        assert(kingTranscripts[i][j].size() == _N);
      }
    }
    kingTranscript.resize(2);
    vector<FieldType> tmp1(2*_K-1), tmp2(2*_K-1);
    for (int i=0; i<_N; i++) {
      for (int j=0; j<2*_K-1; j++) {
        tmp1[j] = kingTranscripts[0][j][i];
        tmp2[j] = kingTranscripts[1][j][i];
      }
      mat_from_2k_to_lambda.MatrixMult(tmp1, _buffer);
      kingTranscript[0][i] = _buffer[0];
      mat_from_2k_to_lambda.MatrixMult(tmp2, _buffer);
      kingTranscript[1][i] = _buffer[0];
    }
  }

  __inline__
  void
  appendKingTranscripts(vector<vector<vector<FieldType>>>& kingTranscripts1,
                        vector<vector<vector<FieldType>>>& kingTranscripts2) {
    assert(kingTranscripts1.size() == 2);
    assert(kingTranscripts2.size() == 2);
    for (int i=1; i<2; i++) {
      kingTranscripts1[i].size() == kingTranscripts1[0].size();
      kingTranscripts2[i].size() == kingTranscripts2[0].size();
    }
    kingTranscripts1[0].
      insert(kingTranscripts1[0].end(),
             kingTranscripts2[0].begin(), kingTranscripts2[0].end());
    kingTranscripts1[1].
      insert(kingTranscripts1[1].end(),
             kingTranscripts2[1].begin(), kingTranscripts2[1].end());
  }

  __inline__
  void VDMMultShares(unsigned int idxStart, VDMt<FieldType>& matrix,
                     vector<FieldType>& inVals,
                     vector<Share<FieldType>> outShares) {
    int nRows = matrix.getNRows();
    int nCols = matrix.getNCols();
    assert(inVals.size() == nCols);
    vector<FieldType> outVals(nRows);
    vector<unsigned int> outIdxs(nCols);
    for (int i=0; i<nCols; i++) {
      outIdxs[i] = idxStart + i;
    }
    matrix.MatrixMult(inVals, outVals, _N-_T);
    outShares.resize(nRows);
    for (int i=0; i<nRows; i++) {
      outShares[i]._value = outVals[i];
      outShares[i]._idxs = outIdxs;
      matrix.getRow(i, outShares[i]._coeffs);
    }
  }

  __inline__
  void makeOkMsg(vector<byte>& msg) {
    // for checkRand
    msg.resize(_dispMsgSize);
    msg[0] = 0;                 // 0: okay
  }

  __inline__
  void makeOpenMsg(vector<byte>& msg, int dealerId, int tid,
                   FieldType dealerShare, FieldType tShare) {
    // for checkRand
    msg.resize(_dispMsgSize);
    msg[0] = 1; msg[1] = dealerId; msg[2] = tid;      // 1: open(d, t)
    _field->elementToBytes(msg.data() + msgPrefixLen, dealerShare);
    _field->elementToBytes(msg.data() + msgPrefixLen + _fieldByteSize, tShare);
  }

  __inline__
  void readOpenMsg(vector<byte>& msg, int& dealerId, int& tid,
                   FieldType& dealerShare, FieldType& tShare) {
    assert(msg.size() == _dispMsgSize);
    assert(msg[0] == 1);
    dealerId = msg[1]; tid = msg[2];
    dealerShare = _field->bytesToElement(msg.data() + msgPrefixLen);
    tShare = _field->bytesToElement(msg.data() + msgPrefixLen + _fieldByteSize);
  }

  __inline__
  void makeAccuseMsg(vector<byte>& msg, int id) {
    // for checkRand
    msg.resize(_dispMsgSize);
    msg[0] = 2; msg[1] = id;    // 2: accuse(i)
  }

  __inline__
  void makeSuspectMsg(vector<byte>& msg, int dealerId) {
    // for localize
    msg.resize(_dispMsgSize);
    msg[0] = 3; msg[1] = dealerId; // 3: suspect dealer(i)
  }

  __inline__
  void makeBadChunkMsg(vector<byte>& msg, int chunkIdx) {
    // for pairlocateBatch
    msg.resize(_dispMsgSize);
    msg[0] = 4; msg[1] = chunkIdx; // 4: badChunk(i)
  }

  __inline__
  void makeRejectMsg(vector<byte>& msg) {
    // for verifier in activeDealer
    msg.resize(_dispMsgSize);
    msg[0] = 5;                 // 5: reject
  }

    __inline__
    void makeOpenMsgTag(vector<byte>& msg, int dealerId, int tid, int idx,
                        FieldType dealerShare, FieldType tShare) {
    // for checkRand
    msg.resize(_dispMsgSize);
    msg[0] = 6;                 // 6: open(d, t, batchIdx)
    msg[1] = dealerId; msg[2] = tid; msg[3] = idx; 
    _field->elementToBytes(msg.data() + msgPrefixLen, dealerShare);
    _field->elementToBytes(msg.data() + msgPrefixLen + _fieldByteSize, tShare);
  }

  __inline__
  void readOpenMsgTag(vector<byte>& msg, int& dealerId, int& tid, int& idx,
                      FieldType& dealerShare, FieldType& tShare) {
    assert(msg.size() == _dispMsgSize);
    assert(msg[0] == 6);
    dealerId = msg[1]; tid = msg[2]; idx = msg[3];
    dealerShare = _field->bytesToElement(msg.data() + msgPrefixLen);
    tShare = _field->bytesToElement(msg.data() + msgPrefixLen + _fieldByteSize);
  }


  __inline__
  void checkDealerRandShares(vector<byte>& kingMsg,
                             vector<vector<byte>>& recBufsD,
                             vector<FieldType>& secrets,
                             vector<vector<FieldType>>& dealerFullShares) {
    // for checkRand and checkDoubleRand
    // TODO: minor: make sure checks cover everything
    secrets.assign(_N, _zero);
    dealerFullShares.resize(_N, vector<FieldType>(_N, _zero));
    for (int dealerId = 0; dealerId < _N; dealerId++) {
      bytesToVector(recBufsD[dealerId], dealerFullShares[dealerId]);
      // -- check 1: are shares consistent?
      if (!checkTShares(dealerFullShares[dealerId])) {
        makeAccuseMsg(kingMsg, dealerId);
        return;
      }
      // -- check 2: shares for Disp(dealer) is zero
      for (int partyId=0; partyId<_N; partyId++) {
        if (_disp.isDisp(dealerId, partyId) &&
            dealerFullShares[dealerId][partyId] != 0) {
          makeAccuseMsg(kingMsg, dealerId);
          return;
        }
      }
      secrets[dealerId] = getSecret(dealerFullShares[dealerId]);
    }
    makeOkMsg(kingMsg);
  }

  __inline__
  void checkDealerRand2TShares(vector<byte>& kingMsg,
                               vector<vector<byte>>& recBufsD,
                               vector<FieldType>& secrets,
                               vector<vector<FieldType>>& dealerFullShares) {
    // for checkDoubleRand
    // TODO: minor: make sure checks cover everything
    dealerFullShares.resize(_N, vector<FieldType>(_N, _zero));
    for (int dealerId = 0; dealerId < _N; dealerId++) {
      bytesToVector(recBufsD[dealerId], dealerFullShares[dealerId]);
      // -- check 1: are shares consistent? N/A
      // -- check 2: shares for Disp(dealer) is zero
      for (int partyId=0; partyId<_N; partyId++) {
        if (_disp.isDisp(dealerId, partyId) &&
            dealerFullShares[dealerId][partyId] != 0) {
          makeAccuseMsg(kingMsg, dealerId);
          return;
        }
      }
      // -- check 3: 2t-shares agree w/ t-shares?
      if (secrets[dealerId] != getSecret(dealerFullShares[dealerId])) {
        makeAccuseMsg(kingMsg, dealerId);
        return;
      }
    }
    makeOkMsg(kingMsg);
  }

  __inline__
  void checkDealerAndT(vector<byte>& kingMsg, vector<int>& TSet,
                       vector<vector<byte>>& recBufsT,
                       vector<vector<FieldType>>& dealerFullShares) {
    // for checkRand and checkDoubleRand
    kingMsg.resize(_dispMsgSize);
    for (int tid : TSet) {
      for (int j=0; j<_N; j++) {
        FieldType jToTShare =
          _field->bytesToElement(recBufsT[tid].data() + (j * _fieldByteSize));
        if (jToTShare != dealerFullShares[j][tid]) {
          makeOpenMsg(kingMsg, j, tid, dealerFullShares[j][tid], jToTShare);
          return;
        }
      }
    }
    makeOkMsg(kingMsg);
    return;
  }

  __inline__
  void checkPartyDecompShares(vector<byte>& kingMsg, // out 
                              vector<vector<FieldType>>& partyFullShares, // out
                              vector<vector<byte>>& recBufsP, // in
                              vector<FieldType>& zFullShare) {   // in
    // for localize
    partyFullShares.resize(_N);
    for (int partyId=0; partyId<_N; partyId++) {
      FieldType partyZSum = _zero;
      bytesToVector(recBufsP[partyId], partyFullShares[partyId]);
      for (int j=0; j<_N; j++) {
        partyZSum += partyFullShares[partyId][j];
      }
      if (partyZSum != zFullShare[partyId]) {
        makeAccuseMsg(kingMsg, partyId);
        return;
      }
    }
    for (int dealerId=0; dealerId < _N; dealerId++) {
      vector<FieldType> fullShare(_N);
      for (int j=0; j<_N; j++) {
        fullShare[j] = partyFullShares[j][dealerId];
      }
      if (!checkTShares(fullShare)) {
        makeSuspectMsg(kingMsg, dealerId);
        return;
      }
    }
    return;
  }

  __inline__
  void checkDealerBatches(int accused, vector<byte>& verifierMsg,
                          vector<byte>& recBufV, vector<FieldType>& lambdaWs,
                          vector<vector<FieldType>>& myNus, FieldType recTau,
                          vector<FieldType>& recSigma) {
    // for flTag
    vector<FieldType> tauAndBatch;
    bytesToVector(recBufV, tauAndBatch);
    vector<FieldType> calcBatch(_authKeySize, _zero);
    FieldType calcTau = _zero;
    FieldType lambdaD = *(_field->GetOne());
    FieldType calcNu = _zero;
    for (int dealer=0; dealer<_N; dealer++) {
      for (int keyIdx=0; keyIdx < _authKeySize; keyIdx++) {
        calcBatch[keyIdx] +=
          tauAndBatch[_N*(dealer+1) + keyIdx] * lambdaD * lambdaWs[dealer];
      }
      calcTau += tauAndBatch[dealer] * lambdaD * lambdaWs[dealer];
    }
    bool error = false;
    if (calcTau != recTau) {
      makeAccuseMsg(verifierMsg, accused);
      return;
    } 
    for (int keyIdx=0; keyIdx < _authKeySize; keyIdx++) {
      if (calcBatch[keyIdx] != recSigma[keyIdx]) {
        makeAccuseMsg(verifierMsg, accused);
        return;
      }
    }
    for (int dealer=0; dealer<_N; dealer++) {
      FieldType dealerTag = myNus[dealer][accused];
      for (int keyIdx=0; keyIdx < _authKeySize; keyIdx++) {
        dealerTag += _keys[keyIdx][accused] * tauAndBatch[_N*(dealer+1)+ keyIdx];
      }
      if (tauAndBatch[dealer] != dealerTag) {
        makeAccuseMsg(verifierMsg, dealer);
        return;
      }
    }
  }

  __inline__
  bool testDealerOSigma(int dealer, int accused,
                        vector<vector<FieldType>>& oSigmaDealerFullShares) {
    // for flSingleDealer
    // check 1: dealer oShares valid
    vector<FieldType> oFullShare(_N);
    for (int party=0; party<_N; party++) {
      oFullShare[party] = oSigmaDealerFullShares[party][0];
    }
    if (getTwistSecret(oFullShare, accused) != _zero) {
      return false;
    }
    
    // check 2: dealer sigmaShares
    vector<FieldType> sigmaFullShare(_N);
    for (int keyId=0; keyId<_N; keyId++) {
      for (int party=0; party<_N; party++) {
        sigmaFullShare[party] = oSigmaDealerFullShares[party][1+keyId];
        if (_disp.isDisp(dealer, party) && sigmaFullShare[party] != _zero) {
          return false;
        }
      }
      if (!checkTShares(sigmaFullShare)) {
        return false;
      }
    }
    return true;
  }

  __inline__
  void checkDealerTagComp(int dealer, int accused, vector<byte> verifierMsg,
                          vector<vector<FieldType>>& tauOSigmaFullShares,
                          vector<vector<FieldType>>& oSigmaDealerFullShares,
                          vector<FieldType>& tauAccusedFullShare,
                          vector<FieldType>& nuMyFullShare) {
    // -- check if received messages match
    for (int party=0; party<_N; party++) {
      // accused vs. party
      FieldType partyTauShare = tauOSigmaFullShares[party][0];
      if (partyTauShare != tauAccusedFullShare[party]) {
        makeOpenMsgTag(verifierMsg, accused, party, -1,
                       tauAccusedFullShare[party], partyTauShare);
        return;
      }
      // dealer vs. party
      FieldType partyOShare = tauOSigmaFullShares[party][1];
      FieldType dealerOShare = oSigmaDealerFullShares[party][0];
      if (partyOShare != dealerOShare) {
        makeOpenMsgTag(verifierMsg, dealer, party, -1,
                       dealerOShare, partyOShare);
        return;
      }
      for (int keyId=0; keyId<_N; keyId++) {
        FieldType partySigmaShare = tauOSigmaFullShares[party][2+keyId];
        FieldType dealerSigmaShare = oSigmaDealerFullShares[party][1+keyId];
        if (partySigmaShare != dealerSigmaShare) {
          makeOpenMsgTag(verifierMsg, dealer, party, keyId,
                         dealerSigmaShare, partySigmaShare);
          return;
        }
      }
    }

    // -- cehck if parties followed protocol
    // TODO: minor: make sure checks cover everything
    // check 1: dealer oShares valid
    // check 2: dealer sigmaShares
    if (!testDealerOSigma(dealer, accused, oSigmaDealerFullShares)) {
      makeAccuseMsg(verifierMsg, dealer);
    }
    // 3. check tauShares, using verifier's nu
    for (int party=0; party<_N; party++) {
      FieldType calcTau = tauOSigmaFullShares[party][1] + nuMyFullShare[party];
      for (int keyId = 0; keyId < _N; keyId++) {
        calcTau += tauOSigmaFullShares[party][2+keyId] *
          _dealtKeyShares[keyId][party][accused];
      }

      if (calcTau != tauOSigmaFullShares[party][0]) {
        makeAccuseMsg(verifierMsg, party);
        return;
      }
    }
  }

  __inline__
  void makeBatchSendBufs(vector<vector<byte>>& sendBufs,
                         vector<unsigned int>& batchRecIdxs,
                         int batchIdx) {
    // for activeDealer and corruptedDealer
    int batchSize = batchRecIdxs.size();
    vector<FieldType> accusedBatchVals(batchSize);
    for (int i=0; i<batchSize; i++) {
      accusedBatchVals[i] = _protoTape.getRecShareVal(batchRecIdxs[i]);
    }
    vector<byte> sendBuf(batchSize * _fieldByteSize);
    _field->elementVectorToByteVector(accusedBatchVals, sendBuf);
    sendBufs.assign(_N, vector<byte>(1));
    vector<FieldType> tagVals;
    _protoTape.getTagVals(batchIdx, tagVals);
    for (int verifyId=0; verifyId<_N; verifyId++) {
      _field->elementToBytes(sendBufs[verifyId].data(), tagVals[verifyId]);
      sendBufs[verifyId].insert(sendBufs[verifyId].end(),
                                sendBuf.begin(), sendBuf.end());
    }
  }

  __inline__
  void verifyBatch(vector<byte>& verifierMsg, vector<byte>& recBuf,
                   int accused, int batchIdx, FieldType broadcastVal,
                   vector<unsigned int>& batchRecIdxs,
                   vector<FieldType>& batchCoeffs) {
    // for activeDealer and corruptedDealer
    int batchSize = batchRecIdxs.size();
    vector<FieldType> accusedBatch(batchSize + 1); // position 0 is the tag
    makeOkMsg(verifierMsg);       // default if no error found
    if (_myId != accused) {
      bytesToVector(recBuf, accusedBatch);
      // -- check 1: tag valid
      FieldType recTag = accusedBatch[0];
      FieldType trueTag = _protoTape.getPadVal(accused, batchIdx);
      for (int i=0; i<batchSize; i++) {
        trueTag += _keys[i][accused] * accusedBatch[i+1];
      }
      if (recTag != trueTag) {
        makeRejectMsg(verifierMsg);
      } else {
        // -- check 2: linear comb equals to accusedVal?
        FieldType sum = _zero;
        for (int i=0; i<batchSize; i++) {
          sum += accusedBatch[i+1] * batchCoeffs[i];
        }
        if (sum != broadcastVal) {
          makeRejectMsg(verifierMsg);
        }
      }
    }
  }

public:
  bool hasOffline() override { return true; }
  bool hasOnline() override { return true; }
  ~ProtocolParty();
  ProtocolParty(int argc, char *argv[]);

  // -- local functions
  void initializeMat();
  void prepareE2Shares(int layerStart, int layerEnd, int nMults,
                       vector<FieldType>& e2Shares);
  int processNonMult(int layerStart, int layerEnd);  
  // -- main protocol functions
  void run() override;
  bool evalSeg(int segStart, int segEnd, int nMults);
  bool runOffline(int segStart, int segEnd);
  void makeRandShares(int nRands, vector<Share<FieldType>> &randShares);
  void makeRandDoubleShares(int nRands, vector<Share<FieldType>> &randDoubleShares);
  bool runOnline(int segStart, int segEnd);
  void inputPhase();
  void computationPhase(int segStart, int segEnd);
  int processMult(int layerStart, int layerEnd, int nMults);
  void batchMultSingle(vector<FieldType>& e2Shares, vector<Share<FieldType>>& e1Shares);
  void refresh(vector<FieldType>& xShares);
  bool verificationPhase(int segStart, int segEnd);
  void outputPhase();

  // -- verification protocols
  bool checkRefresh(int segStart, int segEnd);
  void delinearization(int segStart, int segEnd, int king,
                       vector<Share<FieldType>>& aShares, vector<Share<FieldType>>& bShares,
                       mTranscript<FieldType>& transcript, vector<vector<FieldType>>& kingTranscript,
                       vector<FieldType>& relayTranscript);
  int dimensionReduction(int king , vector<Share<FieldType>>& aShares, vector<Share<FieldType>>& bShares,
                         mTranscript<FieldType>& transcript, vector<vector<FieldType>>& kingTranscript,
                         vector<FieldType>& relayTranscript);
  void compress(int king, int groupSize, vector<Share<FieldType>>& aShares, vector<Share<FieldType>>& bShares,
                vector<mTranscript<FieldType>>& transcripts, vector<vector<vector<FieldType>>>& kingTranscripts,
                vector<vector<FieldType>>& relayTranscripts, mTranscript<FieldType>& transcriptNew,
                vector<vector<FieldType>>& kingTranscript, vector<FieldType>& relayTranscript);
  void randomization(int king, vector<Share<FieldType>>& aShares, vector<Share<FieldType>>& bShares,
                     mTranscript<FieldType>& transcript, vector<vector<FieldType>>& kingTranscript,
                     vector<FieldType>& relayTranscript);
  bool checkSingleMult(int king, Share<FieldType>& aShare, Share<FieldType>& bShare,
                       mTranscript<FieldType>& transcript, vector<vector<FieldType>>& kingTranscript,
                       vector<FieldType>& relayTranscript);
  bool verifySharing();
  bool tag();
  void keyDistribution(vector<FieldType>& keys, vector<vector<FieldType>>& reckeyShares,
                       vector<vector<FieldType>>& sentKeyShares);
  void singleTagComp(vector<vector<FieldType>>& authBatchs, vector<int>& authDealers,
                     vector<vector<FieldType>>& recOShares, vector<vector<FieldType>>& recTauShares,
                     vector<vector<FieldType>>& tauShares,  vector<vector<FieldType>>& dealtOShares,
                     vector<vector<FieldType>>& myNuShares);
  bool checkKey();
  bool checkTag(vector<vector<FieldType>>& authBatchs, vector<int>& authDealers,
                vector<vector<FieldType>>& recOShares, vector<vector<FieldType>>& recTauShares,
                vector<vector<FieldType>>& tauShares,  vector<vector<FieldType>>& dealtOShares,
                vector<vector<FieldType>>& myNuShares);
  void baseSharing(vector<vector<FieldType>>& randomShares);
  bool verifyOutput();
  FieldType challenge();

  // -- fault localization protocols
  void analyzeSharing(vector<FieldType>& fullShare, Share<FieldType>& myShare);
  void checkRand(vector<FieldType>& randShare, Share<FieldType>& share);
  void checkDoubleRand(vector<FieldType>& fullRandShare, vector<FieldType>& fullRand2TShare,
                       Share<FieldType>& randShare, Share<FieldType>& fullShare2T);
  int localize(vector<FieldType>& fullShare, Share<FieldType>& myShare,
               vector<FieldType>& dealerFullShare, vector<FieldType>& myFullShare);
  void activeDealer(int dealerId, vector<FieldType>& dealerFullShare,
                    vector<FieldType>& myFullShare, Share<FieldType>& zShare);
  void corruptedDealer(int dealerId, vector<FieldType>& dealerFullShare, Share<FieldType>& zShare);
  int flTag(int verifier, int accused, FieldType lambda, vector<FieldType>& lambdaWs,
            vector<vector<FieldType>>& dealerTaus, vector<vector<FieldType>>& dealerBatchs,
            vector<vector<FieldType>>& dealerNus, FieldType recTau, vector<FieldType>& recBatch);
  void flSingleDealer(int verifier, int accused, int dealer, FieldType lambda,
                      vector<int>& authDealers, vector<vector<FieldType>>& authBatchs, 
                      vector<vector<FieldType>>& recOShares, vector<vector<FieldType>>& recTauShares,
                      vector<vector<FieldType>>& tauShares,  vector<vector<FieldType>>& dealtOShares,
                      vector<vector<FieldType>>& myNuShares);
  // -- helpers
  // for challenge
  void openTShares(int nShares, bool check, vector<FieldType> &Shares, vector<FieldType> &clears);
  // for compress and dimensionReduction
  void buildPolyVecInd(int king, int groupSize, vector<Share<FieldType>>& aShares,
                       vector<Share<FieldType>>& bShares, vector<mTranscript<FieldType>>& trancripts,
                       vector<vector<vector<FieldType>>>& kingTranscripts, vector<vector<FieldType>>& relayTranscripts);
  void buildPolyVecIndWorker(vector<Share<FieldType>>& aShares, vector<Share<FieldType>>& bShares,
                             vector<vector<Share<FieldType>>>& ASharesEval,
                             vector<vector<Share<FieldType>>>& BSharesEval,
                             int groupSize, int threadId);
  void dotProd(int king, int groupSize, vector<Share<FieldType>>& aShares, vector<Share<FieldType>>& bShares,
               vector<mTranscript<FieldType>>& transcripts, vector<vector<vector<FieldType>>>& kingTranscripts,
               vector<vector<FieldType>>& relayTranscripts);
  // for checkRefresh and checkSingleMult
  void getRefreshTranscripts(int segStart, int segEnd, int king, rTranscript<FieldType>& partyTranscript,
                             vector<vector<FieldType>>& kingTranscript);
  void broadcastKingTranscript(int king, vector<vector<FieldType>>& kingTranscript);
  void broadcastPartyTranscript(rTranscript<FieldType>& myTranscript, vector<vector<FieldType>>& allTranscriptVals);
  // for checkRand and checkDoubleRand
  void solveAccuseRand(int king, int accused, vector<byte>& relayBufs);
  void solveAccuseRand2T(int king, int accused, FieldType secret, vector<byte>& relayBuf);
  void solveOpenRand(int king, vector<byte>& openMsg, vector<vector<byte>>& relayBufs,
                     vector<FieldType>& myShareDecomp, vector<FieldType>& myFullShare);
  // for activeDealer and corruptedDealer
  bool pairLocateBatch(int dealer, int accused, FieldType& accusedBraodcast, Share<FieldType>& share);
  bool pairLocateChunk(int dealer, int accused, FieldType& accusedBraodcast, Share<FieldType>& share);
  bool allLocateBatch(vector<FieldType>& prevBraodcasts, Share<FieldType>& share);
  bool allLocateChunk(vector<FieldType>& prevBraodcasts, Share<FieldType>& share);
  // for localize
  void padZShare(Share<FieldType>& zShare, vector<FieldType>& zFullShare);
  void solveAccuseSum(int king, int accused, FieldType broadcastVal, vector<byte>& relayBuf);
  void solveAccuseVal(int king, int accused, int dealer, FieldType broadcastVal, vector<byte>& relayBuf);
  void broadcastShare(FieldType share, vector<FieldType>& fullShare);
  // for flSingleDealer
  void solveAccuseTag(int verifier, int accused, int dealer, int suspect, vector<byte>& relayBuf, vector<vector<byte>>& relayBufs);
  void solveOpenTag(int verifier, int accused, vector<byte> &openMsg,
                    FieldType myTauShare, vector<FieldType>& recTauFullShare,
                    FieldType myOShare, vector<FieldType>& sentOFullShare,
                    vector<FieldType>& mySigmaShares, vector<vector<FieldType>>& sentSigmaFullShares,
                    vector<byte>& relayBuf, vector<vector<byte>>& relayBufs);
  void checkRelayer(int relayer, int id1, int id2, FieldType share1, FieldType rShare);
};

template <class FieldType>
ProtocolParty<FieldType>::ProtocolParty(int argc, char *argv[])
    : Protocol("MPCHonestMajorityNoTriples", argc, argv) {
  _iterations =
    stoi(getParser().getValueByKey(arguments, "internalIterationsNumber"));
  string fieldType = getParser().getValueByKey(arguments, "fieldType");
  if (fieldType.compare("ZpMersenne31") == 0) {
    _field = new TemplateField<FieldType>(2147483647);
  } else if (fieldType.compare("ZpMersenne61") == 0) {
    _field = new TemplateField<FieldType>(0);
  } else {
    cout << "don't use other fields yet " <<
      " (remember to add TemplatedField functions for them)" << endl;
    abort();
  }
  _zero = *(_field->GetZero());
  _N = stoi(getParser().getValueByKey(arguments, "partiesNumber"));
  // TODO: [ask] what to do for even # of parties
  assert(_N % 2 == 1);
  _T = (_N - 1)/ 2;
  _myId = stoi(getParser().getValueByKey(arguments, "partyID"));
  _inputFile = getParser().getValueByKey(arguments, "inputFile");
  // outputFile = getParser().getValueByKey(arguments, "outputFile");
  string circuitFile =
    getParser().getValueByKey(arguments, "circuitFile");
  _circuit.readCircuit(circuitFile.c_str());
  _circuit.reArrangeCircuit();
  _nGates = _circuit.getNrOfGates();
  _nInputs = _circuit.getNrOfInputGates();
  _nOutputs = _circuit.getNrOfOutputGates();
  _nMults = _circuit.getNrOfMultiplicationGates();
  // -- common computations
  _fieldByteSize = _field->getElementSizeInBytes();
  _nSegs = _N*_N;
  _segSize = (_nMults + _nSegs - 1) / _nSegs;
  // NOTE: the facator 10 is flexible, but
  // (_T+1) * 10 should be smaller than 255
  _authKeySize = _segSize / (_T + 1) / 10;
  _dispMsgSize = msgPrefixLen + 2*_fieldByteSize;
  
  // -- communication channels
  string partiesFile =
      getParser().getValueByKey(arguments, "partiesFile");
  _comm.reset(_N, _myId, _numThreads, partiesFile, &_disp);
  _disp.reset(_N);
  _protoTape.reset(_N, _myId, _authKeySize, _field);

  // -- global storage
  _gateShares.resize(_nGates - _nOutputs);
  _buffer.resize(1);
  initializeMat();
}

template <class FieldType> void ProtocolParty<FieldType>::run() {
  const auto& gates = _circuit.getGates();
  
  for (int i = 0; i < _iterations; i++) {
    int segStart = 0;
    int segEnd = 0;
    for (int curSeg = 0; curSeg < _nSegs; curSeg++){
      if (curSeg == 2) {        // TEST: delete this, test only
        _disp.addDispPairs(0, 1);
        _disp.addDispPairs(2, 3);
        _disp.addCorrParty(0);
      }
      
      if (_disp.isCorrupt(_myId)) {
        segMsg(segStart, segEnd, "I'm corrupt, abort.");
        abort();
      }
      // build segment and eval
      segStart = segEnd;
      int numMultsInSeg = 0;
      while (segEnd < _nGates && numMultsInSeg < _segSize) {
        numMultsInSeg += ((gates[segEnd++].gateType == MULT) ? 1 : 0);
      }
      while (!evalSeg(segStart, segEnd, numMultsInSeg)) {
        segMsg(segStart, segEnd, "failed seg, restart.");
      }
    }
  }
}

template <class FieldType>
bool ProtocolParty<FieldType>::evalSeg(int segStart, int segEnd, int nMults) {
  _protoTape.startNewSeg();
  if (_disp.getNewDisp()) {     // update matrices w/ new disp
    // -- prepare mat for normal sharing
    matForRefresh(_myId, _mat_from_T, _mat_to_T);
    // -- prepare mat for twisted sharing
    matForTwist(_myId, _mat_twist_vec);
    _disp.tMaskP(_myId, _TMask);
  }
  
  if (!runOffline(segStart, segEnd)) {
    segMsg(segStart, segEnd, "failed runOffline");
    _protoTape.clearLastSeg();
    return false;
  }
  segMsg(segStart, segEnd, "finished runOffline");
  if (!runOnline(segStart, segEnd)) {
    segMsg(segStart, segEnd, "failed runOnline");
    _protoTape.clearLastSeg();
    return false;
  }
  segMsg(segStart, segEnd, "finished runOnline");
  // if seg finishes successfully, set false.
  // O/w, set to true inside dispute protocol.
  _disp.setNewDisp(false);
  return true;
}
  
template <class FieldType>
bool ProtocolParty<FieldType>::runOffline(int segStart, int segEnd) {
  // ---- # of random single shares:
  // 1. Padding before reveal
  //    -- randomization -> 1 time of [2]
  //    -- verifySharing (input and computation) -> 2 times
  //    -- checkRefresh -> 1 time
  //    -- localize -> 1 time
  //    in total: 5
  // 2. Random Coins
  //    [1]
  //    -- verifySharing after input and computation phases -> 2 times
  //    -- verifyOutput -> 1 time
  //    -- checkRefresh -> 1 time
  //    -- delinearization -> 1 time
  //    -- compress -> nCompressions + 2 times
  //    -- checkKeys -> 1 time (or 0)
  //    -- checkTag -> 1 time (or 0)
  //    in total: nCompress + 9
  // 3. refresh
  //    [1] per Mult
  //    in total: segSize + 1 (in checkRefresh for random pad)
  _coinOffset = 0;
  _padOffset = 0;
  _refreshOffset = 0;
  int nCompressions = (int)(log(_segSize) / log(_K) + 0.5);
  int nCoinShares = nCompressions + 9;
  int nPadShares = 6;
  int nRefreshShares = _disp.hasDisp() ? _segSize + 1 : 0;
  makeRandShares(nCoinShares, _coinShares); // improve: optimize
  makeRandShares(nPadShares, _padShares);   // TODO: use extension field K
  makeRandShares(nRefreshShares, _refreshShares); // improve: only T holds shares
  // ---- # of random double shares:
  // 1. Compute Mult gates in the circuit
  //    -- in total: numOfMultGates
  // 2. Compress Verification
  //    -- 2 * _K * (nCompressions + 1)
  int numDoubleShares = _segSize + (nCompressions+2)*_K*2;
  _doubleOffset = 0;
  makeRandDoubleShares(numDoubleShares, _doubleShares);

  return true;
}

template <class FieldType>
bool ProtocolParty<FieldType>::runOnline(int segStart, int segEnd) {
  if (segStart == 0) {
    inputPhase();
    if ((!verifySharing()) ||
        (!tag())) {
      return false;
    };
    segMsg(segStart, segEnd, "finished inputPhase");
  }
  
  computationPhase(segStart, segEnd);
  if ((!verificationPhase(segStart, segEnd)) ||
      (!verifySharing()) ||
      (!tag())) {
    return false;
  }
  segMsg(segStart, segEnd, "finished computationPhase");

  if (segEnd == _nGates) {
    // last segment: run outputPhase()
    outputPhase();
    if (!verifyOutput()) {
      return false;
    }
    segMsg(segStart, segEnd, "finished outputPhase");
  }
  return true;
}

template <class FieldType>
void ProtocolParty<FieldType>::computationPhase(int segStart, int segEnd) {
  const vector<int>& layers = _circuit.getLayers();
  auto firstGreater = upper_bound(layers.begin(), layers.end(), segStart);

  int layerStart = segStart;
  for (auto iter = firstGreater;
       iter != layers.end() && *iter < segEnd; iter++) {
    int layerEnd = *iter;
    int nMults = processNonMult(layerStart, layerEnd);
    if (nMults > 0) {
      processMult(layerStart, layerEnd, nMults);
    }
    layerStart = layerEnd;
  }
  // the rest of segment is not a full layer
  int layerEnd = segEnd;
  int nMults = processNonMult(layerStart, layerEnd);
  if (nMults > 0) {
    processMult(layerStart, layerEnd, nMults);
  }
  return;
}

template <class FieldType>
void ProtocolParty<FieldType>::
getRefreshTranscripts(int segStart, int segEnd, int king,
                      rTranscript<FieldType>& partyTranscript,
                      vector<vector<FieldType>>& kingTranscript) {
  FieldType lambda = challenge();
  FieldType lambdaI = lambda;
  int multCount = 0;
  // -- collect receiver transcript: x, xp, r, e, o
  for (int k = segStart; k < segEnd; k++) {
    auto& gate = _circuit.getGates()[k];
    if (gate.gateType == MULT) {
      Share<FieldType> xShare = _gateShares[gate.input1];
      xShare.multByConst(lambdaI);
      Share<FieldType> rShare = _refreshShares[multCount];
      rShare.multByConst(lambdaI);
      FieldType oShare = _protoTape.getRecOShareVal(multCount);
      oShare *= lambdaI;
      partyTranscript.addTo(xShare, rShare, oShare);
      multCount++;
      lambdaI *= lambda;
    }
  }
  
  // -- pad with a random transcript
  vector<Share<FieldType>> pad(1);
  getPadShares(1, pad);
  Share<FieldType> xShare = pad[0]; // record before refresh
  vector<FieldType> tmpForRefresh(1, pad[0]._value);
  refresh(tmpForRefresh);
  Share<FieldType> rShare = _refreshShares[multCount];
  FieldType oShare = _protoTape.getRecOShareVal(multCount);
  partyTranscript.addTo(xShare, rShare, oShare);

  if (_myId == king) {
    lambdaI = lambda;
    // -- collect king transcript
    kingTranscript.assign(2, vector<FieldType>(_N, _zero)); // e, o
    for (int i=0; i<multCount; i++) {
      vector<FieldType> eFullShare, oFullShare;
      _protoTape.getRecO2FullShare(i, eFullShare);
      _protoTape.getDealtOShare(i, oFullShare);
      for (int j=0; j<_N; j++) {
        kingTranscript[0][j] += eFullShare[j] * lambdaI;
        kingTranscript[1][j] += oFullShare[j] * lambdaI;
      }
      lambdaI *= lambda;
    }

    // -- pad with the random transcript
    vector<FieldType> eFullShare, oFullShare;
    _protoTape.getRecO2FullShare(multCount, eFullShare);
    _protoTape.getDealtOShare(multCount, oFullShare);
    for (int j=0; j<_N; j++) {
      kingTranscript[0][j] += eFullShare[j];
      kingTranscript[1][j] += oFullShare[j];
    }
  }
}

template <class FieldType>
int ProtocolParty<FieldType>::
localize(vector<FieldType>& zFullShare, Share<FieldType>& zShare,
         vector<FieldType>& dealerSubFullShare, // out
         vector<FieldType>& mySubFullShare) { // out
  // return suspect dealer ID or -1 
  vector<int> NonTSet, TSet;
  int king = _disp.tAndNonTSet(TSet, NonTSet);
  padZShare(zShare, zFullShare);

  // -- decompose and send to king
  vector<FieldType> myShareDecomp;
  _protoTape.decompose(zShare, myShareDecomp, mySubFullShare);
  vector<byte> sendBuf(_N * _fieldByteSize);
  vector<vector<byte>> recBufs(_N, sendBuf), relayBufs;
  _field->elementVectorToByteVector(myShareDecomp, sendBuf);
  _comm.allToOneStore(sendBuf, relayBufs, recBufs, king); // store relay
  vector<vector<FieldType>> partySubFullShares;
  vector<byte> kingMsg, recMsgBuf(_dispMsgSize);
  if (_myId == king) {
    checkPartyDecompShares(kingMsg, partySubFullShares,
                           recBufs, zFullShare);
  }
  _comm.singleBroadcast(king, kingMsg, recMsgBuf);
  if (recMsgBuf[0] == 2) {
    int accused = recMsgBuf[1];
    solveAccuseSum(king, accused, zFullShare[accused], relayBufs[accused]);
    return -1;
  } // else, msg is inconsistent
  // -- verify w/ broadcast the share is indeed inconsistent
  assert(recMsgBuf[0] == 3);
  int dealer = recMsgBuf[1];
  FieldType dealerSubShare = myShareDecomp[dealer];
  dealerSubFullShare.assign(_N, _zero);
  broadcastShare(dealerSubShare, dealerSubFullShare);
  if (!checkTShares(dealerSubFullShare)) {
    return dealer;
  } // else, someone lied before broadcast
  if (_myId == king) {
    for (int partyId=0; partyId<_N; partyId++) {
      if (dealerSubFullShare[partyId] != partySubFullShares[partyId][dealer]) {
        makeAccuseMsg(kingMsg, partyId);
        break;
      }
    }
  }
  _comm.singleBroadcast(king, kingMsg, recMsgBuf);
  assert(recMsgBuf[0] == 2);
  int accused = recMsgBuf[1];
  solveAccuseVal(king, accused, dealer, dealerSubFullShare[accused],
                 relayBufs[accused]);
  return -1;
}

template <class FieldType>
bool ProtocolParty<FieldType>::
pairLocateChunk(int dealer, int accused, // vvv both in and out
                FieldType& accusedBroadcast, Share<FieldType>& share) {
  
  vector<Share<FieldType>> chunkSubShares; // N chunks
  _protoTape.getChunkSubShares(share, chunkSubShares);
  vector<byte> sendBuf, recBuf(_N * _fieldByteSize);
  if (_myId == accused) {
    sendBuf.resize(_N * _fieldByteSize);
    for (int chunkIdx = 0; chunkIdx < _N; chunkIdx++) {
      _field->elementToBytes(sendBuf.data() + chunkIdx * _fieldByteSize,
                             chunkSubShares[chunkIdx]._value);
    }
  }
  _comm.singleBroadcast(accused, sendBuf, recBuf);
  vector<FieldType> accusedChunkVals(_N);
  bytesToVector(recBuf, accusedChunkVals);
  // ---- check the chunks sum up to accused's bc value
  FieldType sum = _zero;
  for (int chunkIdx=0; chunkIdx<_N; chunkIdx++) {
    sum += accusedChunkVals[chunkIdx];
  }
  if (sum != accusedBroadcast) {
    _disp.addCorrParty(accused);
    return true;
  }
  // ---- dealer points out a wrong chunk
  vector<byte> dealerMsg, recMsgBuf(_dispMsgSize);
  if (_myId == dealer) {
    for (int chunkIdx = 0; chunkIdx < _N; chunkIdx++) {
      if (accusedChunkVals[chunkIdx] !=
          _protoTape.getDealtShareValFor(accused, chunkSubShares[chunkIdx])) {
        makeBadChunkMsg(dealerMsg, chunkIdx);
        break;
      }
    }
  }
  _comm.singleBroadcast(dealer, dealerMsg, recMsgBuf);
  assert(recMsgBuf[0] == 4);
  int chunkIdx = recMsgBuf[1];
  share = chunkSubShares[chunkIdx];
  accusedBroadcast = accusedChunkVals[chunkIdx];
  return false;
}

template <class FieldType>
bool ProtocolParty<FieldType>::
pairLocateBatch(int dealer, int accused,
                FieldType& accusedBroadcast, Share<FieldType>& share) {
  vector<Share<FieldType>> batchSubShares;
  int nBatch = _protoTape.getBatchSubShares(share, batchSubShares);
  vector<byte> sendBuf, recBuf(nBatch * _fieldByteSize);
  if (_myId == accused) {
    sendBuf.resize(nBatch * _fieldByteSize);
    for (int batchIdx = 0; batchIdx < nBatch; batchIdx++) {
      _field->elementToBytes(sendBuf.data() + batchIdx * _fieldByteSize,
                             batchSubShares[batchIdx]._value);
    }
  }
  _comm.singleBroadcast(accused, sendBuf, recBuf);
  vector<FieldType> accusedBatchVals(nBatch);
  bytesToVector(recBuf, accusedBatchVals);
  // check the batches sum up to accused's bc value
  FieldType sum = _zero;
  for (int batchIdx=0; batchIdx<nBatch; batchIdx++) {
    sum += accusedBatchVals[batchIdx];
  }
  if (sum != accusedBroadcast) {
    _disp.addCorrParty(accused);
    return true;
  }
  // ---- dealer points out a wrong batch
  vector<byte> dealerMsg, recMsgBuf(_dispMsgSize);
  if (_myId == dealer) {
    for (int batchIdx = 0; batchIdx < _N; batchIdx++) {
      if (accusedBatchVals[batchIdx] !=
          _protoTape.getDealtShareValFor(accused, batchSubShares[batchIdx])) {
        makeBadChunkMsg(dealerMsg, batchIdx);
        break;
      }
    }
  }
  _comm.singleBroadcast(dealer, dealerMsg, recMsgBuf);
  assert(recMsgBuf[0] == 4);
  int batchIdx = recMsgBuf[1];
  share = batchSubShares[batchIdx];
  accusedBroadcast = accusedBatchVals[batchIdx];
  return false;
}


template <class FieldType>
void ProtocolParty<FieldType>::
activeDealer(int dealerId, vector<FieldType>& dealerSubFullShare,
             vector<FieldType>& mySubFullShare, Share<FieldType>& zShare) {
  // -- dealer accuse bad party
  vector<byte> dealerMsg, recMsgBuf(_dispMsgSize);
  if (_myId == dealerId) {
    for (int partyId=0; partyId<_N; partyId++) {
      if (mySubFullShare[partyId] != dealerSubFullShare[partyId]) {
        makeAccuseMsg(dealerMsg, partyId);
        break;
      }
    }
  }
  _comm.singleBroadcast(dealerId, dealerMsg, recMsgBuf);
  assert(recMsgBuf[0] == 2);
  int accused = recMsgBuf[1];
  if (!_disp.isDisp(dealerId, accused)) {
    _disp.addDispPairs(dealerId, accused);
    return;
  }
  Share<FieldType> dealerSubShare; // Z(i) w/ its linear comb info
  _protoTape.getDealerSubShare(dealerId, zShare, dealerSubShare);
  FieldType accusedVal = dealerSubFullShare[accused];
  //  -- round 1: narrow the size by N
  if (pairLocateChunk(dealerId, accused, accusedVal, dealerSubShare)) {
    return;
  }
  // -- round 2: narrow the size by N
  if (pairLocateChunk(dealerId, accused, accusedVal, dealerSubShare)) {
    return;
  }
  // -- round 3: get a single batch and its value broadcasted by accused
  if (pairLocateBatch(dealerId, accused, accusedVal, dealerSubShare)) {
    return;
  }
  // -- accused send batch and tag to all verifiers
  vector<unsigned int> batchRecIdxs;
  vector<FieldType> batchCoeffs;
  int batchIdx = _protoTape.getBatchInfo(dealerSubShare,
                                         batchRecIdxs, batchCoeffs);
  int batchSize = batchRecIdxs.size();

  vector<vector<byte>> sendBufs;
  if (_myId == accused) {
    makeBatchSendBufs(sendBufs, batchRecIdxs, batchIdx);
  }
  // -- check batch and broadcast accept or reject
  vector<byte> recBuf((batchSize+1) * _fieldByteSize);
  _comm.oneToAll(recBuf, sendBufs, accused, true); // use relay (?)
  vector<byte> verifierMsg;
  verifyBatch(verifierMsg, recBuf, accused, batchIdx, accusedVal,
              batchRecIdxs, batchCoeffs);
  vector<vector<byte>> recMsgBufs(_N, vector<byte>(_dispMsgSize));
  _comm.allBroadcast(verifierMsg, recMsgBufs);
  // -- count for t+1 accept
  int acceptCount = 0;
  for (int verifierId=0; verifierId<_N; verifierId++) {
    if (recMsgBufs[verifierId][0] == 0) {
      acceptCount++;
    }
  }
  if (acceptCount > _T) {
    _disp.addCorrParty(dealerId);
  } else {
    _disp.addCorrParty(accused);
  }
}

template <class FieldType>
bool ProtocolParty<FieldType>::
allLocateBatch(vector<FieldType> &prevBroadcasts, Share<FieldType> &share) {
  vector<Share<FieldType>> batchSubShares;
  int nBatch = _protoTape.getBatchSubShares(share, batchSubShares);
  
  // -- all party broadcast N chunk shares
  vector<byte> sendBuf(nBatch * _fieldByteSize);
  vector<vector<byte>> recBufs(_N, sendBuf);
  for (int batchIdx = 0; batchIdx < nBatch; batchIdx++) {
    _field->elementToBytes(sendBuf.data() + batchIdx * _fieldByteSize,
                           batchSubShares[batchIdx]._value);
  }
  _comm.allBroadcast(sendBuf, recBufs);
  // -- check the chunks sum up to previous bc values
  vector<vector<FieldType>> batchFullShares(nBatch);
  for (int partyId=0; partyId < _N; partyId++) {
    bytesToVector(recBufs[partyId], batchFullShares[partyId]);
    FieldType sum = _zero;
    for (int batchIdx = 0; batchIdx < nBatch; batchIdx++) {
      sum += batchFullShares[partyId][batchIdx];
    }
    if (sum != prevBroadcasts[partyId]) {
      _disp.addCorrParty(partyId);
      return true;
    }
  }
  // -- locate a wrong (i.e. inconsistent) chunk
  int badBatchIdx = -1;
  vector<FieldType> batchFullShare(_N);
  for (int batchIdx=0; batchIdx<nBatch; batchIdx++) {
    for (int partyId=0; partyId<_N; partyId++) {
      batchFullShare[partyId] = batchFullShares[partyId][batchIdx];
    }
    if (!checkTShares(batchFullShare)) {
      badBatchIdx = batchIdx;
      break;
    }
  }
  assert(badBatchIdx > 0);
  share = batchSubShares[badBatchIdx];
  prevBroadcasts = batchFullShare;
  return false;
}

template <class FieldType>
bool ProtocolParty<FieldType>::
allLocateChunk(vector<FieldType> &prevBroadcasts, Share<FieldType> &share) {
  vector<Share<FieldType>> chunkSubShares; // N chunks
  _protoTape.getChunkSubShares(share, chunkSubShares);
  // -- all party broadcast N chunk shares
  vector<byte> sendBuf(_N * _fieldByteSize);
  vector<vector<byte>> recBufs(_N, sendBuf);
  for (int chunkIdx = 0; chunkIdx < _N; chunkIdx++) {
    _field->elementToBytes(sendBuf.data() + chunkIdx * _fieldByteSize,
                           chunkSubShares[chunkIdx]._value);
  }
  _comm.allBroadcast(sendBuf, recBufs);
  // -- check the chunks sum up to previous bc values
  vector<vector<FieldType>> chunkFullShares(_N);
  for (int partyId=0; partyId < _N; partyId++) {
    bytesToVector(recBufs[partyId], chunkFullShares[partyId]);
    FieldType sum = _zero;
    for (int chunkIdx = 0; chunkIdx < _N; chunkIdx++) {
      sum += chunkFullShares[partyId][chunkIdx];
    }
    if (sum != prevBroadcasts[partyId]) {
      _disp.addCorrParty(partyId);
      return true;
    }
  }
  // -- locate a wrong (i.e. inconsistent) chunk
  int badChunkIdx = -1;
  vector<FieldType> chunkFullShare(_N);
  for (int chunkIdx=0; chunkIdx<_N; chunkIdx++) {
    for (int partyId=0; partyId<_N; partyId++) {
      chunkFullShare[partyId] = chunkFullShares[partyId][chunkIdx];
    }
    if (!checkTShares(chunkFullShare)) {
      badChunkIdx = chunkIdx;
      break;
    }
  }
  assert(badChunkIdx > 0);
  share = chunkSubShares[badChunkIdx];
  prevBroadcasts = chunkFullShare;
  return false;
}

template <class FieldType>
void ProtocolParty<FieldType>::
padZShare(Share<FieldType> &zShare,
          vector<FieldType> &zFullShare) {
  vector<Share<FieldType>>  randPad(1);
  getPadShares(1, randPad);
  vector<FieldType> rFullShare;
  broadcastShare(randPad[0]._value, rFullShare);
  if (!checkTShares(rFullShare)) {
    zShare = randPad[0];
    zFullShare = rFullShare;
  } else {
    zShare += randPad[0];
    for (int i=0; i<_N; i++) {
      zFullShare[i] += rFullShare[i];
    }
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
solveAccuseSum(int king, int accused, FieldType broadcastVal, vector<byte> &relayBuf) {
  if (_disp.isDisp(king, accused)) {
    _disp.addDispPairs(king, accused);
    return;
  }
  int relayer = _disp.relayer(king, accused);
  vector<byte> relayerMsg, recMsgBuf(_dispMsgSize);
  if (_myId == relayer) {
    vector<FieldType> partyShares;
    bytesToVector(relayBuf, partyShares);
    assert(partyShares.size() == _N);
    FieldType sum = _zero;
    for (int i=0; i<_N; i++) {
      sum += partyShares[i];
    }
    if (sum != broadcastVal) {
      makeOkMsg(relayerMsg);    // agrees w/ king
    } else {
      makeRejectMsg(relayerMsg); // disagrees w/ king
    }
  }
  _comm.singleBroadcast(relayer, relayerMsg, recMsgBuf);
  if (recMsgBuf[0] == 0) {      // relayer agrees w/ king
    _disp.addDispPairs(relayer, accused);
  } else {                      // relayer disagrees w/ king
    _disp.addDispPairs(relayer, king);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
solveAccuseVal(int king, int accused, int dealer,
               FieldType broadcastVal, vector<byte>& relayBuf) {
  if (_disp.isDisp(king, accused)) {
    _disp.addDispPairs(king, accused);
    return;
  }
  int relayer = _disp.relayer(king, accused);
  vector<byte> relayerMsg, recMsgBuf(_dispMsgSize);
  if (_myId == relayer) {
    vector<FieldType> partyShares;
    bytesToVector(relayBuf, partyShares);
    assert(partyShares.size() == _N);
    if (partyShares[dealer] != broadcastVal) {
      makeOkMsg(relayerMsg);
    } else {
      makeRejectMsg(relayerMsg);
    }
  }
  _comm.singleBroadcast(relayer, relayerMsg, recMsgBuf);
  if (recMsgBuf[0] == 0) {      // relayer agrees w/ king
    _disp.addDispPairs(relayer, accused);
  } else {                      // relayer disagrees w/ king
    _disp.addDispPairs(relayer, king);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
broadcastShare(FieldType share, vector<FieldType>& fullShare) {
  vector<byte> sendBuf(_fieldByteSize);
  vector<vector<byte>> recBufs(_N, sendBuf);
  _field->elementToBytes(sendBuf.data(), share);
  _comm.allBroadcast(sendBuf, recBufs);
  fullShare.resize(_N);
  for (int partyId=0; partyId<_N; partyId++) {
    fullShare[partyId] = _field->bytesToElement(recBufs[partyId].data());
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
corruptedDealer(int dealerId, vector<FieldType>& dealerSubFullShare,
                Share<FieldType>& zShare) {
  Share<FieldType> dealerSubShare; // Z(i) w/ its linear comb info
  _protoTape.getDealerSubShare(dealerId, zShare, dealerSubShare);
  // -- round 1:
  if (allLocateChunk(dealerSubFullShare, dealerSubShare)) {
    return;
  }
  // -- round 2:
  if (allLocateChunk(dealerSubFullShare, dealerSubShare)) {
    return;
  }
  // -- round 3:
  if (allLocateBatch(dealerSubFullShare, dealerSubShare)) {
    return;
  }
  vector<unsigned int> batchRecIdxs;
  vector<FieldType> batchCoeffs;
  int batchIdx = _protoTape.getBatchInfo(dealerSubShare,
                                         batchRecIdxs, batchCoeffs);
  int batchSize = batchRecIdxs.size();
  // -- each party broadcast batch
  vector<vector<byte>> sendBufs, recBufs;
  makeBatchSendBufs(sendBufs, batchRecIdxs, batchIdx);
  recBufs.assign(_N, vector<byte>((batchSize+1)*_fieldByteSize));
  _comm.allToAll(sendBufs, recBufs, true); // use relay
  // -- check each party's batch
  vector<vector<byte>> msgs(_N);
  for (int partyId=0; partyId<_N; partyId++) {
    verifyBatch(msgs[partyId], recBufs[partyId], partyId, batchIdx,
                dealerSubFullShare[partyId], batchRecIdxs, batchCoeffs);
    vector<vector<byte>> recMsgBufs(_N, vector<byte>(_dispMsgSize));
    _comm.allBroadcast(msgs[partyId], recMsgBufs);
    int acceptCount = 0;
    for (int verifierId=0; verifierId<_N; verifierId++) {
      if (recMsgBufs[verifierId][0] == 0) {
        acceptCount++;
      }
    }
    if (acceptCount <= _T) {
      _disp.addCorrParty(partyId);
    }
  }
}


template <class FieldType>
void ProtocolParty<FieldType>::
analyzeSharing(vector<FieldType>& zFullShare, Share<FieldType>& zShare) {
  vector<FieldType> dealerSubFullShare, mySubFullShare;
  int dealerId = localize(zFullShare, zShare, dealerSubFullShare, mySubFullShare);
  if (dealerId < 0) {
    return;
  }

  if(!_disp.isCorrupt(dealerId)) {
    activeDealer(dealerId, dealerSubFullShare, mySubFullShare, zShare);
  }else {
    corruptedDealer(dealerId, dealerSubFullShare, zShare);
  }
  return;
}

template <class FieldType>
void ProtocolParty<FieldType>::
solveAccuseRand(int king, int accused, vector<byte>& relayBuf) {
  // TODO: minor refactor all solveAccused..
  // functions to reuse checks as the dealer
  if (_disp.isDisp(king, accused)) {
    _disp.addDispPairs(king, accused);
    return;
  }
  int relayer = _disp.relayer(king, accused);
  vector<byte> relayerMsg, recMsgBuf(_dispMsgSize);
  if (_myId == relayer) {
    makeRejectMsg(relayerMsg);  // default to reject if no error found
    vector<FieldType> fullShare;
    bytesToVector(relayBuf, fullShare);
    assert(fullShare.size() == _N);
    if (!checkTShares(fullShare)) {
      makeOkMsg(relayerMsg);
    } else {
      for (int partyId=0; partyId<_N; partyId++) {
        if (_disp.isDisp(accused, partyId) && fullShare[partyId] != 0) {
          makeOkMsg(relayerMsg);
          break;
        }
      }
    }
  }
  _comm.singleBroadcast(relayer, relayerMsg, recMsgBuf);
  if (recMsgBuf[0] == 0) {      // relayer agrees w/ king
    _disp.addDispPairs(relayer, accused);
  } else {                      // relayer disagress w/ king
    _disp.addDispPairs(relayer, king);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
solveAccuseRand2T(int king, int accused, FieldType secret,
                  vector<byte> &relayBuf) {
  if (_disp.isDisp(king, accused)) {
    _disp.addDispPairs(king, accused);
    return;
  }
  int relayer = _disp.relayer(king, accused);
  vector<byte> relayerMsg, recMsgBuf(_dispMsgSize);
  if (_myId == relayer) {
    makeRejectMsg(relayerMsg);  // default to reject if no error found
    vector<FieldType> fullShare;
    bytesToVector(relayBuf, fullShare);
    assert(fullShare.size() == _N);

    if (secret != getSecret(fullShare)) {
      makeOkMsg(relayerMsg);
    } else {
      for (int partyId=0; partyId<_N; partyId++) {
        if (_disp.isDisp(accused, partyId) && fullShare[partyId] != 0) {
          makeOkMsg(relayerMsg);
          break;
        }
      }
    }
  }
  _comm.singleBroadcast(relayer, relayerMsg, recMsgBuf);
  if (recMsgBuf[0] == 0) {      // relayer agrees w/ king
    _disp.addDispPairs(relayer, accused);
  } else {                      // relayer disagress w/ king
    _disp.addDispPairs(relayer, king);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
solveOpenRand(int king, vector<byte> &openMsg,
              vector<vector<byte>>& relayBufs, vector<FieldType> &myShareDecomp,
              vector<FieldType> &myFullShare) {
  int dealer, tParty;
  FieldType dealerToTShare, tFromDealerShare;
  readOpenMsg(openMsg, dealer, tParty, dealerToTShare, tFromDealerShare);
  // -- king vs tParty
  vector<byte> sendBuf(_fieldByteSize), recBuf(_fieldByteSize);
  if (_myId == tParty) {
    _field->elementToBytes(sendBuf.data(), myShareDecomp[dealer]);
  }
  _comm.singleBroadcast(tParty, sendBuf, recBuf);
  if (tFromDealerShare != _field->bytesToElement(recBuf.data())) {
    _disp.addDispPairs(king, tParty);
    return;
  } // else: tParty agrees w/ dealer, check dealer next
  // -- king vs dealer
  if (_myId == dealer) {
    _field->elementToBytes(sendBuf.data(), myFullShare[tParty]);
  }
  _comm.singleBroadcast(dealer, sendBuf, recBuf);
  FieldType dealerBroadcastShare = _field->bytesToElement(recBuf.data());
  if (dealerToTShare == dealerBroadcastShare) {
    _disp.addDispPairs(dealer, tParty);
    return;
  } // else: dealer disagrees w/ king
  if (!_disp.isDisp(king, dealer)) {
    _disp.addDispPairs(dealer, king);
    return;
  } // else: let relayer decide
  int relayer = _disp.relayer(king, dealer);
  checkRelayer(relayer, king, dealer, dealerToTShare,
               _field->bytesToElement(relayBufs[dealer].data() +
                                      tParty * _fieldByteSize));
}

template <class FieldType>
void ProtocolParty<FieldType>::
checkRand(vector<FieldType>& fullShare, Share<FieldType>& share) {
  vector<FieldType> myShareDecomp, myFullShare;
  _protoTape.decompose(share, myShareDecomp, myFullShare);
  vector<int> NonTSet, TSet;
  int king = _disp.tAndNonTSet(TSet, NonTSet);
  vector<byte> sendBuf(_N * _fieldByteSize);
  vector<vector<byte>> recBufs(_N, sendBuf), relayBufs;
  _field->elementVectorToByteVector(myFullShare, sendBuf);
  _comm.allToOneStore(sendBuf, relayBufs, recBufs, king); // store relay

  // ---- check 1: dealers sent correct shares?
  vector<vector<FieldType>> dealerFullShares; // each row is a full share
  vector<FieldType> secrets;
  vector<byte> kingMsg, recMsgBuf(_dispMsgSize);
  if (_myId == king) {
    checkDealerRandShares(kingMsg, recBufs, secrets, dealerFullShares);
  }
  _comm.singleBroadcast(king, kingMsg, recMsgBuf);
  if (recMsgBuf[0] != 0) {
    int accused = recMsgBuf[1];
    solveAccuseRand(king, accused, relayBufs[accused]);
    return;
  }
  // ---- check 2: dealer and T agree?
  _field->elementVectorToByteVector(myShareDecomp, sendBuf);
  _comm.TToKing(king, sendBuf, recBufs);
  if (_myId == king) {
    checkDealerAndT(kingMsg, TSet, recBufs, dealerFullShares);
  }
  _comm.singleBroadcast(king, kingMsg, recMsgBuf);
  assert(recMsgBuf[0] != 0);
  solveOpenRand(king, recMsgBuf, relayBufs, myShareDecomp, myFullShare);
  return;
}

template <class FieldType>
void ProtocolParty<FieldType>::
checkDoubleRand(vector<FieldType>& fullRandShare,
                vector<FieldType>& fullRand2TShare,
                Share<FieldType>& randShare, Share<FieldType>& rand2TShare) {
  vector<FieldType> myShareDecomp, myFullShare;
  _protoTape.decompose(randShare, myShareDecomp, myFullShare);
  vector<FieldType> my2TShareDecomp, myFull2TShare;
  _protoTape.decompose2T(rand2TShare, my2TShareDecomp, myFull2TShare);
  vector<int> NonTSet, TSet;
  int king = _disp.tAndNonTSet(TSet, NonTSet);
  vector<byte> sendBuf(_N * _fieldByteSize);
  vector<vector<byte>> recBufs(_N, sendBuf), relayBufs, relayBufs2T;
  _field->elementVectorToByteVector(myFullShare, sendBuf);
  _comm.allToOneStore(sendBuf, relayBufs, recBufs, king); // store relay

  // improve: merge check 1 and check 2?
  // ---- check 1: dealer sent correct t-shares?
  vector<vector<FieldType>> dealerFullShares, relayShares; // king and relayer
  vector<FieldType> secrets;
  vector<byte> kingMsg, recMsgBuf(_dispMsgSize);
  if (_myId == king) {
    checkDealerRandShares(kingMsg, recBufs, secrets, dealerFullShares);
  }
  _comm.singleBroadcast(king, kingMsg, recMsgBuf);
  if (recMsgBuf[0] != 0) {
    int accused = recMsgBuf[1];
    solveAccuseRand(king, accused, relayBufs[accused]);
    return;
  }
  // ---- check 2: dealer sent correct 2t-shares?
  vector<vector<FieldType>> dealerFull2TShares, relay2TShares;
  _field->elementVectorToByteVector(myFull2TShare, sendBuf);
  _comm.allToOneStore(sendBuf, relayBufs2T, recBufs, king); // store relay
  if (_myId == king) {
    checkDealerRand2TShares(kingMsg, recBufs, secrets, dealerFull2TShares);
  }
  _comm.singleBroadcast(king, kingMsg, recMsgBuf);
  if (recMsgBuf[0] != 0) {
    int accused = recMsgBuf[1];
    solveAccuseRand2T(king, accused, secrets[accused], relayBufs2T[accused]);
    return;
  }
  // ---- check 3: dealer and T agree?
  _field->elementVectorToByteVector(myShareDecomp, sendBuf);
  _comm.TToKing(king, sendBuf, recBufs);
  if (_myId == king) {
    checkDealerAndT(kingMsg, TSet, recBufs, dealerFullShares);
  }
  _comm.singleBroadcast(king, kingMsg, recMsgBuf);
  if (recMsgBuf[0] != 0) {
    solveOpenRand(king, recMsgBuf, relayBufs, myShareDecomp, myFullShare);
    return;
  } // else: check 2TShares
  _field->elementVectorToByteVector(my2TShareDecomp, sendBuf);
  _comm.TToKing(king, sendBuf, recBufs);
  if (_myId == king) {
    checkDealerAndT(kingMsg, TSet, recBufs, dealerFull2TShares);
  }
  _comm.singleBroadcast(king, kingMsg, recMsgBuf);
  assert(recMsgBuf[0] != 0);
  solveOpenRand(king, recMsgBuf, relayBufs2T, my2TShareDecomp, myFull2TShare);
  return;
}

template <class FieldType>
void ProtocolParty<FieldType>::
broadcastKingTranscript(int king, vector<vector<FieldType> > &kingTranscript) {
  vector<byte> recBuf(_N * _fieldByteSize);
  vector<byte> sendBuf(_N * _fieldByteSize);
  if (_myId == king) {          // improve: combine to 1 round
    _field->elementVectorToByteVector(kingTranscript[0], sendBuf);
    _comm.singleBroadcast(king, sendBuf, recBuf);
    _field->elementVectorToByteVector(kingTranscript[1], sendBuf);
    _comm.singleBroadcast(king, sendBuf, recBuf);
  } else {
    _comm.singleBroadcast(king, sendBuf, recBuf);
    kingTranscript[0].resize(_N);
    for (int i=0; i<_N; i++) {
      kingTranscript[0][i] =
        _field->bytesToElement(recBuf.data() + i * _fieldByteSize);
    }
    _comm.singleBroadcast(king, sendBuf, recBuf);
    for (int i=0; i<_N; i++) {
      kingTranscript[1].resize(_N);
      kingTranscript[1][i] =
        _field->bytesToElement(recBuf.data() + i * _fieldByteSize);
    }
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
broadcastPartyTranscript(rTranscript<FieldType>& myTranscript,
                         vector<vector<FieldType>>& allTranscriptVals) {
  vector<byte> sendBuf(5 * _fieldByteSize);
  vector<vector<byte>> recBufs(_N, vector<byte>(5 * _fieldByteSize));
  vector<FieldType> values = myTranscript.getValues();
  _field->elementVectorToByteVector(values, sendBuf);
  _comm.allBroadcast(sendBuf, recBufs);

  allTranscriptVals.assign(5, vector<FieldType>(_N, _zero));
  for (int curId=0; curId<_N; curId++) {
    if (_disp.isCorrupt(curId)) { continue; }
    for (int j = 0; j < 5; j++) {
      allTranscriptVals[j][curId] = 
        _field->bytesToElement(recBufs[curId].data() + j * _fieldByteSize);
    }
  }
}

template <class FieldType>
bool ProtocolParty<FieldType>::checkRefresh(int segStart, int segEnd) {
  if (!_disp.hasCorrupt()) {
    return true;
  }
  vector<int> NonTSet, TSet;
  int king = _disp.tAndNonTSet(TSet, NonTSet);

  // -- prepare (virtual) transcripts
  rTranscript<FieldType> myTranscript;         // x, xp, r, e, o
  vector<vector<FieldType>> kingTranscript(2); // e, o
  getRefreshTranscripts(segStart, segEnd, king, myTranscript, kingTranscript);

  // -- broadcast and compare
  vector<vector<FieldType>> allTranscriptVals;
  broadcastPartyTranscript(myTranscript, allTranscriptVals);
  // ---- check 1: is [x] consistent? (else: analyzeSharing)
  if (!checkTShares(allTranscriptVals[0])) {
    analyzeSharing(allTranscriptVals[0], myTranscript._x);
    return false;
  }
  // ---- check 2: is [r] consistent? (else: checkRand)
  if (!checkTShares(allTranscriptVals[2])) {
    checkRand(allTranscriptVals[2], myTranscript._r);
    return false;
  }
  // ---- check 3: king agrees with T? (else: disp pair)
  broadcastKingTranscript(king, kingTranscript);
  for (int tid : TSet) {
    if (allTranscriptVals[3][tid] != kingTranscript[0][tid] ||
        allTranscriptVals[4][tid] != kingTranscript[1][tid]) {
      _disp.addDispPairs(king, tid);
      return false;
    }
  }
  // ---- check 4: some party didn't follow the protocol?
  // TODO: minor
  // 1. king's shares have zeros for nonT
  // 2. consistency of king's shares
  // 3. xp share = x share - o share
  return true;
}

template <class FieldType>
bool ProtocolParty<FieldType>::verificationPhase(int segStart, int segEnd) {
  if (!checkRefresh(segStart, segEnd)) {
    return false;
  }
  vector<int> NonTSet, TSet;
  int king = _disp.tAndNonTSet(TSet, NonTSet);
  vector<Share<FieldType>> aShares, bShares;
  mTranscript<FieldType> myTranscript;         // r, r2, e, e2, c
  vector<vector<FieldType>> kingTranscript(2); // e, e2
  vector<FieldType> relayTranscript(_N);       // e2
  delinearization(segStart, segEnd, king, aShares, bShares,
                  myTranscript, kingTranscript, relayTranscript);
  while(dimensionReduction(king, aShares, bShares, myTranscript,
                           kingTranscript, relayTranscript) >= _K);
  randomization(king, aShares, bShares,
                myTranscript, kingTranscript, relayTranscript);
  return (checkSingleMult(king, aShares[0], bShares[0],
                          myTranscript, kingTranscript, relayTranscript));
}

template <class FieldType>
void ProtocolParty<FieldType>::
delinearization(int segStart, int segEnd, int king,
                vector<Share<FieldType>>& aShares, vector<Share<FieldType>>& bShares,
                mTranscript<FieldType>& transcript,
                vector<vector<FieldType>>& kingTranscript,
                vector<FieldType>& relayTranscript){
  aShares.resize(_segSize); bShares.resize(_segSize);
  FieldType lambda = challenge();
  FieldType lambdaI = lambda;
  // -- collect receiver transcript
  int multCount = 0;
  for (int k = segStart; k < segEnd; k++) {
    auto& gate = _circuit.getGates()[k];
    if(gate.gateType == MULT) {
      aShares[multCount] = _gateShares[gate.input1];
      if (_disp.hasCorrupt()) {
        aShares[multCount] =
          aShares[multCount] - _protoTape.getRecOShare(multCount);
      }
      bShares[multCount] = _gateShares[gate.input2];
      bShares[multCount].multByConst(lambdaI);
      Share<FieldType> rShare = _doubleShares[multCount*2];
      rShare.multByConst(lambdaI);
      Share<FieldType> r2Share = _doubleShares[multCount*2+1];
      r2Share.multByConst(lambdaI);
      FieldType eShare = lambdaI * _protoTape.getRecE1ShareVal(multCount);
      transcript.addTo(rShare, r2Share, eShare,
                       aShares[multCount]._value * bShares[multCount]._value);
      multCount++;
      lambdaI *= lambda;
    }
  }
  aShares.resize(multCount);
  bShares.resize(multCount);

  if (_myId == king) {
    lambdaI = lambda;
    // -- collect king transcript
    kingTranscript.assign(2, vector<FieldType>(_N, _zero)); // e e2
    for (int i=0; i<multCount; i++) {
      vector<FieldType> e1FullShare, e2FullShare;
      _protoTape.getRecE2FullShare(i, e2FullShare);
      _protoTape.getDealtE1Share(i, e1FullShare);
      for (int j=0; j<_N; j++) {
        kingTranscript[0][j] += e1FullShare[j] * lambdaI;
        kingTranscript[1][j] += e2FullShare[j] * lambdaI;
      }
      lambdaI *= lambda;
    }
  }

  if (_disp.isRelayer(king, _myId)) {
    lambdaI = lambda;
    // -- collect relayer transcript
    relayTranscript.assign(_N, _zero);
    for (int i=0; i<multCount; i++) {
      vector<FieldType> relayShares(_N);
      _protoTape.getE2RelayShare(i, relayShares);
      for (int j=0; j<_N; j++) {
        relayTranscript[j] += relayShares[j] * lambdaI;
      }
      lambdaI *= lambda;
    }
  }
}

template <class FieldType>
int ProtocolParty<FieldType>::
dimensionReduction(int king, vector<Share<FieldType>>& aShares,
                   vector<Share<FieldType>>& bShares,
                   mTranscript<FieldType>& transcript,
                   vector<vector<FieldType>>& kingTranscript,
                   vector<FieldType>& relayTranscript){
  int totalLength = aShares.size();
  if (totalLength < _K) { // no need to compress
    return totalLength;
  }
  // -- divide into K groups and pad short groups
  int groupSize = (totalLength + _K -1 )/ _K;
  totalLength = groupSize * _K;
  aShares.resize(totalLength);
  bShares.resize(totalLength);

  // -- one DN mult for each group i to build dShares 0 .. k-2
  vector<mTranscript<FieldType>> transcripts(_K);       // n
  vector<vector<vector<FieldType>>> kingTranscripts(2); // 2 * n * _N
  vector<vector<FieldType>> relayTranscripts;           // n * _N
  dotProd(king, groupSize, aShares, bShares,
          transcripts, kingTranscripts, relayTranscripts);
  
  // build dShares k-1: c - d_0 - ... - d_{k-2}
  transcripts[_K-1] = transcript;
  for (int i=0; i<_K-1; i++) {
    transcripts[_K-1] -= transcripts[i];
  }
  if (_myId == king) {
    kingTranscripts[0][_K-1] = kingTranscript[0];
    kingTranscripts[1][_K-1] = kingTranscript[1];
    for (int i=0; i<_K-1; i++) {
      for (int j=0; j<_N; j++) {
        kingTranscripts[0][_K-1][j] =
          kingTranscripts[0][_K-1][j] - kingTranscripts[0][i][j];
        kingTranscripts[1][_K-1][j] =
          kingTranscripts[1][_K-1][j] - kingTranscripts[1][i][j];
      }
    }
  }
  if (_disp.isRelayer(king, _myId)) {
    relayTranscripts[_K-1] = relayTranscript;
    for (int i=0; i<_K-1; i++) {
      for (int j=0; j<_N; j++) {
        relayTranscripts[_K-1][j] =
          relayTranscripts[_K-1][j] - relayTranscripts[i][j];
      }
    }
  }

  compress(king, groupSize, aShares, bShares,
           transcripts, kingTranscripts, relayTranscripts,
           transcript, kingTranscript, relayTranscript);
  assert(aShares.size() == groupSize);
  assert(bShares.size() == groupSize);
  return groupSize;
}

template <class FieldType>
void ProtocolParty<FieldType>::
compress(int king, int groupSize,
         vector<Share<FieldType>>& aShares, vector<Share<FieldType>>& bShares,
         vector<mTranscript<FieldType>>& transcripts,
         vector<vector<vector<FieldType>>>& kingTranscripts,
         vector<vector<FieldType>>& relayTranscripts,
         mTranscript<FieldType>& transcriptNew,
         vector<vector<FieldType>>& kingTranscriptNew,
         vector<FieldType>& relayTranscriptNew) {
  int totalLength = groupSize * _K;

  // -- compress k groups of vector product  
  buildPolyVecInd(king, groupSize, aShares, bShares, transcripts, kingTranscripts, relayTranscripts);
  
  FieldType lambda = challenge();
  vector<FieldType> beta_lambda(1);
  beta_lambda[0] = lambda;  
  HIMp<FieldType> matrix_for_k_lambda(1, _K, _field);
  matrix_for_k_lambda.InitHIMByVectors(_alpha_k, beta_lambda);
  HIMp<FieldType> matrix_for_2k_lambda(1, 2*_K-1, _field);
  matrix_for_2k_lambda.InitHIMByVectors(_alpha_2km1, beta_lambda);

  // improve: refactor
  vector<Share<FieldType>> ySharesA, ySharesB, buffer(1),
    aSharesNew(groupSize), bSharesNew(groupSize);
  for (int i=0; i<groupSize; i++) {
    ySharesA.clear();
    ySharesB.clear();
    for (int j=i; j<totalLength; j+=groupSize) {
      ySharesA.push_back(aShares[j]);
      ySharesB.push_back(bShares[j]);
    }      
    matrix_for_k_lambda.MatrixMultShares(ySharesA, buffer);
    aSharesNew[i] = buffer[0];
    matrix_for_k_lambda.MatrixMultShares(ySharesB, buffer);
    bSharesNew[i] = buffer[0];
  }
  aShares = aSharesNew;
  bShares = bSharesNew;

  interpTranscripts(matrix_for_2k_lambda, transcripts, transcriptNew);
  if (_myId == king) {
    interpKingTranscripts(matrix_for_2k_lambda,
                          kingTranscripts, kingTranscriptNew);
  }
  if (_disp.isRelayer(king, _myId)) {
    interpRelayTranscripts(matrix_for_2k_lambda,
                           relayTranscripts, relayTranscriptNew);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
randomization(int king, vector<Share<FieldType>>& aShares,
              vector<Share<FieldType>>& bShares,
              mTranscript<FieldType>& transcript,
              vector<vector<FieldType>>& kingTranscript,
              vector<FieldType>& relayTranscript) {
  assert(aShares.size() < _K);
  assert(bShares.size() < _K);
  aShares.resize(_K);
  bShares.resize(_K);

  vector<Share<FieldType>>  randPair(2);
  getPadShares(2, randPair);
  aShares[_K-1] = randPair[0];
  bShares[_K-1] = randPair[1];

  vector<mTranscript<FieldType>> transcripts;
  vector<vector<vector<FieldType>>> kingTranscripts(2);
  vector<vector<FieldType>> relayTranscripts;
  dotProd(king, 1, aShares, bShares,
          transcripts, kingTranscripts, relayTranscripts);

  // build dShares k-1: c - d_0 - ... - d_{k-2}
  transcripts[_K-1] += transcript;
  for (int i = 0; i< _K - 1; i++) {
    transcripts[_K-1] -= transcripts[i];
  }
  if (_myId == king) {
    for (int j = 0; j < _N; j++) {
      kingTranscripts[0][_K - 1][j] += kingTranscript[0][j];
      kingTranscripts[1][_K - 1][j] += kingTranscript[1][j];
    }
    for (int i=0; i<_K-1; i++) {
      for (int j=0; j < _N; j++) {
        kingTranscripts[0][_K-1][j] =
          kingTranscripts[0][_K-1][j] - kingTranscripts[0][i][j];
        kingTranscripts[1][_K-1][j] =
          kingTranscripts[1][_K-1][j] - kingTranscripts[1][i][j];
      }
    }
  }
  if (_disp.isRelayer(king, _myId)) {
    for(int j=0; j<_N; j++) {
      relayTranscripts[_K-1][j] += relayTranscript[j];
    }
    for (int i=0; i<_K-1; i++) {
      for (int j = 0; j < _N; j++) {
        relayTranscripts[_K-1][j] =
          relayTranscripts[_K-1][j] - relayTranscripts[i][j];
      }
    }
  }

  compress(king, 1, aShares, bShares,
           transcripts, kingTranscripts, relayTranscripts,
           transcript, kingTranscript, relayTranscript);
  assert(aShares.size() == 1);
  assert(bShares.size() == 1);
  return;
}

template <class FieldType>
bool ProtocolParty<FieldType>::
checkSingleMult(int king, Share<FieldType>& aShare, Share<FieldType>& bShare,
                mTranscript<FieldType>& transcript,
                vector<vector<FieldType>>& kingTranscript,
                vector<FieldType>& relayTranscript) {
  // -- broadcast and compare
  vector<byte> sendBuf(7 * _fieldByteSize);
  vector<vector<byte>> recBufs(_N, vector<byte>(7 * _fieldByteSize));
  vector<FieldType> sendVec {aShare._value, bShare._value};
  vector<FieldType> tmpVals = transcript.getValues();
  sendVec.insert(sendVec.end(), tmpVals.begin(), tmpVals.end());
  _field->elementVectorToByteVector(sendVec, sendBuf);
  _comm.allBroadcast(sendBuf, recBufs);
  vector<vector<FieldType>> partyTranscripts(5, vector<FieldType>(_N));
  vector<FieldType> aFullShare(_N, _zero), bFullShare(_N, _zero);
  for (int curId = 0; curId < _N; curId++) {
    if (_disp.isCorrupt(curId)) { continue; }
    aFullShare[curId] = _field->bytesToElement(recBufs[curId].data());
    bFullShare[curId] = _field->bytesToElement(recBufs[curId].data() +
                                               _fieldByteSize);
    for (int j=0; j<5; j++) {
      partyTranscripts[j][curId] = 
        _field->bytesToElement(recBufs[curId].data() + (j+2) * _fieldByteSize);
    }
  }
  // ---- check 1: are [a] and [b] consistent? (else: analyzeSharing)
  if (!checkTShares(aFullShare)) {
    analyzeSharing(aFullShare, aShare);
    return false;
  }
  if (!checkTShares(bFullShare)) {
    analyzeSharing(bFullShare, bShare);
    return false;
  }

  // ---- check 2: is [r] consistent, and is r2 the same as r?
  if (!checkTShares(partyTranscripts[0]) ||
      getSecret(partyTranscripts[0]) != getSecret(partyTranscripts[1])) {
    checkDoubleRand(partyTranscripts[0], partyTranscripts[1],
                    transcript._r, transcript._r2);
    return false;
  }

  // ---- check 3: king vs parties
  broadcastKingTranscript(king, kingTranscript);
  for (int partyId =0; partyId<_N; partyId++) {
    // check e1: king vs. T
    if (partyTranscripts[2][partyId] != kingTranscript[0][partyId]) {
      _disp.addDispPairs(partyId, king);
      return false;
    }
    // check e1: king vs. all
    if (partyTranscripts[3][partyId] != kingTranscript[1][partyId]) {
      if (_disp.isDisp(partyId, king)) {
        _disp.addDispPairs(partyId, king);
      } else {
        vector<byte> relayerMsg, recMsgBuf(_dispMsgSize);
        int relayer = _disp.relayer(king, partyId);
        if (_myId == relayer) {
          if (relayTranscript[partyId] == kingTranscript[1][partyId]) {
            makeOkMsg(relayerMsg);
          } else {
            makeRejectMsg(relayerMsg);
          }
        }
        _comm.singleBroadcast(relayer, relayerMsg, recMsgBuf);
        if (recMsgBuf[0] == 0) {
          _disp.addDispPairs(relayer, partyId); // agrees w/ king
        } else {
          _disp.addDispPairs(relayer, king); // disagrees w/ king
        }
      }
      return false;
    }
  }
  // ---- check 4: some party didn't follow the protocol?
  // TODO: minor
  // 1. king's shares not consistent
  // 2. cShare != e1Share + r1Share
  // 3. e2Share != a*bShare - r2Share
  
  return true;
}

template <class FieldType>
void ProtocolParty<FieldType>::
buildPolyVecIndWorker(vector<Share<FieldType>>& aShares,
                      vector<Share<FieldType>>& bShares,
                      vector<vector<Share<FieldType>>>& ASharesEval,
                      vector<vector<Share<FieldType>>>& BSharesEval,
                      int groupSize, int threadId) {
  int totalLength = aShares.size();
  vector<Share<FieldType>> ySharesA, ySharesB;   // tmp vectors
  
  for (int i=threadId; i<groupSize; i+= _numThreads) {
    ySharesA.clear();
    ySharesB.clear();
    for (int j=i; j<totalLength; j+=groupSize) {
      ySharesA.push_back(aShares[j]);
      ySharesB.push_back(bShares[j]);
    }
    _mat_check_km1.MatrixMultShares(ySharesA, ASharesEval[i]);
    _mat_check_km1.MatrixMultShares(ySharesB, BSharesEval[i]);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
buildPolyVecInd(int king, int groupSize, vector<Share<FieldType>>& aShares,
                vector<Share<FieldType>>& bShares,
                vector<mTranscript<FieldType>>& transcripts,
                vector<vector<vector<FieldType>>>& kingTranscripts,
                vector<vector<FieldType>>& relayTranscripts) {

  // batchDegree should be _K
  int batchDegree = (aShares.size() + groupSize - 1) / groupSize;
  
  // interpolate A[i], B[i] and evaluate at A(k), ..., A(2k-2), B(k), ..., B(2k-2)
  vector<vector<Share<FieldType>>>
    ASharesEval(groupSize, vector<Share<FieldType>>(batchDegree-1)),
    BSharesEval(groupSize, vector<Share<FieldType>>(batchDegree-1));
  vector<thread> threads(_numThreads-1);
  // bring up thread workers  
  for (int t=0; t < _numThreads-1; t++) {
    threads[t] =
      thread(&ProtocolParty::buildPolyVecIndWorker, this, ref(aShares),
             ref(bShares), ref(ASharesEval), ref(BSharesEval), groupSize, t+1);
  }
  buildPolyVecIndWorker(aShares, bShares,
                        ASharesEval, BSharesEval, groupSize, 0);
  // join thread workers
  for (int t=0; t < _numThreads-1; t++) {
    threads[t].join();
  }

  // vector-mult to get the rest of dShares
  int nEvals = groupSize * (batchDegree - 1);
  vector<Share<FieldType>> ASharesEvalFlat(nEvals), BSharesEvalFlat(nEvals);
  for (int i=0; i<batchDegree-1; i++) {
    // point i
    for (int j=0; j<groupSize; j++) {
      // poly j
      ASharesEvalFlat[i * groupSize + j] = ASharesEval[j][i];
      BSharesEvalFlat[i * groupSize + j] = BSharesEval[j][i];
    }
  }

  vector<mTranscript<FieldType>> newTranscripts(5);
  vector<vector<vector<FieldType>>> newKingTranscripts(2);
  vector<vector<FieldType>> newRelayTranscripts;
  dotProd(king, groupSize, ASharesEvalFlat, BSharesEvalFlat,
          newTranscripts, newKingTranscripts, newRelayTranscripts);

  // append to passed-in dShares
  transcripts.insert(transcripts.end(),
                     newTranscripts.begin(), newTranscripts.end());
  if (_myId == king) {
    appendKingTranscripts(kingTranscripts, newKingTranscripts);
  }
  if (_disp.isRelayer(king, _myId)) {
    relayTranscripts.insert(relayTranscripts.end(),
                            newRelayTranscripts.begin(),
                            newRelayTranscripts.end());
  }
}

template <class FieldType> void ProtocolParty<FieldType>::inputPhase() {
  vector<int> sizes(_N, 0);
  for (int k = 0; k < _nInputs; k++) {
    assert(_circuit.getGates()[k].gateType == INPUT);
    // get the expected sizes from the other parties
    int inParty = _circuit.getGates()[k].party;
    sizes[inParty]++;
  }
  vector<FieldType> secrets;
  readMyInputs(sizes[_myId], secrets);
  
  vector<vector<FieldType>> myInputFullShares;
  makeTShares(secrets, myInputFullShares);
  _protoTape.recordDealtShares(myInputFullShares);
  vector<vector<byte>> sendBufs, recBufs(_N);
  sendBufs.resize(_N, vector<byte>(sizes[_myId] * _fieldByteSize));
  for (int i = 0; i < _N; i++) {
    _field->elementVectorToByteVector(myInputFullShares[i], sendBufs[i]);
    recBufs[i].resize(sizes[i] * _fieldByteSize, 0);
  }
  for (int i=0; i<_N; i++) {
    if (!_disp.isCorrupt(i)) {
      _comm.oneToAll(recBufs[i], sendBufs, i, false); // no relay
    }
  }
  // -- convert received bytes to shares
  vector<FieldType> inputShares;
  vector<unsigned int> startIdxs(_N, 0);
  for (int i = 0; i < _N; i++) {
    inputShares.resize(sizes[i]);
    for (int j = 0; j < sizes[i]; j++) { // receive actual shares
      inputShares[j] =
        _field->bytesToElement(recBufs[i].data() + (j * _fieldByteSize));
    }      
    // vv record received shares and their idx
    startIdxs[i] = _protoTape.getRecTapeLength();
    _protoTape.recordRecShares(inputShares, i);
  }
  // -- store shares in input gate order
  for (int k = 0; k < _nInputs; k++) {
    assert(_circuit.getGates()[k].gateType == INPUT);
    int inParty = _circuit.getGates()[k].party;
    int inWire = _circuit.getGates()[k].output;
    _gateShares[inWire] = _protoTape.getRecShare(startIdxs[inParty]++);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
makeRandShares(int nRands, vector<Share<FieldType>> &randShares) {
  int nBuckets = (nRands / (_N - _T)) + 1;
  randShares.resize(nBuckets * (_N - _T));
  
  vector<FieldType> secrets(nBuckets);
  for (int i=0; i<nBuckets; i++) {
    secrets[i] = _field->Random();
  }
  vector<vector<FieldType>> fullRandShares(_N, vector<FieldType>(nBuckets));
  makeTShares(secrets, fullRandShares);
  
  _protoTape.recordDealtShares(fullRandShares);
  vector<vector<byte>> sendBufs(_N, vector<byte>(nBuckets*_fieldByteSize));
  for (int i=0; i<_N; i++) {
    _field->elementVectorToByteVector(fullRandShares[i], sendBufs[i]);
  }
  vector<vector<byte>> recBufs(_N, vector<byte>(nBuckets*_fieldByteSize, 0));
  _comm.allToAll(sendBufs, recBufs, false); // no relay
  
  int count = 0;
  vector<FieldType> bufferIn(_N);
  vector<Share<FieldType>> sharesOut(_N - _T);
  for (int i=0; i<nBuckets; i++) { // improve: switch loop order ?
    unsigned int bucketIdxStart = _protoTape.getRecTapeLength();
    for (int j=0; j<_N; j++) {
      bufferIn[j] =
        _field->bytesToElement(recBufs[j].data() + i * _fieldByteSize);
      // record received share and its dealer id
      _protoTape.recordRecShare(bufferIn[j], j);
    }
    VDMMultShares(bucketIdxStart, _mat_vand_trans, bufferIn, sharesOut);
    for (int j=0; j<_N - _T; j++) {
      randShares[count++] = sharesOut[j];
    }
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
makeRandDoubleShares(int nRands, vector<Share<FieldType>> &randSharePairs) {
  int nBuckets = (nRands / (_N - _T)) + 1;
  randSharePairs.resize(nBuckets * (_N - _T) * 2);

  vector<FieldType> secrets;
  vector<vector<FieldType>> dSharesAll, sharesAll, randSharePairsAll;
  make2TShares(nBuckets, secrets, dSharesAll);
  _protoTape.recordDealt2TShares(dSharesAll);
  makeTShares(secrets, sharesAll);
  _protoTape.recordDealtShares(sharesAll);
  // mix single and double shares
  randSharePairsAll.resize(_N, vector<FieldType>(nBuckets*2));
  for (int j=0; j<_N; j++) {
    for (int i = 0; i < nBuckets; i++) {
      randSharePairsAll[j][i*2] = sharesAll[j][i]; // even slots: single share 
      randSharePairsAll[j][i*2+1] = dSharesAll[j][i]; // odd slots: double Share
    }
  }
  vector<vector<byte>> sendBufs(_N, vector<byte>(nBuckets*_fieldByteSize*2));
  for (int i=0; i<_N; i++) {
    _field->elementVectorToByteVector(randSharePairsAll[i], sendBufs[i]);
  }
  vector<vector<byte>> recBufs(_N, vector<byte>(nBuckets*_fieldByteSize*2, 0));
  _comm.allToAll(sendBufs, recBufs, false); // no relay
  int count = 0;
  vector<FieldType> bufferIn(_N), bufferIn2(_N);
  vector<Share<FieldType>> sharesOut(_N-_T), sharesOut2(_N-_T);
  for (int i=0; i<nBuckets; i++) {
    unsigned int bucketIdxStart = _protoTape.getRecTapeLength();
    unsigned int bucketIdxStart2 = _protoTape.getRec2TTapeLength();
    for (int j=0; j<_N; j++) {
      bufferIn[j] =
        _field->bytesToElement(recBufs[j].data() + 2*i*_fieldByteSize);
      _protoTape.recordRecShare(bufferIn[j], j);
      bufferIn2[j] =
        _field->bytesToElement(recBufs[j].data() + (2*i+1)*_fieldByteSize);
      _protoTape.recordRec2TShare(bufferIn2[j], j);
    }
    VDMMultShares(bucketIdxStart, _mat_vand_trans, bufferIn, sharesOut);
    VDMMultShares(bucketIdxStart2, _mat_vand_trans, bufferIn2, sharesOut2);
    for (int j=0; j<_N-_T; j++) {
      randSharePairs[count*2] = sharesOut[j];
      randSharePairs[count*2+1] = sharesOut2[j];
      count++;
    }
  }
}

template <class FieldType>
bool ProtocolParty<FieldType>::
verifySharing() {
  vector<Share<FieldType>> sigma(1);
  getPadShares(1, sigma);
  FieldType lambda = challenge();
  FieldType lambdaI = lambda;

  vector<Share<FieldType>> toVerify;
  _protoTape.makeToVerifyShares(toVerify);
  
  for(Share<FieldType> share : toVerify) {
    share.multByConst(lambdaI);
    sigma[0] += share;
    lambdaI *= lambda;
  }

  vector<byte> sendBuf(1 * _fieldByteSize);
  vector<vector<byte>> recBufs(_N, vector<byte>(1 * _fieldByteSize));
  _field->elementToBytes(sendBuf.data(), sigma[0]._value);
  // _field->elementVectorToByteVector(sigma, sendBuf);
  _comm.allBroadcast(sendBuf, recBufs);
  vector<FieldType> sigmaFullShare(_N, _zero);
  for (int curId = 0; curId < _N; curId++) {
    sigmaFullShare[curId] = 
      _field->bytesToElement(recBufs[curId].data());
  }

  // ---- check 1: is [sigma] consistent? (else analyzeSharing)
  if (!checkTShares(sigmaFullShare)) {
    analyzeSharing(sigmaFullShare, sigma[0]);
    return false;
  }
  return true;
}

template <class FieldType>
bool ProtocolParty<FieldType>::
tag() {
  if (_disp.getNewDisp()) {     // improve: new disp flag per party
    _reckeyShares.resize(_authKeySize);
    _dealtKeyShares.resize(_authKeySize);
    _keys.resize(_authKeySize);
    // improve: batch this
    for (int i=0; i<_authKeySize; i++) {
      keyDistribution(_keys[i], _reckeyShares[i], _dealtKeyShares[i]);
    }
    if (!checkKey()) {
      return false;
    }
    // cout << "[" << _myId << "] ran keyDistribution" << endl;
  }

  vector<vector<FieldType>> authBatchs;
  vector<int> authDealers;
  int nBatch = _protoTape.makeAuthBatch(authBatchs, authDealers);
  vector<vector<FieldType>> recOShares, recTauShares, tauShares,
    dealtOShares, myNuShares;
  singleTagComp(authBatchs, authDealers, recOShares, recTauShares,
                tauShares, dealtOShares, myNuShares);
  FieldType lambda = challenge();
  if (!checkTag(authBatchs, authDealers, recOShares, recTauShares,
                tauShares, dealtOShares, myNuShares)) {
    return false;
  }

  _protoTape.updateAuthOffset();
  return true;
}

template <class FieldType>
void ProtocolParty<FieldType>::
keyDistribution(vector<FieldType>& keys,
                vector<vector<FieldType>>& reckeyShares,
                vector<vector<FieldType>>& sentKeyShares) {
  vector<int> DispSet, NonDispSet;
  int nDisp = _disp.dispAndNonDispSet(_myId, DispSet, NonDispSet);
  int nNonDisp = _N - nDisp;
  vector<vector<byte>> sendBufs(_N, vector<byte>(_N*_fieldByteSize, 0));
  vector<vector<FieldType>> keySharesAll(_N, vector<FieldType>(_N, _zero));
  keys.clear();
  keys.resize(_N, _zero);
  reckeyShares.clear();
  reckeyShares.resize(_N, vector<FieldType>(_N, _zero));
  if (nDisp < _T) {
    vector<FieldType> secrets(nNonDisp);
    for (int i=0; i<nNonDisp; i++) {
      int pid = NonDispSet[i];
      keys[pid] = _field->Random();
      secrets[i] = keys[pid];
    }
    vector<vector<FieldType>> keySharesNonDisp;
    makeTwistTShares(secrets, NonDispSet, keySharesNonDisp);
    for (int i=0; i<nNonDisp; i++) {
      for (int j=0; j<_N; j++) {
        keySharesAll[j][NonDispSet[i]] = keySharesNonDisp[j][i];
      }
    }
    for (int i=0; i<_N; i++) {
      _field->elementVectorToByteVector(keySharesAll[i], sendBufs[i]);
    }
  } // else : all shares (keys) are zero

  
  vector<vector<byte>> recBufs(_N, vector<byte>(_N*_fieldByteSize, 0));
  _comm.allToAll(sendBufs, recBufs, false); // no relay
  for (int i=0; i<_N; i++) {
    // recBufs[i] are keys from i of NonDispSet of i
    if (_disp.isDisp(_myId, i)) {
      continue;                 // no key shares from disp party
    }
    vector<int> DispSetI, NonDispSetI;
    int nDispI = _disp.dispAndNonDispSet(i, DispSet, NonDispSet);
    if (nDispI >= _T) {         // all key Shares from him are zero
      continue;
    }
    for (int pid : NonDispSetI) {
      reckeyShares[i][pid] =
        _field->bytesToElement(recBufs[i].data() + pid * _fieldByteSize);
    }
  }
}

template <class FieldType>
bool ProtocolParty<FieldType>::checkKey() {
  vector<FieldType> padKeys;
  vector<vector<FieldType>> padRecKeyShares;
  vector<vector<FieldType>> padSentKeyShares;
  keyDistribution(padKeys, padRecKeyShares, padSentKeyShares);
  FieldType lambda = challenge();
  for (int i=0; i<_N; i++) {
    for (int j=0; j<_N; j++) {
      FieldType lambdaI = lambda;
      for (int k = 0; k < _authKeySize; k++) {
        padRecKeyShares[i][j] += _reckeyShares[k][i][j] * lambdaI;
        lambdaI *= lambda;
      }
    }
  }

  vector<byte> sendBuf;
  for (int i=0; i<_N; i++) {
    vector<byte> tmp(_N * _fieldByteSize);
    _field->elementVectorToByteVector(padRecKeyShares[i], tmp);
    sendBuf.insert(sendBuf.end(), tmp.begin(), tmp.end());
  }
  vector<vector<byte>> sendBufs(_N, sendBuf);
  vector<vector<byte>> recBufs(_N, vector<byte>(sendBuf.size()));
  _comm.allToAll(sendBufs, recBufs, false); // no relay

  vector<FieldType> keyIJShares(_N);
  vector<HIM<FieldType>> mat_twist_inv_vecI;
  for (int i=0; i<_N; i++) {    // key from i
    if (_disp.isCorrupt(i)) {
      continue;
    }
    matForTwist(i, mat_twist_inv_vecI, true);
    for (int j=0; j<_N; j++) {   // key to j
      if (_disp.isDisp(i, j)) {
        continue;
      }
      vector<int> TwistSetIJ;
      int nDispI = _disp.twistSet(i, j, TwistSetIJ);
      for (int k=0; k<_N; k++) { // shares from k
        keyIJShares[k] =
          _field->bytesToElement(recBufs[k].data() + (i*_N+j)*_fieldByteSize);
      }
      // check that shares are consistent
      vector<FieldType> dispShares(_T+1, _zero), nonDispShares(_T+1);
      for (int k=_T-1; k<2*_T; k++) {
        nonDispShares[k-(_T-1)] = keyIJShares[TwistSetIJ[k]];
      }
      mat_twist_inv_vecI[j].MatrixMult(nonDispShares, dispShares);
      for (int k=0; k<nDispI; k++) {
        assert(dispShares[k] == _zero);
      }
      for (int k=nDispI; k<_T; k++) {
        assert(dispShares[k] == keyIJShares[TwistSetIJ[k]]);
      }
    }
  }
  // TODO: minor: implement dispute control messages
  return true;
}

// improve: refactor
template <class FieldType>
void ProtocolParty<FieldType>::
singleTagComp(vector<vector<FieldType>>& authBatchs, // in
              vector<int>& authDealers,              // in
              vector<vector<FieldType>>& recOShares, // out
              vector<vector<FieldType>>& recTauShares,  // out
              vector<vector<FieldType>>& tauShares,  // out
              vector<vector<FieldType>>& dealtOShares,
              vector<vector<FieldType>>& myNuShares) { // out
  int nBatch = authBatchs.size();
  // -- generate nu tags as verifyer
  // sends _N nu shares per batch, the [i]th nu share is for party i
  myNuShares.assign(_N, vector<FieldType>(nBatch * _N, _zero));
  vector<int> dealerSizes(_N, 0); // how many batches for dealer?
  for (int batchIdx = 0; batchIdx < nBatch; batchIdx++) {
    int dealer = authDealers[batchIdx];
    dealerSizes[dealer]++;
    int batchStart = _N * batchIdx;
    for (int verifyee = 0; verifyee<_N; verifyee++) {
      if (_disp.isDisp(_myId, verifyee) ||_disp.isDisp(dealer, verifyee)) {
        _protoTape.recordSentTag(verifyee, _zero);
        continue;               // it's share is left blank, and tag is zero
      }
      vector<FieldType> nuFullShare;
      FieldType nuTag = makeTwist2TShare(verifyee, nuFullShare);
      for (int recParty=0; recParty<_N; recParty++) {
        myNuShares[recParty][batchStart + verifyee] = nuFullShare[recParty];
      }
      _protoTape.recordSentTag(verifyee, nuTag);
    }
  }
  // -- get nu tags as participant
  vector<vector<byte>>
    sendBufs(_N, vector<byte>(nBatch * _N * _fieldByteSize));
  vector<vector<byte>> recBufs = sendBufs;
  for (int i = 0; i<_N; i++) {
    _field->elementVectorToByteVector(myNuShares[i], sendBufs[i]);
  }
  _comm.allToAll(sendBufs, recBufs, false); // no relay
  vector<vector<FieldType>> recNuShares(_N);
  for (int verifier=0; verifier<_N; verifier++) {
    bytesToVector(recBufs[verifier], recNuShares[verifier]);
  }

  // -- generate o tags as dealer
  int nMyShares = dealerSizes[_myId];
  // sends _N * _N o shares per my batch
  dealtOShares.assign(_N, vector<FieldType>(nMyShares * _N * _N, _zero));
  int count = 0;
  for (int batchIdx = 0; batchIdx < nBatch; batchIdx++) {
    if (_myId != authDealers[batchIdx]) {
      continue;                 // only generates as dealer
    }
    int batchStart = _N * _N * (count++);
    for (int verifier = 0; verifier<_N; verifier++) {
      int verifierStart = batchStart + _N * verifier;
      for (int verifyee = 0; verifyee<_N; verifyee++) {
        if (_disp.isDisp(_myId, verifyee) || _disp.isDisp(verifier, verifyee)) {
          continue;
        }
        vector<FieldType> oFullShare;
        makeTwist2TOShare(verifyee, oFullShare);
        for (int recParty=0; recParty<_N; recParty++) {
          dealtOShares[recParty][verifierStart+verifyee] = oFullShare[recParty];
        }
      }
    }
  }
  // -- get o shares as participant
  sendBufs.assign(_N, vector<byte>(nMyShares * _N * _N * _fieldByteSize));
  for (int i = 0; i<_N; i++) {
    _field->elementVectorToByteVector(dealtOShares[i], sendBufs[i]);
  }
  vector<byte> recBuf;
  recOShares.assign(_N, vector<FieldType>());
  for (int dealer=0; dealer<_N; dealer++) {
    int nDealerShares = dealerSizes[dealer];
    recBuf.assign(nDealerShares * _N * _N * _fieldByteSize, 0);
    if (!_disp.isCorrupt(dealer)) {
      _comm.oneToAll(recBuf, sendBufs, dealer, false);
    }
    bytesToVector(recBuf, recOShares[dealer]);
  }
  // -- compute tag shares as participant
  // each party: per batch _N tagShares
  // the [verifier]th tagShare in a batch is for the verifier
  tauShares.assign(_N, vector<FieldType>(nBatch*_N, _zero));
  vector<int> dealerIdxs(_N, 0);
  for (int batchIdx = 0; batchIdx < nBatch; batchIdx++) {
    int dealer = authDealers[batchIdx];
    int oDStart = _N * _N * (dealerIdxs[dealer]++);
    for (int verifier = 0; verifier<_N; verifier++) {
      int oVStart = oDStart + _N * verifier;
      for (int verifyee = 0; verifyee<_N; verifyee++) {
        if (_disp.isDisp(dealer,verifyee) || _disp.isDisp(verifier, verifyee)) {
          continue;
        }
        // compute tag for (v, i) pair
        FieldType tagShare = recOShares[dealer][oVStart + verifyee];
        tagShare += recNuShares[verifier][_N * batchIdx + verifyee];
        for (int keyIdx=0; keyIdx < _authKeySize; keyIdx++) {
          _reckeyShares[keyIdx][verifier][verifyee] * authBatchs[batchIdx][keyIdx];
        }
        tauShares[verifyee][_N * batchIdx + verifier] = tagShare;
      }
    }
  }
  sendBufs.assign(_N, vector<byte>(nBatch * _N * _fieldByteSize));
  recBufs = sendBufs;
  _comm.allToAll(sendBufs, recBufs, true); // use relay

  // -- gather tag shares as verifyee
  recTauShares.resize(_N);
  for (int party = 0; party<_N; party++) {
    bytesToVector(recBufs[party], recTauShares[party]);
  }
  for (int batchIdx = 0; batchIdx < nBatch; batchIdx++) {
    int batchStart = _N * batchIdx;
    for (int verifier = 0; verifier<_N; verifier++) {
      vector<FieldType> fullShare(_N);
      for (int party = 0; party <_N; party++) {
        fullShare[party] = recTauShares[party][batchStart + verifier];
      }
      FieldType tauTag = getTwistSecret(fullShare, _myId);
      _protoTape.recordRecTag(verifier, tauTag);
    }
  }
}


template <class FieldType>
void ProtocolParty<FieldType>::
baseSharing(vector<vector<FieldType>>& randomShares) {
  // TODO: implement actual baseSharing from BSFO12
}

template <class FieldType>
void ProtocolParty<FieldType>::
solveAccuseTag(int verifier, int accused, int dealer, int suspect,
               vector<byte> &relayBuf, vector<vector<byte>> &relayBufs) {
  int relayer = _disp.relayer(verifier, suspect);
  vector<byte> relayerMsg, recMsgBuf(_dispMsgSize);
  if (_myId == relayer) {
    // NOTE: if verifier in Disp(suspect), then
    // its tau share is equal to oshare
    FieldType sTauShare = _field->bytesToElement(relayBufs[suspect].data() +
                                                 _fieldByteSize);
    FieldType sOShare = _field->bytesToElement(relayBufs[suspect].data() +
                                               2 * _fieldByteSize);
    if (sTauShare != sOShare) {
      makeOkMsg(relayerMsg);    // agrees w/ verifier
    } else if (suspect == dealer) {
      vector<vector<FieldType>> // for verifier
        oSigmaDealerFullShares(_N, vector<FieldType>(_authKeySize + 1)); 
      for (int partyId=0; partyId<_N; partyId++) {
        oSigmaDealerFullShares[partyId][0] =
          _field->bytesToElement(relayBuf.data() + partyId * _fieldByteSize);
        for (int keyId=0; keyId<_N; keyId++) {
          oSigmaDealerFullShares[partyId][keyId] =
            _field->bytesToElement(relayBuf.data() + _N * _fieldByteSize +
                                   (partyId*_authKeySize+keyId)*_fieldByteSize);
        }
      }

      if (!testDealerOSigma(dealer, accused, oSigmaDealerFullShares)) {
        makeOkMsg(relayerMsg);
      }
    }
  }
  _comm.singleBroadcast(relayer, relayerMsg, recMsgBuf);
  if (recMsgBuf[0] == 0) {
    _disp.addDispPairs(relayer, suspect);
  } else {
    _disp.addDispPairs(relayer, verifier);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
checkRelayer(int relayer, int id1, int id2,
             FieldType share1, FieldType rShare) {
  vector<byte> relayerMsg, recMsgBuf(_dispMsgSize);
  if (_myId == relayer) {
    if (rShare == id1) {
      makeOkMsg(relayerMsg);    // agrees w/ id1
    } else {
      makeRejectMsg(relayerMsg);
    }
  }
  _comm.singleBroadcast(relayer, relayerMsg, recMsgBuf);
  if (recMsgBuf[0]) {
    _disp.addDispPairs(relayer, id2);
  } else {
    _disp.addDispPairs(relayer, id1);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
solveOpenTag(int verifier, int accused, vector<byte> &openMsg,
             FieldType myTauShare, vector<FieldType>& recTauFullShare,
             FieldType myOShare, vector<FieldType>& sentOFullShare,
             vector<FieldType>& mySigmaShares,
             vector<vector<FieldType>>& sentSigmaFullShares,
             vector<byte>& relayBuf, vector<vector<byte>>& relayBufs) {
  int id1, id2, keyIdx;
  FieldType vShare1, vShare2;   // from verifier
  readOpenMsgTag(openMsg, id1, id2, keyIdx, vShare1, vShare2);
  vector<byte> sendBuf(_fieldByteSize), recBuf(_fieldByteSize);
  if (id1 == accused) {         // accused tau vs. party tau
    if (_myId == accused) {
      _field->elementToBytes(sendBuf.data(), recTauFullShare[id2]);
    }
    _comm.singleBroadcast(id1, sendBuf, recBuf);
    if (_field->bytesToElement(recBuf.data()) != vShare1) {
      _disp.addDispPairs(accused, verifier);
      return;
    } // else check id2

    if(_myId == id2) {
      _field->elementToBytes(sendBuf.data(), myTauShare);
    }
    _comm.singleBroadcast(id2, sendBuf, recBuf);
    FieldType pShare2 = _field->bytesToElement(recBuf.data());
    if (pShare2 == vShare2) {
      _disp.addDispPairs(id2, accused);
    } // else check relay
    
    int relayer = _disp.relayer(verifier, id2);
    checkRelayer(relayer, verifier, id2, vShare2,
                 _field->bytesToElement(relayBufs[id2].data()));
    return;
  }

  FieldType dShare1, pShare2, rShare1, rShare2;
  if (keyIdx < 0) {             // dealer o vs. party o
    dShare1 = sentOFullShare[id2];
    pShare2 = myOShare;
    rShare1 = _field->bytesToElement(relayBuf.data() + id2 * _fieldByteSize);
    rShare2 = _field->bytesToElement(relayBufs[id2].data() + _fieldByteSize);
  } else {                      // dealer sigma vs. party sigma
    dShare1 = sentSigmaFullShares[id2][keyIdx];
    pShare2 = mySigmaShares[keyIdx];
    rShare1 = _field->bytesToElement(relayBuf.data() + _N * _fieldByteSize +
                                     (id2*_authKeySize+keyIdx)*_fieldByteSize);
    rShare2 = _field->bytesToElement(relayBufs[id2].data() +
                                     (2+keyIdx)*_fieldByteSize);
  }
  if (_myId == id1) {
    _field->elementToBytes(sendBuf.data(), dShare1);
  }
  _comm.singleBroadcast(id1, sendBuf, recBuf);
  dShare1 = _field->bytesToElement(recBuf.data());
  if (dShare1 != vShare1) {
    if (!_disp.isDisp(verifier, id1)) {
      _disp.addDispPairs(verifier, id1);
      return;
    }
    // else, check relayer
    int relayer = _disp.relayer(verifier, id1);
    checkRelayer(relayer, verifier, id1, vShare1, rShare1);
    return;
  } // else, check id2
  if (_myId == id2) {
    _field->elementToBytes(sendBuf.data(), pShare2);
  }
  _comm.singleBroadcast(id2, sendBuf, recBuf);
  pShare2 = _field->bytesToElement(recBuf.data());
  if (pShare2 == vShare2) {
    _disp.addDispPairs(id1, id2);
    return;
  } // else, check relayer
  int relayer = _disp.relayer(verifier, id2);
  checkRelayer(relayer, verifier, id2, vShare2, rShare2);
}



template <class FieldType>
void ProtocolParty<FieldType>::
flSingleDealer(int verifier, int accused, int dealer, FieldType lambda,
               vector<int>& authDealers, vector<vector<FieldType>>& authBatchs, 
               vector<vector<FieldType>>& recOShares,
               vector<vector<FieldType>>& recTauShares,
               vector<vector<FieldType>>& tauShares,
               vector<vector<FieldType>>& dealtOShares,
               vector<vector<FieldType>>& myNuShares) {
  FieldType lambdaW = lambda;
  int nBatch = authBatchs.size();
  
  FieldType dTauShare = _zero;
  FieldType dRecOShare = _zero;
  vector<FieldType> dSigmaShares(_authKeySize, _zero);
  int dCount = 0;
  // -- all parties prepare shares for [tau'(d)], [sigma'(d)], and [o'(d)]
  for (int batchIdx=0; batchIdx<nBatch; batchIdx++) {
    if (authDealers[batchIdx] == dealer) {
      dTauShare += tauShares[accused][verifier] * lambdaW;
      dRecOShare += recOShares[dealer][_N*_N*dCount + _N*verifier + accused];
      for (int keyIdx=0; keyIdx < _authKeySize; keyIdx++) {
        dSigmaShares[keyIdx] += authBatchs[batchIdx][keyIdx] * lambdaW;
      }
      lambdaW *= lambda;
      dCount++;
    } // else, ignore
  }
  vector<byte> sendBuf(2 * _fieldByteSize), tmp(_authKeySize * _fieldByteSize);
  _field->elementToBytes(sendBuf.data(), dTauShare);
  _field->elementToBytes(sendBuf.data() + _fieldByteSize, dRecOShare);
  _field->elementVectorToByteVector(dSigmaShares, tmp);
  sendBuf.insert(sendBuf.end(), tmp.begin(), tmp.end());
  vector<vector<byte>> relayBufs,
    recBufs(_N, vector<byte>((2+_authKeySize) * _fieldByteSize));
  _comm.allToOneStore(sendBuf, relayBufs, recBufs, verifier);
  vector<vector<FieldType>> tauOSigmaFullShares; // for verifier
  if (_myId == verifier) {
    tauOSigmaFullShares.resize(_N, vector<FieldType>(_authKeySize+2));
    for (int partyId=0; partyId<_N; partyId++) {
      bytesToVector(recBufs[partyId], tauOSigmaFullShares[partyId]);
    }
  }

  lambdaW = lambda;
  // -- dealer prepare full Shares of [sigma'(d)] and [o'(d)]
  vector<FieldType> dSentOFullShare;
  vector<vector<FieldType>> dSentSigmaFullShares;
  if (_myId == dealer) {
    dSentOFullShare.resize(_N, _zero);
    dSentSigmaFullShares.resize(_N, vector<FieldType>(_authKeySize, _zero));
    for (int i=0; i<dCount; i++) {
      vector<vector<FieldType>> sentFullShares;
      _protoTape.getDealerBatchFullShares(i, sentFullShares);
      for (int partyId=0; partyId<_N; partyId++) {
        dSentOFullShare[partyId] +=
          dealtOShares[partyId][i*_N*_N + verifier*_N + accused];
        for (int keyIdx=0; keyIdx<_authKeySize; keyIdx++) {
          dSentSigmaFullShares[partyId][keyIdx] +=
            sentFullShares[partyId][keyIdx] * lambdaW;
        }
      }
      lambdaW *= lambda;
    }
    sendBuf.resize(_N);
    _field->elementVectorToByteVector(dSentOFullShare, sendBuf);
    vector<byte> tmp(_authKeySize);
    for (int partyId=0; partyId<_N; partyId++) {
      _field->elementVectorToByteVector(dSentSigmaFullShares[partyId], tmp);
      sendBuf.insert(sendBuf.end(), tmp.begin(), tmp.end());
    }
  }
  vector<byte> recBuf, relayBuf;
  if (_myId == verifier) {
    recBuf.resize((_N+_N*_authKeySize) * _fieldByteSize);
  }
  _comm.oneToOneStore(dealer, verifier, sendBuf, relayBuf, recBuf);
  vector<vector<FieldType>> // for verifier
    oSigmaDealerFullShares(_N, vector<FieldType>(_authKeySize + 1));   
  if (_myId == verifier) {
    for (int partyId=0; partyId<_N; partyId++) {
      oSigmaDealerFullShares[partyId][0] =
        _field->bytesToElement(recBuf.data() + partyId * _fieldByteSize);
      for (int keyId=0; keyId<_N; keyId++) {
        oSigmaDealerFullShares[partyId][keyId] =
          _field->bytesToElement(recBuf.data() + _N * _fieldByteSize +
                                 (partyId*_authKeySize+keyId)*_fieldByteSize);
      }
    }
  }
  
  lambdaW = lambda;
  // -- accused prepare received full share for [tau'(d)]
  vector<FieldType> dRecTauFullShare;
  if (_myId == accused) {
    dRecTauFullShare.resize(_N, _zero);
    for (int batchIdx=0; batchIdx<nBatch; batchIdx++) {
      if (authDealers[batchIdx] == dealer) {
        for(int partyId=0; partyId<_authKeySize; partyId++) {
          dRecTauFullShare[partyId] += lambdaW * 
            recTauShares[partyId][batchIdx*_N + verifier];
        }
        lambdaW *= lambda;
      }
    }
    sendBuf.resize(_N * _fieldByteSize);
    _field->elementVectorToByteVector(dRecTauFullShare, sendBuf);
  }
  recBuf.resize(_N * _fieldByteSize); // NOTE: no relay should happen
  _comm.oneToOneStore(accused, verifier, sendBuf, relayBuf, recBuf);
  vector<FieldType> tauAccusedFullShare; // for verifier
  if (_myId == verifier) {
    bytesToVector(recBuf, tauAccusedFullShare);
  }
  
  lambdaW = lambda;
  // -- verifier checks everyting
  vector<byte> verifierMsg, recMsgBuf(_dispMsgSize);
  if (_myId == verifier) {
    // -- verifier prepare full share for [nu'(d)]
    vector<FieldType> nuMyFullShare;
    nuMyFullShare.resize(_N, _zero);
    for (int batchIdx=0; batchIdx<nBatch; batchIdx++) {
      if (authDealers[batchIdx] == dealer) {
        for(int partyId=0; partyId<_authKeySize; partyId++) {
          nuMyFullShare[partyId] += lambdaW * 
            myNuShares[partyId][batchIdx*_N + accused];
        }
        lambdaW *= lambda;
      }
    }
    checkDealerTagComp(dealer, accused, verifierMsg,
                       tauOSigmaFullShares, oSigmaDealerFullShares,
                       tauAccusedFullShare, nuMyFullShare);
  }
  _comm.singleBroadcast(verifier, verifierMsg, recMsgBuf);
  if (recMsgBuf[0] == 2) {
    int suspect = recMsgBuf[1];
    if (!_disp.isDisp(verifier, suspect)) {
      _disp.addDispPairs(suspect, verifier);
      return;
    }
    solveAccuseTag(verifier, accused, dealer, suspect, relayBuf, relayBufs);
    return;
  } // else, should be open message (special format)
  solveOpenTag(verifier, accused, recMsgBuf, dTauShare, dRecTauFullShare,
               dRecOShare, dSentOFullShare, dSigmaShares, dSentSigmaFullShares,
               relayBuf, relayBufs);
}

template <class FieldType>
int ProtocolParty<FieldType>::
flTag(int verifier, int accused, FieldType lambda, vector<FieldType>& lambdaWs,
      vector<vector<FieldType>>& dealerTaus,
      vector<vector<FieldType>>& dealerBatchs,
      vector<vector<FieldType>>& dealerNus,
      FieldType recTau, vector<FieldType>& recBatch){
  if (_myId == accused) {
    vector<FieldType> tauAndBatch(_N);
    for (int dealer=0; dealer<_N; dealer++) {
      tauAndBatch[dealer] = dealerTaus[dealer][verifier];
    }
    for (int dealer=0; dealer<_N; dealer++) {
      tauAndBatch.insert(tauAndBatch.end(), dealerBatchs[dealer].begin(),
                         dealerBatchs[dealer].end());
    }
    vector<byte> sendBuf((_N + _N *_authKeySize) * _fieldByteSize);
    _field->elementVectorToByteVector(tauAndBatch, sendBuf);
    _comm.write(sendBuf, verifier);
  }

  vector<byte> verifierMsg, recMsgBuf(_dispMsgSize);
  if (_myId == verifier) {
    vector<byte> recBuf((_N + _N *_authKeySize) * _fieldByteSize);
    _comm.read(recBuf, accused);
    checkDealerBatches(accused, verifierMsg, recBuf, lambdaWs, dealerNus,
                       recTau, recBatch);
  }
  _comm.singleBroadcast(verifier, verifierMsg, recMsgBuf);
  assert(recMsgBuf[0] == 2);
  if (recMsgBuf[1] == accused) {
    _disp.addDispPairs(verifier, accused);
    return -1;
  } else {
    return recMsgBuf[1];
  }
}

template <class FieldType>
bool ProtocolParty<FieldType>::
checkTag(vector<vector<FieldType>>& authBatchs, vector<int>& authDealers,
         vector<vector<FieldType>>& recOShares,
         vector<vector<FieldType>>& recTauShares,
         vector<vector<FieldType>>& tauShares,
         vector<vector<FieldType>>& dealtOShares,
         vector<vector<FieldType>>& myNuShares) {
  // TODO: use baseSharing to prepare a pad batch for every dealer
  int nBatch = authDealers.size();
  int authBatchStart = _protoTape.getPrevNumAuthBatches();
  FieldType lambda = challenge();
  
  vector<FieldType> lambdaWs(_N, lambda);
  vector<vector<FieldType>>
    dealerTaus(_N, vector<FieldType>(_N, _zero)),
    dealerBatchs(_N, vector<FieldType>(_authKeySize, _zero)),
    dealerNus(_N, vector<FieldType>(_N, _zero));
  for (int batchIdx = 0; batchIdx < nBatch; batchIdx++) {
    int dealer = authDealers[batchIdx];
    // -- computes combined sigma and tau as verifyee
    vector<FieldType> tags;
    _protoTape.getTagVals(authBatchStart + batchIdx, tags);
    for (int verifier = 0; verifier < _N; verifier++) {
      dealerTaus[dealer][verifier] += lambdaWs[dealer] * tags[verifier];
    }
    for (int keyIdx=0; keyIdx < _authKeySize; keyIdx++) {
      dealerBatchs[dealer][keyIdx] +=
        lambdaWs[dealer] * authBatchs[batchIdx][keyIdx];
    }
    // -- computes combined nu as verifier
    for (int verifyee = 0; verifyee < _N; verifyee++) {
      dealerNus[dealer][verifyee] += lambdaWs[dealer] *
        _protoTape.getPadVal(verifyee, authBatchStart + batchIdx);
    }
    lambdaWs[dealer] *= lambda;
  }
  FieldType lambdaD = *(_field->GetOne());
  vector<FieldType> myTaus(_N, _zero), myNus(_N, _zero),
    myBatch(_authKeySize, _zero);
  for (int dealer=0; dealer<_N; dealer++) {
    for (int verifier = 0; verifier < _N; verifier++) {
      myTaus[verifier] +=
        dealerTaus[dealer][verifier] * lambdaD * lambdaWs[dealer];
    }
    for (int keyIdx=0; keyIdx < _authKeySize; keyIdx++) {
      myBatch[keyIdx] +=
        dealerBatchs[dealer][keyIdx] * lambdaD * lambdaWs[dealer];
    }
    for (int verifyee = 0; verifyee < _N; verifyee++) {
      myNus[verifyee] +=
        dealerNus[dealer][verifyee] * lambdaD * lambdaWs[dealer];
    }
    lambdaD *= lambda;
  }
  // -- gether taus and sigmas as verifier
  vector<vector<byte>> sendBufs(_N, vector<byte>(_fieldByteSize));
  vector<vector<byte>> recBufs = sendBufs;
  for (int verifier = 0; verifier < _N; verifier++) {
    _field->elementToBytes(sendBufs[verifier].data(), myTaus[verifier]);
  }
  _comm.allToAll(sendBufs, recBufs, false);
  vector<FieldType> recTaus(_N);
  for (int verifyee = 0; verifyee < _N; verifyee++) {
    recTaus[verifyee] = _field->bytesToElement(recBufs[verifyee].data());
  }
  sendBufs.assign(_N, vector<byte>(_authKeySize * _fieldByteSize));
  recBufs = sendBufs;
  for (int verifier = 0; verifier < _N; verifier++) {
    _field->elementVectorToByteVector(myBatch, sendBufs[verifier]);
  }
  _comm.allToAll(sendBufs, recBufs, false);
  vector<vector<FieldType>> recBatchs(_N);
  for (int verifyee = 0; verifyee < _N; verifyee++) {
    bytesToVector(recBufs[verifyee], recBatchs[verifyee]);
  }
  // -- check received tau and sigma
  vector<byte> verifierMsg;
  vector<vector<byte>> recMsgBufs(_N, vector<byte>(_dispMsgSize));
  makeOkMsg(verifierMsg);
  for (int verifyee = 0; verifyee < _N; verifyee++) {
    FieldType calcTag = myNus[verifyee];
    for (int keyIdx=0; keyIdx < _authKeySize; keyIdx++) {
      calcTag += recBatchs[verifyee][keyIdx] * _keys[keyIdx][verifyee];
    }
    if (calcTag != recTaus[verifyee]) {
      assert(!_disp.isDisp(_myId, verifyee));
      makeAccuseMsg(verifierMsg, verifyee);
      break;
    }
  }
  // -- let each verifier broadcast their msg
  _comm.allBroadcast(verifierMsg, recMsgBufs);
  for (int verifier = 0; verifier < _N; verifier++) {
    if (recMsgBufs[verifier][0] != 0) {
      int accused = recMsgBufs[verifier][1];
      int dealer =
        flTag(verifier, accused, lambda, lambdaWs, dealerTaus, dealerBatchs,
              dealerNus, recTaus[accused], recBatchs[accused]);
      if (dealer >= 0) {
        flSingleDealer(verifier, accused, dealer, lambda, authDealers, authBatchs, 
                       recOShares, recTauShares, tauShares, dealtOShares, myNuShares);
      }
      return false;
    }
  }
  return true;
}

// TODO: impl (with analyzeSharing)
template <class FieldType>
bool ProtocolParty<FieldType>::verifyOutput() {
  return true;
}


template <class FieldType>
void ProtocolParty<FieldType>::initializeMat() {
  vector<FieldType> beta_zero(1, _zero); // a sinlge zero element
  vector<FieldType> alpha_n(_N); // N distinct non-zero field elements
  vector<FieldType> alpha_2n(_N*2);
  for (int i = 0; i < _N; i++) {
    alpha_n[i] = _field->GetElement(i + 1);
    alpha_2n[i] = alpha_n[i];
    alpha_2n[i+_N] = _field->GetElement(i+_N+1);
  }
  _alpha_2km1.resize(_K*2 - 1);
  for (int i=0; i< _K*2 - 1; i++) {
    _alpha_2km1[i] = _field->GetElement(i+1);
  }
  _alpha_k.resize(_K);
  for (int i=0; i< _K; i++) {
    _alpha_k[i] = _field->GetElement(i+1);
  }

  // Interpolate first T+1 positions (deg = T)
  _mat_check_tp1.allocate(_N - (1+_T), _T + 1, _field); 
  _mat_check_tp1.InitHIMVectorAndsizes(alpha_n, _T + 1, _N - _T - 1);

  // Interpolate first K positions (deg = K-1)
  _mat_check_km1.allocate(2*_K-1 - _K, _K, _field);
  _mat_check_km1.InitHIMVectorAndsizes(_alpha_2km1, _K, 2*_K-1 - _K);

  // Interpolate from all _N positions to 0
  _mat_use_n.allocate(1, _N, _field);
  _mat_use_n.InitHIMByVectors(alpha_n, beta_zero);

  matFor2TTwist(_mat_use_n_twist_vec);
  
  _mat_vand_trans.allocate(_N, _N, _field);
  _mat_vand_trans.InitVDMTranspose();
}

template <class FieldType>
int ProtocolParty<FieldType>::processNonMult(int layerStart, int layerEnd) {
  int count = 0;
  for (int k = layerStart; k < layerEnd; k++) {
    if (_circuit.getGates()[k].gateType == ADD) {
      _gateShares[_circuit.getGates()[k].output] +=
          _gateShares[_circuit.getGates()[k].input2];
    }
    else if (_circuit.getGates()[k].gateType == SUB)
    {
      _gateShares[_circuit.getGates()[k].output] -=
          _gateShares[_circuit.getGates()[k].input2];
    } else if (_circuit.getGates()[k].gateType == SCALAR) {
      long scalar(_circuit.getGates()[k].input2);
      FieldType e = _field->GetElement(scalar);
      _gateShares[_circuit.getGates()[k].output] =
        _gateShares[_circuit.getGates()[k].input1];
      _gateShares[_circuit.getGates()[k].output].multByConst(e);
    }
    else if (_circuit.getGates()[k].gateType == SCALAR_ADD) {
      cout << "Scalar Add not implemented" << endl; // improve: 
      abort();
      // long scalar(_circuit.getGates()[k].input2);
      // FieldType e = _field->GetElement(scalar);
      // _gateShares[_circuit.getGates()[k].output] =
      //     _gateShares[_circuit.getGates()[k].input1] + e;
    } else if (_circuit.getGates()[k].gateType == MULT) {
      count++;
    }
  }
  return count;
}

template <class FieldType>
void ProtocolParty<FieldType>::
refresh(vector<FieldType>& xShares) {
  vector<int> NonTSet, TSet;
  int king = _disp.tAndNonTSet(TSet, NonTSet);
  int nShares = xShares.size();

  vector<FieldType> eShares(nShares);
  for (int i=0; i<nShares; i++) {
    eShares[i] = _refreshShares[_refreshOffset++]._value + xShares[i];
  }

  vector<byte> eSharesBytes(nShares * _fieldByteSize);
  _field->elementVectorToByteVector(eShares, eSharesBytes);
  vector<vector<byte>> recBufs, sendBufs;
  _comm.TToKing(king, eSharesBytes, recBufs); // rec from T only
  
  if (_myId == king) {          // king prepare oShares
    // improve: refactor
    vector<vector<FieldType>> oSharesAll(_N, vector<FieldType>(nShares, _zero));
    vector<FieldType> TShares(_T + 1), nonTShares(_N-_T, _zero);
    for (int i=0; i<nShares; i++) { // ith share
      vector<FieldType> eSharesAll(_N, _zero);
      for (int j=0; j<_T+1; j++) {
        int pt = TSet[j];
        TShares[j] =
          _field->bytesToElement(recBufs[pt].data() + i * _fieldByteSize);
        eSharesAll[pt] = TShares[j];
      }
      _protoTape.recordO2Share(eSharesAll);
      _mat_from_T.MatrixMult(TShares, nonTShares);
      for (int j=0; j<_N-_T-1; j++) {
        int pnt = NonTSet[j];
        if (_disp.isCorrupt(pnt)) {     // corr's oshare is its eshare
          oSharesAll[pnt][i] = nonTShares[j];
        } else {                // non-T's oshare is zero
          nonTShares[j] = _zero;
        }
      }
      _mat_to_T.MatrixMult(nonTShares, TShares);
      for (int j=0; j<_T+1; j++) {
        int pt = TSet[j];
        oSharesAll[pt][i] = TShares[j];
      }
    }
    _protoTape.recordDealtOShares(oSharesAll);
    
    sendBufs.resize(_N, vector<byte>(nShares * _fieldByteSize));
    for (int i=0; i<_N; i++) {
      if (_TMask[i]) {
        _field->elementVectorToByteVector(oSharesAll[i], sendBufs[i]);
      }
    }
  }

  vector<byte> oSharesBytes(nShares * _fieldByteSize);
  _comm.kingToT(king, oSharesBytes, sendBufs); // send to T only
  for (int i=0; i<nShares; i++) {
    FieldType oShare =
      _field->bytesToElement(oSharesBytes.data() + i * _fieldByteSize);
    xShares[i] = xShares[i] - oShare;
    _protoTape.recordRecOShare(oShare, king);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
prepareE2Shares(int layerStart, int layerEnd, int nMults,
                vector<FieldType>& e2Shares) {
  e2Shares.resize(nMults);
  int multCounter = 0;
  vector<FieldType> xShares(nMults);
  vector<FieldType> yShares(nMults);
  for (int k = layerStart; k < layerEnd; k++) {
    auto &gate = _circuit.getGates()[k];
    if (gate.gateType == MULT) {
      xShares[multCounter] = _gateShares[gate.input1]._value;
      yShares[multCounter] = _gateShares[gate.input2]._value;
      multCounter++;
    }
  }                                                  

  if (_disp.hasCorrupt()) {
    refresh(xShares);
  }

  for (int i=0; i<nMults; i++) {
    e2Shares[i] = xShares[i] * yShares[i] -
      _doubleShares[ 2*(_doubleOffset + i) + 1 ]._value;
  }
}

template <class FieldType>
int ProtocolParty<FieldType>::
processMult(int layerStart, int layerEnd, int nMults) {
  // collect e2 shares, refreshed if neccessary
  vector<FieldType> e2Shares;
  prepareE2Shares(layerStart, layerEnd, nMults, e2Shares);

  // -- batch mult
  vector<Share<FieldType>> e1Shares;
  // batchMultSpread(e2Shares, e1Shares);
  batchMultSingle(e2Shares, e1Shares);

  // -- store mult results to gate wire array
  int curOffset = _doubleOffset;
  int multCounter = 0;
  for (int k = layerStart; k < layerEnd; k++) {
    auto &gate = _circuit.getGates()[k];
    if (gate.gateType == MULT) {
      _gateShares[gate.output] =
        _doubleShares[2 * curOffset] + e1Shares[multCounter];
      curOffset++;
      multCounter++;
    }
  }

  _doubleOffset += nMults;
  return nMults;
}

template <class FieldType>
void ProtocolParty<FieldType>::
batchMultSingle(vector<FieldType>& e2Shares, vector<Share<FieldType>>& e1Shares) {
  vector<int> NonTSet, TSet;
  int king = _disp.tAndNonTSet(TSet, NonTSet);
  int nMults = e2Shares.size();
  
  vector<byte> sendBuf(nMults * _fieldByteSize);
  _field->elementVectorToByteVector(e2Shares, sendBuf);
  vector<vector<byte>> recBufs, sendBufs, relayBufs;
  _comm.allToOneStore(sendBuf, relayBufs, recBufs, king); // store relay
  vector<vector<FieldType>> relayVecs(_N);
  for (int partyId=0; partyId<_N; partyId++) {
    bytesToVector(relayBufs[partyId], relayVecs[partyId]);
  }
  _protoTape.recordE2Relay(relayVecs);

  if (_myId == king) {          // king reconstruct
    vector<FieldType> e2FullShare(_N), eSharesT(_T+1), eSharesNonT(_N-_T,_zero);
    vector<vector<FieldType>> eFullShares(_N, vector<FieldType>(nMults, _zero));
    for (int i = 0; i < nMults; i++) {
      for (int j = 0; j < _N; j++) {
        e2FullShare[j] =
          _field->bytesToElement(recBufs[j].data() + (i * _fieldByteSize));
      }
      _protoTape.recordE2Share(e2FullShare);
      eSharesNonT[_N-_T-1] = getSecret(e2FullShare); // secret = e, rest zero
      _mat_to_T.MatrixMult(eSharesNonT, eSharesT);
      for (int j=0; j < _T+1; j++) {
        int pt = TSet[j];
        eFullShares[pt][i] = eSharesT[j];
      }
    }
    _protoTape.recordDealtE1Shares(eFullShares);

    sendBufs.resize(_N, vector<byte>(nMults * _fieldByteSize, 0));
    for (int j=0; j < _T+1; j++) {
      int pt = TSet[j];
      _field->elementVectorToByteVector(eFullShares[pt], sendBufs[pt]);
    }
  }
  vector<byte> e1SharesByte(nMults * _fieldByteSize);
  _comm.kingToT(king, e1SharesByte, sendBufs); // send to T only
  // -- convert byte to elements
  vector<FieldType> e1Vals(nMults);
  for (int i = 0; i < nMults; i++) {
    e1Vals[i] =               // improve: add vec function
      _field->bytesToElement(e1SharesByte.data() + (i * _fieldByteSize));
  }
  unsigned int idxStart = _protoTape.getRecTapeLength();
  _protoTape.recordRecE1Shares(e1Vals, king);
  e1Shares.resize(nMults);
  for (int i=0; i<nMults; i++) {
    e1Shares[i] = _protoTape.getRecShare(idxStart+i);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
dotProd(int king, int groupSize,
        vector<Share<FieldType>>& a, vector<Share<FieldType>>& b,
        vector<mTranscript<FieldType>>& transcripts,
        vector<vector<vector<FieldType>>>& kingTranscripts,
        vector<vector<FieldType>>& relayTranscripts) {
  assert( a.size() == b.size() );
  int totalLength = a.size();
  int nMults = (totalLength + groupSize - 1) / groupSize;  
  
  // -- generate the 2t-sharings for xy - r
  // -- also output the transcripts of r2, e2
  vector<FieldType> e2Shares(nMults, _zero);
  vector<Share<FieldType>> e1Shares(nMults);
  for (int group_i = 0; group_i < nMults; group_i++) {
    int group_end = ((group_i + 1) * groupSize > totalLength) ?
      totalLength : (group_i + 1) * groupSize;
    for (int i = group_i * groupSize; i < group_end; i++) {
      e2Shares[group_i] += a[i]._value * b[i]._value;
    }
    e2Shares[group_i] = e2Shares[group_i] -
      _doubleShares[(_doubleOffset + group_i) * 2 + 1]._value;
  }

  batchMultSingle(e2Shares, e1Shares);

  // -- compute xy - r + [r]_t = t-sharing of xy
  // -- and output the transcripts
  transcripts.resize(nMults);
  for (int group_i = 0; group_i < nMults; group_i++) {
    Share<FieldType> rShare = _doubleShares[(_doubleOffset + group_i)*2];
    Share<FieldType> r2Share = _doubleShares[(_doubleOffset + group_i) * 2 + 1];
    transcripts[group_i].addTo(rShare, r2Share, e1Shares[group_i]._value,
                               e2Shares[group_i] + r2Share._value); // ab=e2+r2
  }
  if (_myId == king) {          // build kingTranscripts
    kingTranscripts.assign(2, vector<vector<FieldType>>(nMults));
    int dealtE1IdxStart = _protoTape.getDealtE1Length() - nMults;
    int recE2IdxStart = _protoTape.getRecE2Length() - nMults;
    for (int i=0; i<nMults; i++) {
      _protoTape.getDealtE1Share(dealtE1IdxStart+i, kingTranscripts[0][i]);
      _protoTape.getRecE2FullShare(recE2IdxStart+i, kingTranscripts[1][i]);
    }
  }
  if (_disp.isRelayer(king, _myId)) {
    relayTranscripts.assign(nMults, vector<FieldType>());
    int relayE2IdxStart = _protoTape.getRelayTapeLength() - nMults;
    for (int i=0; i<nMults; i++) {
      _protoTape.getE2RelayShare(relayE2IdxStart+i, relayTranscripts[i]);
    }
  }
  _doubleOffset += nMults;
}

template <class FieldType>
void ProtocolParty<FieldType>::
openTShares(int nShares, bool check,
            vector<FieldType> &shares, vector<FieldType> &clears) {
  vector<byte> sharesByte(nShares * _fieldByteSize);
  _field->elementVectorToByteVector(shares, sharesByte);
  vector<vector<byte>> sendBufs(_N, sharesByte),
    recBufs(_N, vector<byte>(nShares * _fieldByteSize));
  _comm.allToAll(sendBufs, recBufs, true); // use relay

  vector<FieldType> x1(_N);
  for (int i = 0; i < nShares; i++) {
    for (int j = 0; j < _N; j++) {
      x1[j] = _field->bytesToElement(recBufs[j].data() + (i*_fieldByteSize));
    }
    if (check) {
      clears[i] = getSecretCheck(x1);
    } else {
      clears[i] = getSecret(x1);
    }
  }
}

template <class FieldType>
FieldType ProtocolParty<FieldType>::challenge() {
  // TODO: implement proper challenge w/ *broadcast*
  vector<FieldType> lambdaShare(1, _coinShares[_coinOffset++]._value);
  vector<FieldType> lambda(1);
  openTShares(1, true, lambdaShare, lambda);
  vector<vector<byte>> sendEmpty(_N), recEmpty(_N); // sync for security
  _comm.allToAll(sendEmpty, recEmpty);
  return lambda[0];
}

template <class FieldType> void ProtocolParty<FieldType>::outputPhase() {
  vector<FieldType> outShareAll(_N); // vector for the shares of my outputs
  vector<vector<FieldType>> sendElems(_N);
  vector<vector<byte>> sendBufs(_N);
  for (int k = _nGates - _nOutputs; k < _nGates; k++) {
    assert(_circuit.getGates()[k].gateType == OUTPUT);
    // send to party (which need this gate) your share for this gate
    int outParty = _circuit.getGates()[k].party;
    int outWire = _circuit.getGates()[k].input1;
    sendElems[outParty].push_back(_gateShares[outWire]._value);
  }

  for (int i = 0; i < _N; i++) {
    sendBufs[i].resize(sendElems[i].size() * _fieldByteSize);
    _field->elementVectorToByteVector(sendElems[i], sendBufs[i]);
  }

  vector<vector<byte>> recBufs;
  for (int i=0; i<_N; i++) {
    if (!_disp.isCorrupt(i)) {
      _comm.allToOne(sendBufs[i], recBufs, i, true); // use relay
    }
  }

  int counter = 0;
  for (int k = _nGates - _nOutputs; k < _nGates; k++) {
    assert(_circuit.getGates()[k].gateType == OUTPUT);
    int outParty = _circuit.getGates()[k].party;
    int outWire = _circuit.getGates()[k].input1;
    if (outParty == _myId) {
      for (int i = 0; i < _N; i++) {
        outShareAll[i] =
          _field->bytesToElement(recBufs[i].data() + counter*_fieldByteSize);
      }
      _protoTape.recordOutShare(outShareAll);
      if (flag_print_output)
        cout << "the result for " << outWire
             << " is : " << _field->elementToString(getSecret(outShareAll)) << '\n';
      counter++;
    }
  }
}

template <class FieldType> ProtocolParty<FieldType>::~ProtocolParty() {
  delete _field;
}
#endif /* PROTOCOLPARTY_H_ */
