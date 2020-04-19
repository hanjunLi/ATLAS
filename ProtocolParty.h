// Hanjun Li <lihanjun1212@gmail.com>

#ifndef PROTOCOLPARTY_H_
#define PROTOCOLPARTY_H_

#include <stdlib.h>
// #include "ProtocolTimer.h"
// #include <chrono>
// #include <libscapi/include/infra/Measurement.hpp>
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
#include <libscapi/include/primitives/Mersenne.hpp>
#include <libscapi/include/cryptoInfra/Protocol.hpp>
#include <libscapi/include/comm/MPCCommunication.hpp>
#include <libscapi/include/primitives/Matrix.hpp>
#include <libscapi/include/circuits/ArithmeticCircuit.hpp>

#define flag_print_output false

using namespace std;
// using namespace std::chrono;

// TODO: implement localization and analyzeSharing protocols
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
  int _fieldByteSize, _segSize, _nSegs, _authKeySize;
  
  // -- communication channels and dispute sets
  Communication _comm;
  Dispute _disp;
  
  // -- various matrices set in constructor
  vector<FieldType> _alpha_k;    // K distinct non-zero field elements
  vector<FieldType> _alpha_2km1; // 2K-1 distinct non-zero field element
  HIM<FieldType> _mat_use_n;     // all _N points to 0
  vector<HIM<FieldType>> _mat_use_n_twist_vec; // all _N points (and 0) to i
  HIM<FieldType> _mat_check_tp1; // first t+1 points to rest N-t-1
  HIM<FieldType> _mat_check_km1; // first k points to rest 2k-1 - k
  VDMTranspose<FieldType> _mat_vand_trans; // for extracting random shares
  vector<FieldType> _buffer;    // for interp to a sinlge element

  // -- local storage: global
  vector<FieldType> _gateShares;         // my share of each gate
  vector<vector<FieldType>> _deltShares; // shares I delt: accumulated
  vector<vector<FieldType>> _kingShares; // shares I reconstructed as king
  vector<FieldType> _recShares;          // shares I received
  vector<int> _dealerIds;                // dealers for recShares
  vector<vector<FieldType>> _keys;       // [i]-auth keys for nonDips (party j)
  vector<vector<vector<FieldType>>> _keyShares; // [i]-auth key shares received from j of k
  vector<vector<FieldType>> _recTags;   // tags from i of jth share
  vector<vector<FieldType>> _sentTags;  // tags for i about jth share
  int _verifyOffset = 0;     // upto which recShares are verified
  int _authOffset = 0;       // upto which recShares are authenticated

  // -- local storage: per seg
  vector<FieldType> _coinShares;    // for public random coins
  int _coinOffset = 0;              // next share for coin
  vector<FieldType> _padShares;     // for padding transcripts
  int _padOffset = 0;               // next share for padding
  vector<FieldType> _refreshShares; // used for refresh only
  int _refreshOffset = 0;           // next share for refresh
  vector<FieldType> _doubleShares;  // used for DN mult only
  vector<int> _doubleKingIds;       // 1 king per pair
  int _doubleOffset = 0;            // next random double share to use

  vector<int> _e1ShareIdx;       // e1 transcripts (ptrs)
  vector<FieldType> _recOShares; // refresh transcripts (direct)
  vector<vector<FieldType>> _deltOShares;

  HIM<FieldType> _mat_to_T;    // for calculating share of T set
  HIM<FieldType> _mat_from_T;  // for calculating share of nonT set
  vector<HIM<FieldType>> _mat_twist_vec;
  bool _newDisp = true;         // flag new Disp pair added
  vector<bool> _TMask;

  // -- simple local helper functions
  __inline__
  void segMsg(int segStart, int segEnd, string msg) {
    // cout << "[" << _myId << "] Seg: "
    //      << segStart << " - " << segEnd << ": "
    //      << msg  << endl;
  }

  __inline__
  FieldType getSecret(vector<FieldType> fullShares) {
    assert(fullShares.size() == _N);
    _mat_use_n.MatrixMult(fullShares, _buffer);
    return _buffer[0];
  }

  __inline__
  void getPadShares(int nRands, vector<FieldType> &result) {
      result.assign(_padShares.begin() + _padOffset,
                    _padShares.begin() + _padOffset + nRands);
      _padOffset += nRands;
  }

  __inline__
  bool checkTShares(vector<FieldType> &fullShares) {
    vector<FieldType> restTShares(_N-1-_T); // = T
    vector<FieldType> firstTp1(fullShares.begin(), fullShares.begin() + _T+1);
    _mat_check_tp1.MatrixMult(firstTp1, restTShares);
    for (int i = 0; i < _N-_T-1; i++) {
      if (restTShares[i] != fullShares[_T+1+i]) {
        return false;
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
  void recordDeltShares(vector<vector<FieldType>>& fullSharesVec,
                        vector<vector<FieldType>>& deltShares) {
    int nShares = fullSharesVec[0].size();
    vector<FieldType> fullShares(_N);
    for (int i=0; i<nShares; i++) {
      for (int j=0; j<_N; j++) {
        fullShares[j] = fullSharesVec[j][i];
      }
      deltShares.push_back(fullShares);
    }
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
  void matForRefresh(int pid, HIM<FieldType>& mat_from_T, HIM<FieldType>& mat_to_T) {
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
  int makeAuthBatch(vector<vector<FieldType>>& batchToAuth, vector<int>& batchDealer) {
    auto authStart = _recShares.begin() + _authOffset;
    vector<FieldType> toAuth(authStart, _recShares.end());
    auto authDealerStart = _dealerIds.begin() + _authOffset;
    vector<int> toAuthDealer(authDealerStart, _dealerIds.end());
  
    int nToAuth = _recShares.size() - _authOffset;
    // -- compute [key] * [dealer batch] shares
    vector<vector<FieldType>> curBatch(_N);
    for (int i=0; i<nToAuth; i++) {
      int dealer = toAuthDealer[i];
      curBatch[dealer].push_back(toAuth[i]);
      if (curBatch[dealer].size() == _authKeySize) { // store current batch
        batchToAuth.push_back(curBatch[dealer]);
        curBatch[dealer].clear();
        batchDealer.push_back(dealer);
      }
    }
    for (int i=0; i<_N; i++) {
      if (curBatch[i].size() > 0) {
        curBatch[i].resize(_authKeySize);
        batchToAuth.push_back(curBatch[i]);
        batchDealer.push_back(i);
      }
    }
    return batchToAuth.size();
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
  void delinearization(int segStart, int segEnd, vector<FieldType>& aShares,
                       vector<FieldType>& bShares, vector<FieldType>& transcript);
  
  // -- main protocol functions
  void run() override;
  bool runOffline(int segStart, int segEnd);
  bool prepareSeg();
  void makeRandShares(int nRands, vector<FieldType> &randShares);
  void makeRandDoubleShares(int nRands, vector<FieldType> &randDoubleShares, vector<int>& kingIds);
  bool runOnline(int segStart, int segEnd);
  void inputPhase();
  void computationPhase(int segStart, int segEnd);
  int processMult(int layerStart, int layerEnd, int nMults);
  void batchMultSingle(vector<FieldType>& e2Shares, vector<FieldType>& e1Shares);
  void refresh(vector<FieldType>& xShares);
  bool verificationPhase(int segStart, int segEnd);
  void outputPhase();

  // -- verification protocols
  bool checkRefresh(int segStart, int segEnd);
  int dimensionReduction(vector<FieldType>& aShares, vector<FieldType>& bShares,
                         vector<FieldType>& transcript);
  void compress(vector<FieldType>& aShares, vector<FieldType>& bShares,
                vector<vector<FieldType>>& transcripts, int groupSize,
                vector<FieldType>& aSharesNew, vector<FieldType>& bSharesNew,
                vector<FieldType>& transcriptNew);
  void randomization(vector<FieldType>& aShares, vector<FieldType>& bShares,
                     vector<FieldType>& transcript, vector<FieldType>& ultimateTranscript);
  bool checkSingleMult(vector<FieldType>& ultimateTranscript);
  bool verifySharing();
  bool tag();
  void keyDistribution(vector<FieldType>& keys, vector<vector<FieldType>>& keyShares);
  void singleTagComp();
  bool checkKey();
  bool checkTag();
  bool verifyOutput();   // TODO: [ask] public output / private output
  FieldType challenge();

  // -- helpers
  void buildPolyVecInd(vector<FieldType>& aShares, vector<FieldType>& bShares,
                       vector<vector<FieldType>>& trancripts, int groupSize);
  void buildPolyVecIndWorker(vector<FieldType>& aShares, vector<FieldType>& bShares,
                             vector< vector<FieldType> >& ASharesEval,
                             vector< vector<FieldType> >& BSharesEval,
                             int groupSize, int threadId);
  void multVec(vector<FieldType>& aShares, vector<FieldType>& bShares,
               vector<vector<FieldType>>& transcripts, int groupSize);
  void openTShares(int nShares, bool check, vector<FieldType> &Shares, vector<FieldType> &clears);
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
    cout << "don't use other fields yet" << endl;
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
  _authKeySize = _segSize / (_T + 1) / 10; // 10 is flexible
  
  // -- communication channels
  string partiesFile =
      getParser().getValueByKey(arguments, "partiesFile");
  _comm.reset(_N, _myId, _numThreads, partiesFile, &_disp);
  _disp.reset(_N);
  _sentTags.resize(_N);
  _recTags.resize(_N);
  // -- global storage
  _gateShares.resize(_nGates - _nOutputs);
  _buffer.resize(1);
  initializeMat();
}

template <class FieldType> void ProtocolParty<FieldType>::run() {
  int segStart = 0;
  int segEnd = 0;
  const auto& gates = _circuit.getGates();
  
  for (int i = 0; i < _iterations; i++) {
    for (int curSeg = 0; curSeg < _nSegs; curSeg++){
      // if (curSeg == 2) {        // TODO: delete this, test only
      //   _disp.addDispPairs(0, 1);
      //   _disp.addDispPairs(2, 3);
      //   _disp.addCorrParty(0);
      // }
      
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
      while (!runOffline(segStart, segEnd)) {
        segMsg(segStart, segEnd, "failed runOffline");
      }
      segMsg(segStart, segEnd, "finished runOffline");
      while (!runOnline(segStart, segEnd)) {
        segMsg(segStart, segEnd, "failed runOnline");
      }
      segMsg(segStart, segEnd, "finished runOnline");
      _disp.setNewDisp(false);
    }

    
  }
}

  
template <class FieldType>
bool ProtocolParty<FieldType>::runOffline(int segStart, int segEnd) {
  if (_disp.getNewDisp()) {     // update matrices w/ new disp
    // -- prepare mat for normal sharing
    matForRefresh(_myId, _mat_from_T, _mat_to_T);
    // -- prepare mat for twisted sharing
    matForTwist(_myId, _mat_twist_vec);
    _disp.tMaskP(_myId, _TMask);
    
  }
  
  return prepareSeg();
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
bool ProtocolParty<FieldType>::checkRefresh(int segStart, int segEnd) {
  if (!_disp.hasCorrupt()) {
    return true;
  }
  FieldType x, xp, r, e, o;
  // -- pad with a random transcript (x, xp, r, e, o)
  vector<FieldType> pad(1), transcript(5);
  getPadShares(1, pad);         // TODO: [ask] why pad in K?
  transcript[0] = pad[0];
  refresh(pad);
  transcript[1] = pad[0];
  transcript[2] = _refreshShares[_refreshOffset]; // TODO: [ask] why T?
  transcript[3] = x + r;
  transcript[4] = *(_recOShares.rbegin()); // the one in last refresh
  // -- combine transcripts with public coin
  FieldType lambda = challenge();
  FieldType lambdaI = lambda;
  int multCount = 0;
  for (int k = segStart; k < segEnd; k++) {
    auto& gate = _circuit.getGates()[k];
    if (gate.gateType == MULT) {
      int inWire = gate.input1;
      transcript[0] += _gateShares[inWire] * lambdaI;
      transcript[1] +=
        (_gateShares[inWire] -_recOShares[multCount]) * lambdaI;
      transcript[2] += _refreshShares[multCount] * lambdaI;
      transcript[3] +=
        (_gateShares[inWire] + _refreshShares[multCount]) * lambdaI;
      transcript[4] += _recOShares[multCount] * lambdaI;
      multCount++;
      lambdaI *= lambda;
    }
  }
  // TODO: proper check w/ analyzeSharing
  
  vector<FieldType> check(1, transcript[1]);
  openTShares(1, true, check, _buffer);
  return true;
}


template <class FieldType>
bool ProtocolParty<FieldType>::verificationPhase(int segStart, int segEnd) {
  if (!checkRefresh(segStart, segEnd)) {
    return false;
  }
  vector<FieldType> aShares, bShares, transcript; // [r, r2, e, e2, c]
  delinearization(segStart, segEnd, aShares, bShares, transcript);
  while(dimensionReduction(aShares, bShares, transcript) >= _K);
  vector<FieldType> ultimateTranscript(7);
  randomization(aShares, bShares, transcript, ultimateTranscript);
  return (checkSingleMult(ultimateTranscript));
}

template <class FieldType>
void ProtocolParty<FieldType>::
delinearization(int segStart, int segEnd, vector<FieldType>& aShares,
                 vector<FieldType>& bShares, vector<FieldType>& transcript){
  aShares.resize(_segSize);
  bShares.resize(_segSize);
  transcript.resize(5, _zero);
  FieldType lambda = challenge();
  FieldType lambdaI = lambda;
  int multCount = 0;
  for (int k = segStart; k < segEnd; k++) {
    auto& gate = _circuit.getGates()[k];
    if(gate.gateType == MULT) {
      aShares[multCount] = _gateShares[gate.input1];
      if (_disp.hasCorrupt()) {
        aShares[multCount] = aShares[multCount] - _recOShares[multCount];
      }
      bShares[multCount] = _gateShares[gate.input2] * lambdaI;
      transcript[0] += _doubleShares[multCount*2] * lambdaI;
      transcript[1] += _doubleShares[multCount*2+1] * lambdaI;
      transcript[2] += _recShares[_e1ShareIdx[multCount]] * lambdaI;
      transcript[3] += (aShares[multCount] * bShares[multCount] -
                        _doubleShares[multCount*2+1]) * lambdaI;
      transcript[4] += _gateShares[gate.output] * lambdaI;      
      multCount++;
      lambdaI *= lambda;
    }
  }
  aShares.resize(multCount);
  bShares.resize(multCount);
}

template <class FieldType>
int ProtocolParty<FieldType>::
dimensionReduction(vector<FieldType>& aShares, vector<FieldType>& bShares,
                   vector<FieldType>& transcript){
  
  // improve: comment on transcript format
  int totalLength = aShares.size();
  if (totalLength < _K) { // no need to compress
    return totalLength;
  }
  // -- divide into K groups and pad short groups
  int groupSize = (totalLength + _K -1 )/ _K;
  totalLength = groupSize * _K;
  aShares.resize(totalLength, _zero);
  bShares.resize(totalLength, _zero);

  // -- one DN mult for each group i to build dShares 0 .. k-2
  vector<vector<FieldType>> transcripts(5);
  multVec(aShares, bShares, transcripts, groupSize);
  
  // build dShares k-1: c - d_0 - ... - d_{k-2}
  transcripts[0][_K-1] = transcript[0];
  transcripts[1][_K-1] = transcript[1];
  transcripts[2][_K-1] = transcript[2];
  transcripts[3][_K-1] = transcript[3];
  transcripts[4][_K-1] = transcript[4];
  for (int i = 0; i< _K - 1; i++) {
    transcripts[0][_K-1] = transcripts[0][_K-1] - transcripts[0][i];
    transcripts[1][_K-1] = transcripts[1][_K-1] - transcripts[1][i];
    transcripts[2][_K-1] = transcripts[2][_K-1] - transcripts[2][i];
    transcripts[3][_K-1] = transcripts[3][_K-1] - transcripts[3][i];
    transcripts[4][_K-1] = transcripts[4][_K-1] - transcripts[4][i];    
  }

  vector<FieldType> aSharesNew(groupSize);
  vector<FieldType> bSharesNew(groupSize);
  compress(aShares, bShares, transcripts, groupSize,
           aSharesNew, bSharesNew, transcript);
    
  aShares = aSharesNew;
  bShares = bSharesNew;
  return groupSize;
}

template <class FieldType>
void ProtocolParty<FieldType>::
compress(vector<FieldType>& aShares, vector<FieldType>& bShares,
         vector<vector<FieldType>>& transcripts, int groupSize,
         vector<FieldType>& aSharesNew, vector<FieldType>& bSharesNew,
         vector<FieldType>& transcriptNew) {
  // improve: comment transcripts format
  int totalLength = groupSize * _K;
  
  // -- compress k groups of vector product  
  buildPolyVecInd(aShares, bShares, transcripts, groupSize);
  
  FieldType lambda = challenge();
  vector<FieldType> beta_lambda(1);
  beta_lambda[0] = lambda;  
  HIM<FieldType> matrix_for_k_lambda(1, _K, _field);
  matrix_for_k_lambda.InitHIMByVectors(_alpha_k, beta_lambda);
  HIM<FieldType> matrix_for_2k_lambda(1, 2*_K-1, _field);
  matrix_for_2k_lambda.InitHIMByVectors(_alpha_2km1, beta_lambda);

  // improve: refactor
  vector<FieldType> ySharesA;   // tmp vector
  vector<FieldType> ySharesB;   // tmp vector
  for (int i=0; i<groupSize; i++) {
    ySharesA.clear();
    ySharesB.clear();
    for (int j=i; j<totalLength; j+=groupSize) {
      ySharesA.push_back(aShares[j]);
      ySharesB.push_back(bShares[j]);
    }      
    matrix_for_k_lambda.MatrixMult(ySharesA, _buffer);
    aSharesNew[i] = _buffer[0];
    matrix_for_k_lambda.MatrixMult(ySharesB, _buffer);
    bSharesNew[i] = _buffer[0];
  }

  matrix_for_2k_lambda.MatrixMult(transcripts[0], _buffer);
  transcriptNew[0] = _buffer[0];
  matrix_for_2k_lambda.MatrixMult(transcripts[1], _buffer);
  transcriptNew[1] = _buffer[0];
  matrix_for_2k_lambda.MatrixMult(transcripts[2], _buffer);
  transcriptNew[2] = _buffer[0];
  matrix_for_2k_lambda.MatrixMult(transcripts[3], _buffer);
  transcriptNew[3] = _buffer[0];
  matrix_for_2k_lambda.MatrixMult(transcripts[4], _buffer);
  transcriptNew[4] = _buffer[0];
  return;
}

template <class FieldType>
void ProtocolParty<FieldType>::
randomization(vector<FieldType>& aShares, vector<FieldType>& bShares,
              vector<FieldType>& transcript,
              vector<FieldType>& ultimateTranscript) {
  assert(aShares.size() < _K);
  assert(bShares.size() < _K);
  aShares.resize(_K, _zero);
  bShares.resize(_K, _zero);

  vector<FieldType>  randPair(2);
  getPadShares(2, randPair);
  aShares[_K-1] = randPair[0];
  bShares[_K-1] = randPair[1];

  vector<vector<FieldType>> transcripts(5);
  multVec(aShares, bShares, transcripts, 1);
  
  // build dShares k-1: c - d_0 - ... - d_{k-2}
  transcripts[0][_K-1] += transcript[0];
  transcripts[1][_K-1] += transcript[1];
  transcripts[2][_K-1] += transcript[2];
  transcripts[3][_K-1] += transcript[3];
  transcripts[4][_K-1] += transcript[4];
  for (int i = 0; i< _K - 1; i++) {
    transcripts[0][_K-1] = transcripts[0][_K-1] - transcripts[0][i];
    transcripts[1][_K-1] = transcripts[1][_K-1] - transcripts[1][i];
    transcripts[2][_K-1] = transcripts[2][_K-1] - transcripts[2][i];
    transcripts[3][_K-1] = transcripts[3][_K-1] - transcripts[3][i];
    transcripts[4][_K-1] = transcripts[4][_K-1] - transcripts[4][i];    
  }

  vector<FieldType> aSharesNew(1);
  vector<FieldType> bSharesNew(1);
  compress(aShares, bShares, transcripts, 1,
           aSharesNew, bSharesNew, transcript);
  ultimateTranscript.resize(7);
  ultimateTranscript[0] = aSharesNew[0];
  ultimateTranscript[1] = bSharesNew[0];
  copy(transcript.begin(), transcript.end(), ultimateTranscript.begin() + 2);
  return;
}

template <class FieldType>
bool ProtocolParty<FieldType>::
checkSingleMult(vector<FieldType> &ultimateTranscript) {
  // TODO: impl proper check-single-mult protocol
  // 1. check procedure
  // 2. update dispute, corrupt sets
  
  // -- currently only check t-shares a, b, r1, e1, c
  //    satisfying a * b = c, r1 + e1 = c
  vector<FieldType> sharesToCheck
    {ultimateTranscript[0] * ultimateTranscript[1], // a * b
     ultimateTranscript[2], ultimateTranscript[4],  // r1, e1
     ultimateTranscript[6]};                        // c

  vector<FieldType> sharesOpen(4);
  openTShares(4, false, sharesToCheck, sharesOpen);

  assert(sharesOpen[0] == sharesOpen[3]);
  assert(sharesOpen[1] + sharesOpen[2] == sharesOpen[3]);
  return true;
}

template <class FieldType>
void ProtocolParty<FieldType>::
buildPolyVecIndWorker(vector<FieldType>& aShares, vector<FieldType>& bShares,
                      vector< vector<FieldType> >& ASharesEval,
                      vector< vector<FieldType> >& BSharesEval,
                      int groupSize, int threadId) {
  int totalLength = aShares.size();
  vector<FieldType> ySharesA;   // tmp vector
  vector<FieldType> ySharesB;   // tmp vector
  
  for (int i=threadId; i<groupSize; i+= _numThreads) {
    ySharesA.clear();
    ySharesB.clear();
    for (int j=i; j<totalLength; j+=groupSize) {
      ySharesA.push_back(aShares[j]);
      ySharesB.push_back(bShares[j]);
    }
    _mat_check_km1.MatrixMult(ySharesA, ASharesEval[i]);
    _mat_check_km1.MatrixMult(ySharesB, BSharesEval[i]);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
buildPolyVecInd(vector<FieldType>& aShares, vector<FieldType>& bShares,
                vector<vector<FieldType>>& transcripts, int groupSize) {

  // batchDegree should be _K
  int batchDegree = (aShares.size() + groupSize - 1) / groupSize;
  
  // interpolate A[i], B[i] and evaluate at A(k), ..., A(2k-2), B(k), ..., B(2k-2)
  vector< vector<FieldType> > ASharesEval(groupSize, vector<FieldType>(batchDegree-1));
  vector< vector<FieldType> > BSharesEval(groupSize, vector<FieldType>(batchDegree-1));
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
  vector<FieldType> ASharesEvalFlat(groupSize * (batchDegree - 1));
  vector<FieldType> BSharesEvalFlat(groupSize * (batchDegree - 1));
  for (int i=0; i<batchDegree-1; i++) {
    // point i
    for (int j=0; j<groupSize; j++) {
      // poly j
      ASharesEvalFlat[i * groupSize + j] = ASharesEval[j][i];
      BSharesEvalFlat[i * groupSize + j] = BSharesEval[j][i];
    }
  }

  vector<vector<FieldType>> newTranscripts(5);
  multVec(ASharesEvalFlat, BSharesEvalFlat, newTranscripts, groupSize);

  // append to passed-in dShares
  transcripts[0].insert(transcripts[0].end(), newTranscripts[0].begin(), newTranscripts[0].end());
  transcripts[1].insert(transcripts[1].end(), newTranscripts[1].begin(), newTranscripts[1].end());
  transcripts[2].insert(transcripts[2].end(), newTranscripts[2].begin(), newTranscripts[2].end());
  transcripts[3].insert(transcripts[3].end(), newTranscripts[3].begin(), newTranscripts[3].end());
  transcripts[4].insert(transcripts[4].end(), newTranscripts[4].begin(), newTranscripts[4].end());
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
  
  vector<vector<FieldType>> inputSharesAll;
  makeTShares(secrets, inputSharesAll);
  recordDeltShares(inputSharesAll, _deltShares);
  vector<vector<byte>> sendBufs, recBufs(_N);
  sendBufs.resize(_N, vector<byte>(sizes[_myId] * _fieldByteSize));
  for (int i = 0; i < _N; i++) {
    _field->elementVectorToByteVector(inputSharesAll[i], sendBufs[i]);
    recBufs[i].resize(sizes[i] * _fieldByteSize, 0);
  }
  for (int i=0; i<_N; i++) {
    if (!_disp.isCorrupt(i)) {
      _comm.oneToAll(recBufs[i], sendBufs, i, false);
    }
  }
  // -- convert received bytes to shares
  vector<vector<FieldType>> inputShares(_N);
  for (int i = 0; i < _N; i++) {
    inputShares[i].resize(sizes[i], _zero);
    for (int j = 0; j < sizes[i]; j++) { // receive actual shares
      inputShares[i][j] =
        _field->bytesToElement(recBufs[i].data() + (j * _fieldByteSize));
    }      
    // vv record received shares and their dealer
    vector<int> idVec(sizes[i], i);
    _recShares.insert(_recShares.end(), inputShares[i].begin(), inputShares[i].end());
    _dealerIds.insert(_dealerIds.end(), idVec.begin(), idVec.end());
  }
  // -- store shares in input gate order
  vector<int> counters(_N, 0);
  for (int k = 0; k < _nInputs; k++) {
    assert(_circuit.getGates()[k].gateType == INPUT);
    int inParty = _circuit.getGates()[k].party;
    int inWire = _circuit.getGates()[k].output;
    _gateShares[inWire] = inputShares[inParty][counters[inParty]++];
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
makeRandShares(int nRands, vector<FieldType> &randShares) {
  int nBuckets = (nRands / (_N - _T)) + 1;
  randShares.resize(nBuckets * (_N - _T));
  
  vector<FieldType> secrets(nBuckets);
  for (int i=0; i<nBuckets; i++) {
    secrets[i] = _field->Random();
  }
  vector<vector<FieldType>> randSharesAll(_N, vector<FieldType>(nBuckets));
  makeTShares(secrets, randSharesAll);
  
  recordDeltShares(randSharesAll, _deltShares);
  vector<vector<byte>> sendBufs(_N, vector<byte>(nBuckets*_fieldByteSize));
  for (int i=0; i<_N; i++) {
    _field->elementVectorToByteVector(randSharesAll[i], sendBufs[i]);
  }
  vector<vector<byte>> recBufs(_N, vector<byte>(nBuckets*_fieldByteSize, 0));
  _comm.allToAll(sendBufs, recBufs, false); // no relay
  
  int count = 0;
  vector<FieldType> bufferIn(_N), bufferOut(_N);
  for (int i=0; i<nBuckets; i++) { // improve: switch loop order
    for (int j=0; j<_N; j++) {
      bufferIn[j] =
        _field->bytesToElement(recBufs[j].data() + i * _fieldByteSize);
      // record received share and its dealer id
      _recShares.push_back(bufferIn[j]);
      _dealerIds.push_back(j);
    }
    _mat_vand_trans.MatrixMult(bufferIn, bufferOut, _N - _T);
    for (int j=0; j<_N - _T; j++) {
      randShares[count++] = bufferOut[j];
    }
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
makeRandDoubleShares(int nRands, vector<FieldType> &randSharePairs, vector<int>& kingIds) {
  int nBuckets = (nRands / (_N - _T)) + 1;
  randSharePairs.resize(nBuckets * (_N - _T) * 2);
  kingIds.resize(nBuckets * (_N - _T));

  vector<FieldType> secrets;
  vector<vector<FieldType>> dSharesAll, sharesAll, randSharePairsAll;
  make2TShares(nBuckets, secrets, dSharesAll);
  makeTShares(secrets, sharesAll);
  recordDeltShares(sharesAll, _deltShares);
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
  int kingId = 0;
  vector<FieldType> bufferIn(_N), bufferIn2(_N), bufferOut(_N), bufferOut2(_N);
  for (int i=0; i<nBuckets; i++) {
    for (int j=0; j<_N; j++) {
      bufferIn[j] =
        _field->bytesToElement(recBufs[j].data() + 2*i*_fieldByteSize);
      bufferIn2[j] =
        _field->bytesToElement(recBufs[j].data() + (2*i+1)*_fieldByteSize);
      _recShares.push_back(bufferIn[j]);
      _dealerIds.push_back(j);
    }
    _mat_vand_trans.MatrixMult(bufferIn, bufferOut, _N - _T);
    _mat_vand_trans.MatrixMult(bufferIn2, bufferOut2, _N - _T);    
    for (int j=0; j<_N-_T; j++) {
      if (!_disp.isCorrupt(j)) {
        randSharePairs[count*2] = bufferOut[j];
        randSharePairs[count*2+1] = bufferOut2[j];
        while (_disp.isCorrupt(kingId)) {
          kingId = (kingId+1) % _N;
        }
        kingIds[count] = kingId;
        count++;
        kingId = (kingId+1) % _N;
      }
    }
  }
}

template <class FieldType>
bool ProtocolParty<FieldType>::
verifySharing() {
  vector<FieldType> sigma(1);
  getPadShares(1, sigma);
  FieldType lambda = challenge();
  FieldType lambdaI = lambda;
  auto verifyStart = _recShares.begin() + _verifyOffset;
  vector<FieldType> toVerify(verifyStart, _recShares.end());
  for(auto share : toVerify) {
    sigma[0] += share * lambdaI;
    lambdaI *= lambda;
  }
  // TODO: impl proper check w/ analyze-sharing
  vector<FieldType> result(1);
  openTShares(1, true, sigma, result);
  _verifyOffset = _recShares.size();
  return true;
}

template <class FieldType>
bool ProtocolParty<FieldType>::
tag() {
  if (_disp.getNewDisp()) {     // improve: new disp flag per party
    _keyShares.resize(_authKeySize);
    _keys.resize(_authKeySize);
    // improve: batch this
    for (int i=0; i<_authKeySize; i++) {
      keyDistribution(_keys[i], _keyShares[i]);
    }
    if (!checkKey()) {
      return false;
    }
    // cout << "[" << _myId << "] ran keyDistribution" << endl;
  }
  
  singleTagComp();
  if (!checkTag()) {
    return false;
  }
  _authOffset = _recShares.size();
  return true;
}

template <class FieldType>
void ProtocolParty<FieldType>::
keyDistribution(vector<FieldType>& keys, vector<vector<FieldType>>& keyShares) {
  vector<int> DispSet, NonDispSet;
  int nDisp = _disp.dispAndNonDispSet(_myId, DispSet, NonDispSet);
  int nNonDisp = _N - nDisp;
  vector<vector<byte>> sendBufs(_N, vector<byte>(_N*_fieldByteSize, 0));
  vector<vector<FieldType>> keySharesAll(_N, vector<FieldType>(_N, _zero));
  keys.clear();
  keys.resize(_N, _zero);
  keyShares.clear();
  keyShares.resize(_N, vector<FieldType>(_N, _zero));
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
      keyShares[i][pid] =
        _field->bytesToElement(recBufs[i].data() + pid * _fieldByteSize);
    }
  }
}

template <class FieldType>
bool ProtocolParty<FieldType>::checkKey() {
  vector<FieldType> padKeys;
  vector<vector<FieldType>> padKeyShares;
  keyDistribution(padKeys, padKeyShares);
  FieldType lambda = challenge();
  for (int i=0; i<_N; i++) {
    for (int j=0; j<_N; j++) {
      FieldType lambdaI = lambda;
      for (int k = 0; k < _authKeySize; k++) {
        padKeyShares[i][j] += _keyShares[k][i][j] * lambdaI;
        lambdaI *= lambda;
      }
    }
  }

  vector<byte> sendBuf;
  for (int i=0; i<_N; i++) {
    vector<byte> tmp(_N * _fieldByteSize);
    _field->elementVectorToByteVector(padKeyShares[i], tmp);
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
  // TODO: implement proper check w/ fault localization
  return true;
}

// improve: refactor
template <class FieldType>
void ProtocolParty<FieldType>::
singleTagComp() {
  vector<vector<FieldType>> authBatchs;
  vector<int> authDealers;
  int nAuthBatchs = makeAuthBatch(authBatchs, authDealers);
  // -- compute and distribute v and o shares per batch
  vector<vector<FieldType>>
    tagSharesAll(_N, vector<FieldType>(nAuthBatchs * _N, _zero)),
    zeroSharesAll(_N, vector<FieldType>(nAuthBatchs * _N, _zero));
  for (int i=0; i<nAuthBatchs; i++) {
    int dealer = authDealers[i];
    for (int j=0; j<_N; j++) {  // verify dealer's share for j
      if (_disp.isDisp(_myId, j) || _disp.isDisp(dealer, j)) {
        continue;               // leave this slot to all zero
      }
      // -- generate 2t twist share of tag(i,j)
      vector<int> inIds = allNButOne(j);
      vector<FieldType> inShares(_N, _zero);
      for (int k=0; k<_N-1; k++) {
        int pid = inIds[k];
        if (_disp.isDisp(dealer, pid)) {
          continue;             // disp(dealer) shares are zero
        }
        inShares[k] = _field->Random();
        tagSharesAll[pid][ i*_N+j ] = inShares[k];
      }
      _mat_use_n_twist_vec[j].MatrixMult(inShares, _buffer);
      // _sentTags[j][_authOffset+i] = _buffer[0];
      _sentTags[j].push_back(_buffer[0]);
      // -- generate 2t twist share of zero
      inIds = allNButOne(dealer);
      inShares.clear();
      inShares.resize(_N, _zero);
      for (int k=0; k<_N-1; k++) {
        int pid = inIds[k];
        if (_disp.isDisp(dealer, pid) || pid == j) {
          continue;             // j's share is fixed to zero
        }
        inShares[k] = _field->Random();
        zeroSharesAll[pid][i*_N+j] = inShares[k];
      }
      _mat_use_n_twist_vec[dealer].MatrixMult(inShares, _buffer);
      zeroSharesAll[dealer][i*_N+j] = _buffer[0];
    }
  }
  vector<vector<byte>> recBufs(_N), sendBufs(_N);
  for (int i=0; i<_N; i++) {
    vector<byte> tmp (nAuthBatchs * _N * _fieldByteSize);
    _field->elementVectorToByteVector(tagSharesAll[i], tmp);
    sendBufs[i].insert(sendBufs[i].end(), tmp.begin(), tmp.end());
    _field->elementVectorToByteVector(zeroSharesAll[i], tmp);
    sendBufs[i].insert(sendBufs[i].end(), tmp.begin(), tmp.end());
    recBufs[i].resize(nAuthBatchs*_N*2*_fieldByteSize, 0);
  }
  _comm.allToAll(sendBufs, recBufs, false); // no relay

  // -- compute dot product share + v share + o share per batch
  vector<vector<FieldType>>
    tauSharesAll(_N, vector<FieldType>(_N*nAuthBatchs, _zero));
  for (int i=0; i<_N; i++) {    // from i
    for (int j=0; j<_N; j++) {  // to j
      for (int k=0; k<nAuthBatchs; k++) { // kth share to auth
        FieldType dotProd = _zero;
        for (int l = 0; l<_authKeySize; l++) {
          dotProd += authBatchs[k][l] * _keyShares[l][i][j];
        }
        FieldType vShare = _field->
          bytesToElement(recBufs[i].data()+(k*_N+j)*_fieldByteSize);
        FieldType oShare = _field->
          bytesToElement(recBufs[i].data()+((nAuthBatchs+k)*_N+j)*_fieldByteSize);
        tauSharesAll[j][k*_N+i] = dotProd + vShare + oShare;
      }
    }
  }

  // -- distribut and store tag tag shares
  for (int i=0; i<_N; i++) {
    sendBufs[i].resize(nAuthBatchs*_N*_fieldByteSize);
    _field->elementVectorToByteVector(tauSharesAll[i], sendBufs[i]);
    recBufs[i].clear();
    recBufs[i].resize(nAuthBatchs*_N*_fieldByteSize, 0);
  }
  _comm.allToAll(sendBufs, recBufs, true); // use relay

  vector<FieldType> tauShareAll(_N);
  vector<int> inIds = allNButOne(_myId);
  for (int i=0; i<_N; i++) {    // from i
    for (int k=0; k<nAuthBatchs; k++) { // kth share to auth
      for (int j=0; j<_N; j++) {  // j'th share
        tauShareAll[j] = _field->
          bytesToElement(recBufs[j].data()+(k*_N+i)*_fieldByteSize);
      }
      vector<FieldType> inShares(_N, _zero);
      for (int j=0; j<_N-1; j++) {
        int pid = inIds[j];
        inShares[j] = tauShareAll[pid];
      }
      _mat_use_n_twist_vec[_myId].MatrixMult(inShares, _buffer);
      // _recTags[i][_authOffset+k] = _buffer[0];
      _recTags[i].push_back(_buffer[0]);
    }
  }
}

// TODO: impl
template <class FieldType>
bool ProtocolParty<FieldType>::checkTag() {
  return true;
}

// TODO: impl
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

template <class FieldType> bool ProtocolParty<FieldType>::prepareSeg() {
  // ---- # of random single shares:
  // 1. Padding before reveal
  //    -- randomization -> 1 time of [2]
  //    -- verifySharing (input and computation) -> 2 times
  //    -- checkRefresh -> 1 time
  //    in total: 5
  // 2. Random Coins
  //    [1]
  //    -- verifySharing after input and computation phases -> 2 times
  //    -- verifyOutput -> 1 time
  //    -- checkRefresh -> 1 time
  //    -- delinearization -> 1 time
  //    -- compress -> nCompressions + 2 times
  //    -- checkKeys -> 1 time (or 0)
  //    in total: nCompress + 8
  // 3. refresh
  //    [1] per Mult
  //    in total: segSize + 1 (in checkRefresh for random pad)
  _coinOffset = 0;
  _padOffset = 0;
  _refreshOffset = 0;
  int nCompressions = (int)(log(_segSize) / log(_K) + 0.5);
  int nCoinShares = nCompressions + 8;
  int nPadShares = 5;
  int nRefreshShares = _disp.hasDisp() ? _segSize + 1 : 0;
  makeRandShares(nCoinShares, _coinShares); // improve: optimize
  makeRandShares(nPadShares, _padShares);
  makeRandShares(nRefreshShares, _refreshShares);
  // ---- # of random double shares:
  // 1. Compute Mult gates in the circuit
  //    -- in total: numOfMultGates
  // 2. Compress Verification
  //    -- 2 * _K * (nCompressions + 1)
  int numDoubleShares = _segSize + (nCompressions+2)*_K*2;
  _doubleOffset = 0;
  makeRandDoubleShares(numDoubleShares, _doubleShares, _doubleKingIds);

  // ---- allocate transcript buffers for current segment
  _e1ShareIdx.clear();
  _recOShares.clear();
  _deltOShares.clear();
  return true;
}

template <class FieldType>
int ProtocolParty<FieldType>::processNonMult(int layerStart, int layerEnd) {
  int count = 0;
  for (int k = layerStart; k < layerEnd; k++) {
    if (_circuit.getGates()[k].gateType == ADD) {
      _gateShares[_circuit.getGates()[k].output] =
          _gateShares[_circuit.getGates()[k].input1] +
          _gateShares[_circuit.getGates()[k].input2];
    }
    else if (_circuit.getGates()[k].gateType == SUB)
    {
      _gateShares[_circuit.getGates()[k].output] =
          _gateShares[_circuit.getGates()[k].input1] -
          _gateShares[_circuit.getGates()[k].input2];
    } else if (_circuit.getGates()[k].gateType == SCALAR) {
      long scalar(_circuit.getGates()[k].input2);
      FieldType e = _field->GetElement(scalar);
      _gateShares[_circuit.getGates()[k].output] =
          _gateShares[_circuit.getGates()[k].input1] * e;
    } else if (_circuit.getGates()[k].gateType == SCALAR_ADD) {
      long scalar(_circuit.getGates()[k].input2);
      FieldType e = _field->GetElement(scalar);
      _gateShares[_circuit.getGates()[k].output] =
          _gateShares[_circuit.getGates()[k].input1] + e;
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
  auto start = _refreshShares.begin() + _refreshOffset;
  vector<FieldType> eShares(start, start + nShares);
  _refreshOffset += nShares;
  for (int i=0; i<nShares; i++) {
    eShares[i] += xShares[i];
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
      _kingShares.push_back(eSharesAll); // record collected shares as king
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
    recordDeltShares(oSharesAll, _deltOShares);
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
    _recOShares.push_back(oShare); // record received oShare
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
      xShares[multCounter] = _gateShares[gate.input1];
      yShares[multCounter] = _gateShares[gate.input2];
      multCounter++;
    }
  }                                                  

  if (_disp.hasCorrupt()) {
    refresh(xShares);
  }

  for (int i=0; i<nMults; i++) {
    e2Shares[i] = xShares[i] * yShares[i] -
      _doubleShares[ 2*(_doubleOffset + i) + 1 ];
  }
}

template <class FieldType>
int ProtocolParty<FieldType>::
processMult(int layerStart, int layerEnd, int nMults) {
  // collect e2 shares, refreshed if neccessary
  vector<FieldType> e2Shares;
  prepareE2Shares(layerStart, layerEnd, nMults, e2Shares);

  // -- batch mult
  vector<FieldType> e1Shares;
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
batchMultSingle(vector<FieldType>& e2Shares, vector<FieldType>& e1Shares) {
  int nMults = e2Shares.size();
  vector<vector<FieldType>> e2SharesLoad(_N);
  vector<int> loadSizes(_N, 0);
  for (int i=0; i<nMults; i++) {
    int king = _doubleKingIds[i + _doubleOffset];
    e2SharesLoad[king].push_back(e2Shares[i]);
    loadSizes[king]++;
  }
  int maxLoad = loadSizes[0];
  for (int i=0; i<_N; i++) {
    maxLoad = max(maxLoad, loadSizes[i]);
  }
  
  vector<vector<byte>> recBufs(_N, vector<byte>(maxLoad * _fieldByteSize)),
    e2SharesByte(_N, vector<byte>(maxLoad * _fieldByteSize));
  for (int i=0; i<_N; i++) {
    if (!_disp.isCorrupt(i)) {
      _field->elementVectorToByteVector(e2SharesLoad[i], e2SharesByte[i]);
    }
  }
  _comm.allToAll(e2SharesByte, recBufs, true);
  // -- all non-corr parties may be king
  int myLoad = e2SharesLoad[_myId].size();
  vector<FieldType> secrets(myLoad), e2ShareAll(_N);
  for (int i=0; i<myLoad; i++) {
    for (int j = 0; j < _N; j++) {
      e2ShareAll[j] =
        _field->bytesToElement(recBufs[j].data() + (i * _fieldByteSize));
    }
    _kingShares.push_back(e2ShareAll); // record collected shares as king
    secrets[i] = getSecret(e2ShareAll);
  }
  vector<vector<FieldType>> eSharesAll(_N, vector<FieldType>(myLoad, _zero));
  makeTShares(secrets, eSharesAll, true); // nonT shares all set to zero
  recordDeltShares(eSharesAll, _deltShares);
  vector<vector<byte>> e1SharesByte(_N, vector<byte>(maxLoad * _fieldByteSize, 0)),
    sendBufs(_N, vector<byte>(maxLoad * _fieldByteSize, 0));
  for (int i=0; i < _N; i++) {
    if (_TMask[i]) {
      _field->elementVectorToByteVector(eSharesAll[i], sendBufs[i]);
    }
  }
  _comm.allToT(sendBufs, e1SharesByte);
  // _comm.allToAll(sendBufs, e1SharesByte, false);
  
  // for (int i=0; i<_N; i++) {
  //   e1SharesByte[i].resize(e2SharesLoad[i].size() * _fieldByteSize);
  //   if (!_disp.isCorrupt(i)) {
  //     _comm.kingToT(i, e1SharesByte[i], sendBufs);
  //   }
  // }
  // -- convert byte to elements
  vector<int> idx(_N, 0);
  e1Shares.resize(nMults);
  for (int i = 0; i < nMults; i++) {
    int king = _doubleKingIds[i + _doubleOffset];
    e1Shares[i] = _field->bytesToElement(e1SharesByte[king].data() +
                                         (idx[king]++) * _fieldByteSize);
    _e1ShareIdx.push_back(_recShares.size()); // record e1ShareIdx
    // record received shares and their dealers
    _dealerIds.push_back(king);
    _recShares.push_back(e1Shares[i]);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
multVec(vector<FieldType>& a, vector<FieldType>& b,
          vector<vector<FieldType>>& transcripts, int groupSize) {
  // assert( a.size() == b.size() );
  int totalLength = a.size();
  int nMults = (totalLength + groupSize - 1) / groupSize;
  transcripts.clear();
  transcripts.resize(5, vector<FieldType>(nMults));
  
  // -- generate the 2t-sharings for xy - r
  // -- also output the transcripts of r2, e2
  for (int group_i = 0; group_i < nMults; group_i++) {
    int group_end = ((group_i + 1) * groupSize > totalLength) ?
      totalLength : (group_i + 1) * groupSize;
    for (int i = group_i * groupSize; i < group_end; i++) {
      transcripts[3][group_i] += a[i] * b[i];
    }
    transcripts[1][group_i] =  _doubleShares[(_doubleOffset + group_i)*2 +1];    
    transcripts[3][group_i] = transcripts[3][group_i] - transcripts[1][group_i];
  }

  // batchMultSpread(transcripts[3], transcripts[2]);
  batchMultSingle(transcripts[3], transcripts[2]);
  // -- compute xy - r + [r]_t = t-sharing of xy
  // -- also output the transcripts of r1, e1
  for (int group_i = 0; group_i < nMults; group_i++) {
    transcripts[0][group_i] = _doubleShares[(_doubleOffset + group_i)*2];
    transcripts[4][group_i] = transcripts[0][group_i] + transcripts[2][group_i];
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
  vector<FieldType> lambdaShare(1, _coinShares[_coinOffset++]);
  vector<FieldType> lambda(1);
  openTShares(1, true, lambdaShare, lambda);
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
    sendElems[outParty].push_back(_gateShares[outWire]);
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
      _kingShares.push_back(outShareAll); // record collected shares as king
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
