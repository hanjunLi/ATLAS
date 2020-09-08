// author: Hanjun Li <lihanjun1212@gmail.com>
#ifndef TAPE_H_
#define TAPE_H_

#include <stdlib.h>
#include <algorithm>

#include <libscapi/include/primitives/Matrix.hpp>
#include <vector>
#include <iostream>
#include <libscapi/include/primitives/Mersenne.hpp>
#include <libscapi/include/infra/Common.hpp>

#include "Share.h"

using namespace std;

// short helpers
__inline__
void recordIdxs(unsigned int idx, int num,
                vector<unsigned int>& tape){
  for (int i=0; i<num; i++) {
    tape.push_back(idx + i);
  }
}

// internal data class
template <class FieldType>
class SegTape{

public:
  // relayer for king stores its relay share
  // TODO: maintain indices for _relayE2Share
  vector<vector<FieldType>> _relayE2Shares; // [j]th e2 share from [i] to king

  vector<vector<FieldType>> _recFullShares; // received shares as king
  vector<unsigned int> _recE2Idxs; // received e2 full shares as king
                                // points to recFullShares
  vector<unsigned int> _recO2Idxs; // received O2 full shares as king
                                // points to recFullShares
  vector<unsigned int> _recOutIdxs; // received Output full shares as king
                                // points to recFullShares
  vector<unsigned int> _dealtE1Idxs; // dealt e1 shares as king
                                // points to dealtShares
  vector<unsigned int> _recE1Idxs; // received e1 shares
                                // points to recShares
  vector<unsigned int> _dealtOIdxs; // dealt O shares as king
                                // points to dealtShares
  vector<unsigned int> _recOIdxs; // received O shares
                                // points to recShares
  // len = rl2T
  vector<FieldType> _rec2TShares;      // received 2T shares
  vector<unsigned short> _dealer2TIds; // dealer2TIds[i] dealt 2T share [i]
  vector<unsigned int> _dealt2TIdxs;   // -1 if dealer is not me
  unsigned int _rL2T = 0;
  // len = dl2T
  vector<vector<FieldType>> _dealt2TShares; // dealt 2T shares
  unsigned int _dL2T = 0;                   // current length of dl

  void reset(int N) {
    // _relayE1Share.clear();
    _relayE2Shares.clear(); _relayE2Shares.resize(N);
    _recFullShares.clear();
    _recE2Idxs.clear(); _recO2Idxs.clear(); _recOutIdxs.clear();    
    _dealtE1Idxs.clear(); _recE1Idxs.clear();
    _dealtOIdxs.clear(); _recOIdxs.clear();

    _rec2TShares.clear(); _dealer2TIds.clear(); _dealt2TIdxs.clear();
    _rL2T = 0;
    _dealt2TShares.clear();
    _dL2T = 0;
  }
};

template <class FieldType>
class PartyTape{

  // NOTE:
  // unsigned int has max range = 4,294,967,295 ~ 4 billion
  // should not use larger than 4 billion circuits
  // unless we change to long long
  
private:
  // received shares tape (len = rl)
  vector<FieldType> _recShares; // received shares in the protocol
  vector<unsigned short> _dealerIds;  // dealerIds[i] dealt share [i]
  vector<unsigned int> _dealtIdxs;    // -1 if dealer is not me
  vector<unsigned int> _batchIdxs;    // authentication batch
  vector<unsigned int> _segRecBegins; // recShares[segRecBegins[i]] is
                                // the first share in seg i
  vector<unsigned int> _segBatchBegins; // segBatchEnds[i] is
                                // the first in seg i

  // verification tags tape (len = N x rl)
  vector<vector<FieldType>> _recTauTags; // tau from i of jth batch
  vector<vector<FieldType>> _sentNuTags; // nu for i about jth batch

  unsigned int _rL = 0;           // current lenth of rl
  unsigned int _verifyOffset = 0; // next recShare to verify
  unsigned int _authOffset = 0;   // next recShare to authenticate
  unsigned int _nextAuthBatch = 0; // next authentication batch idx
  unsigned int _authKeySize = 0;

  // dealt shares tape (len = dl)
  vector<vector<FieldType>> _dealtShares; // dealt shares in the protocol
  vector<unsigned int> _segDealtBegins;   // similar to recSegBegins
  unsigned int _dL = 0;                   // current length of dl

  unsigned short _N;                          // N parties
  unsigned short _myId;                       // my party ID
  TemplateField<FieldType>* _field = nullptr; // passed from outer protocol
  SegTape<FieldType> _segTape;                // record per-seg tape

public:
  // PartyTape() {}
  // ~PartyTape() {}

  void updateAuthOffset() { _authOffset = _rL; }
  void updateVerifyOffset() { _verifyOffset = _rL; }
  unsigned int getRecTapeLength() { return _rL; }
  unsigned int getDealtTapeLength() { return _dL; }
  unsigned int getRelayTapeLength() { return _segTape._relayE2Shares.size(); }
  unsigned int getNumAuthBatches() { return _nextAuthBatch; }
  unsigned int getPrevNumAuthBatches() { return _segBatchBegins.back(); }
  FieldType getRecShareVal(unsigned int idx) { return _recShares[idx]; }
  Share<FieldType> getRecShare(unsigned int idx) {
    Share<FieldType> result;
    result._value = _recShares[idx];
    result._idxs = vector<unsigned int>(1, idx);
    result._coeffs = vector<FieldType>(1, *(_field->GetOne()));
    return result;
  }
  
  // void getDealtShare(unsigned int idx, vector<FieldType>& dealtFullShare) {
  //   dealtFullShare = _dealtShares[idx];
  // }

  void reset(int N, int myId, int keySize, TemplateField<FieldType>* field) {
    // for received shares
    _recShares.clear(); _dealerIds.clear();
    _dealtIdxs.clear(); _batchIdxs.clear();
    _segRecBegins.clear(); _segBatchBegins.clear();
    _recTauTags.clear(); _recTauTags.resize(N);
    _sentNuTags.clear(); _sentNuTags.resize(N);
    
    _rL = 0; _verifyOffset = 0; _authOffset = 0;
    _nextAuthBatch = 0; _authKeySize = keySize;

    // for dealt shares
    _dealtShares.clear(); _segDealtBegins.clear();
    _dL = 0;

    // common
    _N = N; _myId = myId; _field = field;
    _segTape.reset(_N);
  }

  void startNewSeg() {
    _segRecBegins.push_back(_rL);
    _segBatchBegins.push_back(_nextAuthBatch);
    _segDealtBegins.push_back(_dL);
    _segTape.reset(_N);
  }

  void clearLastSeg() {
    if (_segRecBegins.size() > 0) {
      unsigned int prevSegEnd = _segRecBegins.back();
      _segRecBegins.pop_back();
      _recShares.resize(prevSegEnd);
      _dealerIds.resize(prevSegEnd);
      _batchIdxs.resize(prevSegEnd);
      _verifyOffset = min(prevSegEnd, _verifyOffset);
      _authOffset = min(prevSegEnd, _authOffset);
      _batchIdxs.resize(_authOffset);
      _rL = prevSegEnd;
    }

    if (_segBatchBegins.size() > 0) {
      _nextAuthBatch = _segBatchBegins.back();
      _segBatchBegins.pop_back();
      for (int i=0; i<_N; i++) {
        _recTauTags[i].resize(_nextAuthBatch);
        _sentNuTags[i].resize(_nextAuthBatch);
      }
    }

    if (_segDealtBegins.size() > 0) {
      unsigned int prevSegEnd = _segDealtBegins.back();
      _segDealtBegins.pop_back();
      _dealtShares.resize(prevSegEnd);
      _dL = prevSegEnd;
    }

    _segTape.reset(_N);
  }

  // NOTE: currently not used
  void recordRecShare(FieldType share, unsigned short dealer) {
    _recShares.push_back(share);
    _dealerIds.push_back(dealer);
    _rL++;
    if (dealer == _myId) {
      _dealtIdxs.push_back(_dL);
      _dL++;
    } else {
      _dealtIdxs.push_back(-1);
    }
  }

  void recordRecShares(vector<FieldType>& shares, unsigned short dealer) {
    int nNewShares = shares.size();
    _recShares.insert(_recShares.end(), shares.begin(), shares.end());
    _dealerIds.insert(_dealerIds.end(), nNewShares, dealer);
    _rL += nNewShares;
    if (dealer == _myId) {
      for (int i=0; i<nNewShares; i++) { _dealtIdxs.push_back(_dL + i); }
      _dL += nNewShares;
    } else {
      _dealtIdxs.insert(_dealtIdxs.end(), nNewShares, -1);
    }
  }

  // NOTE: currently not used
  void recordDealtShare(vector<FieldType>& fullShare) {
    assert(fullShare.size() == _N);
    _dealtShares.push_back(fullShare);
    // _dL is updated by recordRecShare(s)
  }

  void recordDealtShares(vector<vector<FieldType>>& fullShares) {
    // NOTE: transpse fullShares before record
    assert(fullShares.size() == _N);
    for (int i=1; i<_N; i++) {
      assert(fullShares[i].size() == fullShares[0].size());
    }
    int nShares = fullShares[0].size();
    vector<FieldType> tmp(_N);
    for (int idx=0; idx<nShares; idx++) {
      for (int id=0; id<_N; id++) {
        tmp[id] = fullShares[id][idx];
      }
      _dealtShares.push_back(tmp);
    }
    // _dL is updated by recordRecShare(s)
  }

  void recordSentTag(unsigned short toId, FieldType nu) {
    _sentNuTags[toId].push_back(nu);
  }

  void recordRecTag(unsigned short fromId, FieldType tau) {
    _recTauTags[fromId].push_back(tau);
  }
  
  // void linearCombChunks(const vector<unsigned int>& shareIdxs,
  //                       const vector<int>& coeffs, unsigned int chunkSize,
  //                       unsigned int start, unsigned int end,
  //                       vector<FieldType>& chunks) { // output
  //   // TODO: fix
  //   int nChunks = (end - start + chunkSize -1) / chunkSize;
  //   FieldType _zero = *(_field->GetZero());
  //   chunks.assign(nChunks, _zero);

  //   int nShares = shareIdxs.size();
  //   for (int i=0; i<nShares; i++) {
  //     int curIdx = shareIdxs[i];
  //     int curChunkIdx = (curIdx - start) / chunkSize;
  //     chunks[i] += _field->GetElement(coeffs[i]) * _recShares[curIdx];
  //   }
  // }

  void decompose(const Share<FieldType>& share, // in 
                 vector<FieldType>& myShareDecomp, // out
                 vector<FieldType>& myFullShare) { // out
    myShareDecomp.assign(_N, *(_field->GetZero()));
    myFullShare.assign(_N, *(_field->GetZero()));
    int length = share._idxs.size();
    for (int i=0; i<length; i++) {
      unsigned int recIdx = share._idxs[i];
      FieldType coeff = share._coeffs[i];
      unsigned short dealerId = _dealerIds[recIdx];
      myShareDecomp[dealerId] += _recShares[recIdx] * coeff;
      if (dealerId == _myId) {
        unsigned int dealtIdx = _dealtIdxs[recIdx];
        for (int j=0; j<_N; j++) {
          myFullShare[j] += _dealtShares[dealtIdx][j] * coeff;
        }
      }
    }
  }

  void getDealerSubShare(unsigned short dealer, const Share<FieldType>& share,
                         Share<FieldType>& dealerSubShare) {
    int length = share._idxs.size();
    dealerSubShare._value = *(_field->GetZero());
    dealerSubShare._idxs.clear();
    dealerSubShare._coeffs.clear();
    for (int i=0; i<length; i++) {
      unsigned int recIdx = share._idxs[i];
      FieldType coeff = share._coeffs[i];
      if (_dealerIds[recIdx] == dealer) {
        dealerSubShare._value += _recShares[recIdx] * coeff;
        dealerSubShare._idxs.push_back(recIdx);
        dealerSubShare._coeffs.push_back(coeff);
      }
    }
  }

  void getChunkSubShares(const Share<FieldType>& share,
                         vector<Share<FieldType>>& chunkSubShares) {
    unsigned int length = share._idxs.size();
    assert(length > 0); // not an empty share. should not happen
    chunkSubShares.clear();
    chunkSubShares.resize(_N);

    unsigned int minBatchIdx = _batchIdxs[share._idxs[0]];
    unsigned int maxBatchIdx = minBatchIdx;
    for (unsigned int recIdx : share._idxs) {
      unsigned int curBatchIdx = _batchIdxs[recIdx];
      minBatchIdx = (curBatchIdx < minBatchIdx) ? curBatchIdx : minBatchIdx;
      maxBatchIdx = (curBatchIdx > maxBatchIdx) ? curBatchIdx : maxBatchIdx;
    }

    unsigned int inc = (maxBatchIdx - minBatchIdx) / _N +1;
    for (int i=0; i<length; i++) {
      unsigned int recIdx = share._idxs[i];
      FieldType coeff = share._coeffs[i];
      unsigned int curBatchIdx = _batchIdxs[recIdx];
      unsigned short chunkIdx = (curBatchIdx - minBatchIdx) / inc;
      chunkSubShares[chunkIdx]._value += _recShares[recIdx] * coeff;
      chunkSubShares[chunkIdx]._idxs.push_back(recIdx);
      chunkSubShares[chunkIdx]._coeffs.push_back(coeff);
    }
  }

  unsigned int getBatchSubShares(const Share<FieldType>& share,
                                 vector<Share<FieldType>>& batchSubShares) {
    unsigned int length = share._idxs.size();
    assert(length > 0); // not an empty share. should not happen
    
    unsigned int minBatchIdx = _batchIdxs[share._idxs[0]];
    unsigned int maxBatchIdx = minBatchIdx;
    for (unsigned int recIdx : share._idxs) {
      unsigned int curBatchIdx = _batchIdxs[recIdx];
      minBatchIdx = (curBatchIdx < minBatchIdx) ? curBatchIdx : minBatchIdx;
      maxBatchIdx = (curBatchIdx > maxBatchIdx) ? curBatchIdx : maxBatchIdx;
    }
    unsigned int nBatches = maxBatchIdx - minBatchIdx + 1;
    batchSubShares.clear();
    batchSubShares.resize(nBatches);

    for (int i=0; i<length; i++) {
      unsigned int recIdx = share._idxs[i];
      FieldType coeff = share._coeffs[i];
      unsigned int curBatchIdx = _batchIdxs[recIdx];
      unsigned short chunkIdx = (curBatchIdx - minBatchIdx);
      batchSubShares[chunkIdx]._value += _recShares[recIdx] * coeff;
      batchSubShares[chunkIdx]._idxs.push_back(recIdx);
      batchSubShares[chunkIdx]._coeffs.push_back(coeff);
    }
    return nBatches;
  }

  unsigned int getBatchInfo(const Share<FieldType>& share,
                            vector<unsigned int>& batchRecIdxs,
                            vector<FieldType>& batchCoeffs) {
    unsigned int length = share._idxs.size();
    assert(length > 0); // not an empty share
    unsigned int batchIdx = _batchIdxs[share._idxs[0]];
    for (unsigned int recIdx : share._idxs) {
      assert(_batchIdxs[recIdx] == batchIdx ); // all same batch
    }
    assert(batchIdx < _segBatchBegins.back()); // is a batch from prev seg

    batchRecIdxs.clear();
    batchCoeffs.clear();
    auto gtPtr = upper_bound(_segBatchBegins.begin(),
                             _segBatchBegins.end(), batchIdx);
    unsigned int segEnd = *(gtPtr);
    unsigned int segBegin = *(gtPtr-1);
    unsigned int curCoeffIdx = 0;
    for (unsigned int recIdx=segBegin; recIdx<segEnd; recIdx++) {
      if (_batchIdxs[recIdx] == batchIdx) {
        batchRecIdxs.push_back(recIdx);
        if (curCoeffIdx < length && recIdx == share._idxs[curCoeffIdx]) {
          batchCoeffs.push_back(share._coeffs[curCoeffIdx]);
          curCoeffIdx++;
        } else {
          batchCoeffs.push_back(*(_field->GetZero()));
        }
      }
    }
    return batchIdx;
  }

  void getDealerBatchFullShares(unsigned int batchIdx,
                                vector<vector<FieldType>>& fullShares) {
    fullShares.assign(_N, vector<FieldType>(_authKeySize));
    for (int keyIdx=0; keyIdx<_authKeySize; keyIdx++) {
      unsigned int absIdx = batchIdx * _authKeySize + keyIdx;
      for (int recId=0; recId<_N; recId++) {
        fullShares[recId][keyIdx] = _dealtShares[absIdx][recId];
      }
    }
  }

  void getTagVals(unsigned int batchIdx,
                  vector<FieldType>& tagVals) {
    tagVals.assign(_N, *(_field->GetZero()));
    for (unsigned short fromId =0; fromId < _N; fromId++) {
      tagVals[fromId] = _recTauTags[fromId][batchIdx];
    }
  }

  FieldType getPadVal(unsigned short partyId, unsigned int batchIdx) {
    return _sentNuTags[partyId][batchIdx];
  }

  FieldType getDealtShareValFor(unsigned short accused,
                                const Share<FieldType> share) {
    for (unsigned int recIdx : share._idxs) {
      assert(_dealtIdxs[recIdx] > 0); // I'm the dealer
    }

    FieldType result = *(_field->GetZero());
    unsigned int length = share._idxs.size();
    for (int i=0; i<length; i++) {
      unsigned int recIdx = share._idxs[i];
      FieldType coeff = share._coeffs[i];
      result += _dealtShares[_dealtIdxs[recIdx]][accused] * coeff;
    }
    return result;
  }

  void decompose2T(const Share<FieldType>& share, // in
                   vector<FieldType>& myShareDecomp, // out
                   vector<FieldType>& myFullShare) { // out
    myShareDecomp.assign(_N, *(_field->GetZero()));
    myFullShare.assign(_N, *(_field->GetZero()));
    int length = share._idxs.size();
    for (int i=0; i<length; i++) {
      unsigned int rec2TIdx = share._idxs[i];
      FieldType coeff = share._coeffs[i];
      unsigned short dealer2TId = _segTape._dealer2TIds[rec2TIdx];
      myShareDecomp[dealer2TId] += _segTape._rec2TShares[rec2TIdx] * coeff;
      if (dealer2TId == _myId) {
        unsigned int dealt2TIdx = _segTape._dealt2TIdxs[rec2TIdx];
        for (int j=0; j<_N; j++) {
          myFullShare[j] += _segTape._dealt2TShares[dealt2TIdx][j] * coeff;
        }
      }
    }
  }

  int makeAuthBatch (vector<vector<FieldType>>& batchToAuth,
                    vector<int>& batchDealer) {
    auto authStart = _recShares.begin() + _authOffset;
    auto authDealerStart = _dealerIds.begin() + _authOffset;
    int nToAuth = _rL - _authOffset;
    // -- compute [key] * [dealer batch] shares
    vector<unsigned int> dealerBatchIdx(_N);
    for (int i=0; i<_N; i++) {
      dealerBatchIdx[i] = _nextAuthBatch++; // assign batch idx for each dealer
    }
    vector<vector<FieldType>> curBatch(_N);
    for (int i=0; i<nToAuth; i++) {
      int dealer = authDealerStart[i];
      curBatch[dealer].push_back(authStart[i]);
      _batchIdxs.push_back(dealerBatchIdx[dealer]); // record batch idx
      if (curBatch[dealer].size() == _authKeySize) { // store current batch
        batchToAuth.push_back(curBatch[dealer]);
        dealerBatchIdx[dealer] = _nextAuthBatch++; // update batch idx
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

  void makeToVerifyShares(vector<Share<FieldType>>& toVerify) {
    int nShares = _rL - _verifyOffset;
    toVerify.resize(nShares);
    for (int i=0; i<nShares; i++) {
      toVerify[i]._value = _recShares[_verifyOffset + i];
      toVerify[i]._idxs = vector<unsigned int>(1, _verifyOffset + i);
      toVerify[i]._coeffs = vector<FieldType>(1, *(_field->GetOne()));
    }
  }
  
  unsigned int getDealtE1Length() { return _segTape._dealtE1Idxs.size(); }
  unsigned int getRecE2Length() { return _segTape._recE2Idxs.size(); }
  // unsigned int getDealtOLength() { return _dealtOIdxs.size(); }
  // unsigned int getRecO2Length() { return _recO2Idxs.size(); }
  unsigned int getRec2TTapeLength() { return _segTape._rL2T; }

  void getRecO2FullShare(int nth, vector<FieldType>& O2FullShare) {
    O2FullShare = _segTape._recFullShares[_segTape._recO2Idxs[nth]];
  }
  void getRecE2FullShare(int nth, vector<FieldType>& E2FullShare) {
    E2FullShare = _segTape._recFullShares[_segTape._recE2Idxs[nth]];
  }
  void getDealtOShare(int nth, vector<FieldType>& OFullShare) {
    OFullShare = _dealtShares[_segTape._dealtOIdxs[nth]];
  }
  void getDealtE1Share(int nth, vector<FieldType>& E1FullShare) {
    E1FullShare = _dealtShares[_segTape._dealtE1Idxs[nth]];
  }
  FieldType getRecOShareVal(int nth) {
    return _recShares[_segTape._recOIdxs[nth]];
  }
  Share<FieldType> getRecOShare(int nth) {
    Share<FieldType> result;
    unsigned int recIdx = _segTape._recOIdxs[nth];
    result._value = _recShares[recIdx];
    result._idxs = vector<unsigned int>(1, recIdx);
    result._coeffs = vector<FieldType>(1, *(_field->GetOne()));
    return result;
  }
  FieldType getRecE1ShareVal(int nth) {
    return _recShares[_segTape._recE1Idxs[nth]];
  }
  Share<FieldType> getRecE1Share(int nth) {
    Share<FieldType> result;
    unsigned int recIdx = _segTape._recE1Idxs[nth];
    result._value = _recShares[recIdx];
    result._idxs = vector<unsigned int>(1, recIdx);
    result._coeffs = vector<FieldType>(1, *(_field->GetOne()));
    return result;
  }

  void recordE2Share(vector<FieldType>& full2TShare) {
    _segTape._recE2Idxs.push_back(_segTape._recFullShares.size());
    _segTape._recFullShares.push_back(full2TShare);
  }

  void recordE2Relay(vector<vector<FieldType>>& relayShares) {
    for (int i = 0; i < _N; i++) {
      if (relayShares[i].size() > 0) {
        _segTape._relayE2Shares[i].insert(_segTape._relayE2Shares[i].end(),
                                         relayShares[i].begin(),
                                         relayShares[i].end());
      }
    }
  }

  void getE2RelayShare(unsigned int idx, vector<FieldType>& relayShares) {
    relayShares.resize(_N);
    for (int i = 0; i < _N; i++) {
      if (_segTape._relayE2Shares[i].size() == 0) {
        relayShares[_N] = *(_field->GetZero());
      } else {
        assert(idx < _segTape._relayE2Shares[i].size() > idx);
        relayShares[_N] = _segTape._relayE2Shares[i][idx];
      }
    }
  }

  void recordO2Share(vector<FieldType>& full2TShare) {
    _segTape._recO2Idxs.push_back(_segTape._recFullShares.size());
    _segTape._recFullShares.push_back(full2TShare);
  }

  void recordOutShare(vector<FieldType>& fullTShare) {
    _segTape._recOutIdxs.push_back(_segTape._recFullShares.size());
    _segTape._recFullShares.push_back(fullTShare);
  }

  void recordDealtE1Shares(vector<vector<FieldType>>& fullShares) {
    int idxStart = _dL;
    recordDealtShares(fullShares);
    recordIdxs(idxStart, fullShares[0].size(), _segTape._dealtE1Idxs);
  }

  void recordDealtOShares(vector<vector<FieldType>>& fullShares) {
    int idxStart = _dL;
    recordDealtShares(fullShares);
    recordIdxs(idxStart, fullShares[0].size(), _segTape._dealtOIdxs);
  }

  void recordRecE1Shares(vector<FieldType>& shares, unsigned short dealer) {
    int idxStart = _rL;
    recordRecShares(shares, dealer);
    recordIdxs(idxStart, shares.size(), _segTape._recE1Idxs);
  }

  // void recordE1Relay(vector<vector<FieldType>>& relayShares) {
  //   for (int i = 0; i < _N; i++) {
  //     _segTape._relayE1Share.
  //       insert(_segTape._relayE1Share.end(),
  //              relayShares[i].begin(), relayShares[i].end());
  //   }
  // }

  // NOTE: currently not used
  void recordRecOShares(vector<FieldType>& shares, unsigned short dealer) {
    int idxStart = _rL;
    recordRecShares(shares, dealer);
    recordIdxs(idxStart, shares.size(), _segTape._recOIdxs);
  }

  void recordRecOShare(FieldType& share, unsigned short dealer) {
    int idx = _rL;
    recordRecShare(share, dealer);
    _segTape._recOIdxs.push_back(idx);
  }

  void recordDealt2TShares(vector<vector<FieldType>>& full2TShares) {
    // NOTE: transpse fullShares before record
    assert(full2TShares.size() == _N);
    for (int i=1; i<_N; i++) {
      assert(full2TShares[i].size() == full2TShares[0].size());
    }
    int nShares = full2TShares[0].size();
    vector<FieldType> tmp(_N);
    for (int idx=0; idx<nShares; idx++) {
      for (int id=0; id<_N; id++) {
        tmp[id] = full2TShares[id][idx];
      }
      _segTape._dealt2TShares.push_back(tmp);
    }
  }

  void recordRec2TShare(FieldType& share2T, unsigned short dealer) {
    _segTape._rec2TShares.push_back(share2T);
    _segTape._dealer2TIds.push_back(dealer);
    _segTape._rL2T++;
    if (dealer == _myId) {
      _segTape._dealt2TIdxs.push_back(_segTape._dL2T++);
    } else {
      _segTape._dealt2TIdxs.push_back(-1);
    }
  }
};


#endif
