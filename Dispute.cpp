#include "Dispute.h"


// TODO: optimize
void Dispute::reset(int N) {
  _N = N;
  _T = (N - 1)/ 2;
  _nCorr = 0;
  _newDisp = true;
  _hasDisp = false;
  _disp.clear();
  _disp.resize(N, vector<bool>(N, false));
  _corr.clear();
  _corr.resize(N, false);
  _dispSize.clear();
  _dispSize.resize(N, 0);
}


void Dispute::addDispPairs(int p1, int p2) {
  assert(p1 < _N);
  assert(p2 < _N);
  assert(!_disp[p1][p2]);
  assert(!_disp[p2][p1]);

  _newDisp = true;
  _hasDisp = true;
  _disp[p1][p2] = true;
  _disp[p2][p1] = true;

  _dispSize[p1]++;
  _dispSize[p2]++;
    
  if (_dispSize[p1] > _T) {
    addCorrParty(p1);
  }

  if (_dispSize[p2] > _T) {
    addCorrParty(p2);
  }
}

void Dispute::addCorrParty(int p) {
  _corr[p] = true;
  _nCorr++;

  _newDisp = true;
  _hasDisp = true;
  _dispSize[p] = _N;
  for(int i=0; i<_N; i++) {
    _disp[p][i] = true;
    if (!_disp[i][p]) {
      _disp[i][p] = true;
      _dispSize[i]++;
    }
  }
}

// TODO: remember for repeated use
int Dispute::tSet(vector<int> &TSet) {
  int king = 0;
  while (_corr[king]) {
    king++;
  }

  TSet.resize(_T + 1);
  TSet[0] = king;
  int count = 1;
  for (int i=0; i<_N; i++) {
    if (count >= _T+1) {
      break;
    }
    if (!_disp[king][i] && king != i) {
      TSet[count++] = i;
    }
  }

  assert(count == _T+1);
  return king;
}

int Dispute::nonTSet(vector<int> &NonTSet) {
  NonTSet.resize(_N - _T - 1);
  vector<bool> TMask;
  int king = tMask(TMask);
  int nonTCount = 0;
  for (int i=0; i<_N; i++) {
    if (!TMask[i]) {
      NonTSet[nonTCount++] = i;
    }
  }
  return king;
}

int Dispute::tAndNonTSet(vector<int> &TSet, vector<int>& NonTSet) {
  int king = tSet(TSet);
  nonTSet(NonTSet);
  return king;  

}

void Dispute::tAndNonTSetP(int p, vector<int>& TSet, vector<int>& NonTSet) {
  assert(!_corr[p]);
  vector<int> NonDispSet;
  int nDisp = dispAndNonDispSet(p, NonTSet, NonDispSet);
  TSet.resize(_T+1);
  TSet[0] = p;
  NonTSet.resize(_N - _T - 1);
  int tCount = 1;
  int nonTCount = nDisp;
  for (int i=0; i<_N-nDisp; i++) {
    int curId = NonDispSet[i];
    if (curId == p) {
      continue;
    }
    if (tCount < _T+1) {
      TSet[tCount++] = curId;
    } else {
      NonTSet[nonTCount++] = curId;
    }
  }
}


int Dispute::tMask(vector<bool> &TMask) {
  TMask.clear();
  TMask.resize(_N, false);
  vector<int> TSet;
  int king = tSet(TSet);
  for (int t : TSet) {
    TMask[t] = true;
  }
  return king;
}

void Dispute::tMaskP(int p, vector<bool> &TMask) {
  TMask.clear();
  TMask.resize(_N, false);
  if (_corr[p]) {
    return;
  }
  vector<int> TSetP, tmp;
  tAndNonTSetP(p, TSetP, tmp);
  for (int i : TSetP) {
    TMask[i] = true;
  }
}

// TODO: use memory and optimize
void Dispute::tMaskPRev(int p, vector<bool> &TMaskRev) {
  TMaskRev.clear();
  TMaskRev.resize(_N, false);
  if (_corr[p]) {
    return;
  }
  TMaskRev[p] = true;
  vector<bool> TMaskI;
  for (int i=0; i<_N; i++) {
    if ( i == p || _disp[p][i]) {
      continue;                 // not in TSet of disp parties
    }
    tMaskP(i, TMaskI);
    TMaskRev[i] = TMaskI[p];
  }
}

// TODO: remember for repeated use
bool Dispute::
rlyeeVecs(int p, vector<int>& relay, vector<bool>& relayerMask,
                 vector<vector<int>>& relayerLoad) {
  assert(!_corr[p]);
  relay.clear();
  relay.resize(_N, -1);         // -1 means no relay needed
  relayerMask.clear();
  relayerMask.resize(_N, false); // false means not a relayer
  relayerLoad.clear();
  relayerLoad.resize(_N);       // empty vector means not a relayer
  if (_dispSize[p] == 0) {
    return false;               // not a relayee
  } // else: is a relayee
  
  for (int i=0; i<_N; i++) {
    if ((!_corr[i]) && _disp[p][i]) { // requires a relay
      int r = 0;
      while (_disp[p][r] || _disp[i][r]) { r++; }
      relay[i] = r;
      relayerMask[r] = true;
      relayerLoad[r].push_back(i);
    }
  }
  return true;
}

// TODO: memorize and reuse
bool Dispute::
rlyerVecs(int p, vector<bool>& relayeeMask, vector<vector<int>>& relayeeLoad) {
  assert(!_corr[p]);
  relayeeMask.clear();
  relayeeMask.resize(_N, false);
  relayeeLoad.clear();
  relayeeLoad.resize(_N);
  bool isRelayer = false;
  vector<int> relay;
  vector<bool> relayerMask;
  vector<vector<int>> relayerLoad;
  for (int i=0; i<_N; i++) {
    if ((!_corr[i]) &&  _dispSize[i] > 0) {     // is an relayee
      rlyeeVecs(i, relay, relayerMask, relayerLoad);
      if (relayerMask[p]) {
        isRelayer = true;
        relayeeMask[i] = true;
        relayeeLoad[i] = relayerLoad[p];
      }
    }
  }
  return isRelayer;
}

void Dispute::corrSet(vector<int>& corrs) {
  corrs.resize(_nCorr);
  int count = 0;
  for (int i=0; i<_N; i++) {
    if (!_corr[i]) {
      corrs[count++] = i;
    }
  }
}

int Dispute::dispAndNonDispSet(int p, vector<int>& DispSet,
                               vector<int>& NonDispSet) {
  DispSet.clear();
  DispSet.resize(_dispSize[p]);
  NonDispSet.clear();
  NonDispSet.resize(_N - _dispSize[p]);
  int dispCount = 0;
  int nonDispCount = 0;
  for (int i=0; i < _N; i++) {
    if (_disp[p][i]) {          // i is disp
      DispSet[dispCount++] = i;
    } else {                    // i is nonDisp
      NonDispSet[nonDispCount++] = i;
    }
  }
  return dispCount;
}

int Dispute::twistSet(int p, int target, vector<int>& TwistSet) {
  assert(!_disp[p][target]);
  vector<int> DispSet, NonDispSet;
  int nDisp = dispAndNonDispSet(p, DispSet, NonDispSet);
  TwistSet = DispSet;
  TwistSet.resize(_N-1);
  int idx = nDisp;
  for (int npid : NonDispSet) {
    if (npid != target) {
      TwistSet[idx++] = npid;
    }
  }
  return nDisp;
}

Dispute::Dispute() {
}

Dispute::~Dispute(){
}
