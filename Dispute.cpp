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

void Dispute::tAndNonTSetP(int p, vector<int>& TSet, vector<int>& NonTSet) {
  vector<int> rotate(_N);
  for (int i=0; i<_N; i++) {
    rotate[i] = (i + p) % _N;
  }
  TSet.clear();
  TSet.insert(TSet.end(), rotate.begin(), rotate.begin()+_T+1);
  NonTSet.clear();
  NonTSet.insert(NonTSet.end(), rotate.begin()+_T+1, rotate.end());
}


void Dispute::tMaskP(int p, vector<bool> &TMask) {
  TMask.clear();
  TMask.resize(_N, false);
  for (int i=0; i<_T+1; i++) {
    TMask[ (i + p)%_N ] = true;
  }
}

// TODO: use memory and optimize
void Dispute::tMaskPRev(int p, vector<bool> &TMaskRev) {
  TMaskRev.clear();
  TMaskRev.resize(_N, false);
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

Dispute::Dispute() {
}

Dispute::~Dispute(){
}
