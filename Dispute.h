// author: Hanjun Li <lihanjun1212@gmail.com>
#ifndef DISPUTE_H_
#define DISPUTE_H_

#include <assert.h>
#include <vector>
#include <unordered_set>

using namespace std;

class Dispute{

private:
  int _N, _T, _nCorr;
  bool _hasDisp = false;
  bool _newDisp = false;
  vector<int> _dispSize;
  vector<vector<bool>> _disp;   // bit vector for disputes
  vector<bool> _corr;           // bit vector for corrupted
  
public:
  Dispute();
  ~Dispute();

  void reset(int N);
  void addDispPairs(int p1, int p2);
  void addCorrParty(int p);

  int tSet(vector<int>& TSet);
  int nonTSet(vector<int>& NonTSet);
  int tAndNonTSet(vector<int>& TSet, vector<int>& NonTSet);
  void corrSet(vector<int>& corrs);
  int dispAndNonDispSet(int p, vector<int>& DispSet, vector<int>& NonDispSet);
  int twistSet(int p, int target, vector<int>& TwistSet);
  
  bool rlyeeVecs(int p, vector<int>& relay, vector<bool>& relayerMask,
                 vector<vector<int>>& relayerLoad);
  bool rlyerVecs(int p, vector<bool>& relayeeMask,
                 vector<vector<int>>& relayeeLoad);
  int tMask(vector<bool>& TMask);
  
  void corrMask(vector<bool>& cMask) { cMask = _corr;}
  void nonDispMask(int p, vector<bool> &dMask) {
    dMask.resize(_N);
    for (int i=0; i<_N; i++) { dMask[i] = !_disp[p][i]; }
  }
  bool hasCorrupt() {return (_nCorr > 0);}
  bool hasDisp() {return _hasDisp;};
  bool hasDisp(int p) {return (_dispSize[p] > 0);};
  bool isCorrupt(int p) {return _corr[p];}
  bool isDisp(int myId, int p) {return _disp[myId][p];}
  void setNewDisp(bool newDisp) {_newDisp = newDisp;}
  bool getNewDisp() {return _newDisp;}

  // void getDisp(int p, vector<bool>& disp) {disp = _disp[p];};
};


#endif
