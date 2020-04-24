// author: Hanjun Li <lihanjun1212@gmail.com>
#ifndef DISPUTE_H_
#define DISPUTE_H_

#define NDEBUG

#include <assert.h>
#include <vector>
#include <unordered_set>

using namespace std;

class Dispute{
  // dummy class

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

  void tMaskP(int p, vector<bool>& TMask);
  void tMaskPRev(int p, vector<bool>& TMaskRev);
  void tAndNonTSetP(int p, vector<int>& TSet, vector<int>& NonTSet);
  void reset(int N);
};


#endif
