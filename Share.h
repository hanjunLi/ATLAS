// author: Hanjun Li <lihanjun1212@gmail.com>
#ifndef SHARE_H_
#define SHARE_H_

#include <stdlib.h>
#include <algorithm>

#include <libscapi/include/primitives/Matrix.hpp>
#include <vector>
#include <iostream>
#include <libscapi/include/primitives/Mersenne.hpp>
#include <libscapi/include/infra/Common.hpp>

template <class FieldType>
class Share{
public:
  FieldType _value;
  vector<unsigned int> _idxs;
  vector<FieldType> _coeffs;
  // bool degT = true;
  
  // Share();
  // ~Share();

  Share<FieldType>& operator+=(const Share<FieldType>& s2) {
    _value += s2._value;
    vector<unsigned int> newIdxs;
    vector<FieldType> newCoeffs;

    int nTerms1 = _idxs.size();
    int nTerms2 = s2._idxs.size();
    int pos2 = 0;
    for (int i=0; i<nTerms1; i++) {
      unsigned int idx1 = _idxs[i];
      FieldType coeff1 = _coeffs[i];
      while (pos2 < nTerms2 && s2._idxs[pos2] < idx1) {
        newIdxs.push_back(s2._idxs[pos2]);
        newCoeffs.push_back(s2._coeffs[pos2]);
        pos2++;
      }
      newIdxs.push_back(idx1);
      if (pos2 < nTerms2 && s2._idxs[pos2] == idx1) { // merge
        newCoeffs.push_back(coeff1 + s2._coeffs[pos2]);
      } else {
        newCoeffs.push_back(coeff1);
      }
    }

    _idxs = newIdxs;
    _coeffs = newCoeffs;
    return *this;
  }

  Share<FieldType>& operator-=(const Share<FieldType>& s2) {
    _value = _value - s2._value;
    vector<unsigned int> newIdxs;
    vector<FieldType> newCoeffs;

    int nTerms1 = _idxs.size();
    int nTerms2 = s2._idxs.size();
    int pos2 = 0;
    for (int i=0; i<nTerms1; i++) {
      unsigned int idx1 = _idxs[i];
      FieldType coeff1 = _coeffs[i];
      while (pos2 < nTerms2 && s2._idxs[pos2] < idx1) {
        newIdxs.push_back(s2._idxs[pos2]);
        // CHECK: FieldType(0) works for Zp or ?
        newCoeffs.push_back(FieldType(0) - s2._coeffs[pos2]);
        pos2++;
      }
      newIdxs.push_back(idx1);
      if (pos2 < nTerms2 && s2._idxs[pos2] == idx1) { // merge
        newCoeffs.push_back(coeff1 - s2._coeffs[pos2]);
      } else {
        newCoeffs.push_back(coeff1);
      }
    }

    _idxs = newIdxs;
    _coeffs = newCoeffs;
    return *this;
  }

  const Share<FieldType> operator+(const Share<FieldType>& s2) const {
    return Share<FieldType>(*this) += s2;
  }

  const Share<FieldType> operator-(const Share<FieldType>& s2) const {
    return Share<FieldType>(*this) -= s2;
  }


  void multByConst(FieldType c) {
    _value = _value * c;
    int nTerms = _idxs.size();
    for (int i=0; i<nTerms; i++) {
      _coeffs[i] *= c;
    }
  }
};

template <class FieldType>
class rTranscript {             // for refresh
public:
  // x, r
  Share<FieldType> _x, _r;
  // xp, e, o
  FieldType _xp, _e, _o;

  // void init(); default for members should just work

  vector<FieldType> getValues() {
    // x, xp, r, e, o
    vector<FieldType> values { _x._value, _xp, _r._value, _e, _o };
    return values;
  }

  void addTo(Share<FieldType>& x, Share<FieldType>& r, FieldType o) {
    _x += x; _r += r;
    _xp += (x._value - o);
    _e += (x._value + r._value);
    _o += o;
  }
};

template <class FieldType>
class mTranscript {
public:
  Share<FieldType> _r, _r2;
  FieldType _e, _e2, _z;

  vector<FieldType> getValues() {
    // r, r2, e, e2, z
    vector<FieldType> values { _r._value, _r2._value, _e, _e2, _z};
    return values;
  }

  void addTo(Share<FieldType>& r, Share<FieldType>& r2, FieldType e, FieldType xy) {
    _r += r; _r2 += r2; _e += e;
    _e2 += (xy - r2._value);
    _z += (e + r._value);
  }

  mTranscript<FieldType>& operator-=(const mTranscript<FieldType>& t2) {
    _r -= t2._r; _r2 -= t2._r2;
    _e = _e - t2._e; _e2 = _e2 - t2._e2; _z = _z - t2._z;
    return *this;
  }

  mTranscript<FieldType>& operator+=(const mTranscript<FieldType>& t2) {
    _r += t2._r; _r2 += t2._r2;
    _e += t2._e; _e2 += t2._e2; _z += t2._z;
    return *this;
  }
  
  const mTranscript<FieldType>
  operator-(const mTranscript<FieldType>& t2) const {
    return mTranscript<FieldType>(*this) -= t2;
  }

  const mTranscript<FieldType>
  operator+(const mTranscript<FieldType>& t2) const {
    return mTranscript<FieldType>(*this) += t2;
  }
};


#endif
