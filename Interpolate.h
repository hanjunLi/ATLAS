// author: Hanjun Li <hanjunl@andrew.cmu.edu>

#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <stdlib.h>

#include <libscapi/include/primitives/Matrix.hpp>
#include <libscapi/include/cryptoInfra/Protocol.hpp>
// #include <libscapi/include/circuits/ArithmeticCircuit.hpp>
#include <libscapi/include/infra/Measurement.hpp>
#include <vector>
#include <bitset>
#include <iostream>
#include <fstream>
#include <chrono>
#include <libscapi/include/primitives/Mersenne.hpp>
#include "ProtocolTimer.h"
#include <libscapi/include/comm/MPCCommunication.hpp>
#include <libscapi/include/infra/Common.hpp>
// #include <libscapi/include/primitives/Prg.hpp>
// #include "HashEncrypt.h"
#include <emmintrin.h>
// TODO: clean up unused includes / functions


using namespace std;

template <class FieldType>
class Interpolate{

private:
  void trimZeroes(vector<FieldType>& a);

  void addToPolynomial(vector<FieldType>& p1, // input
                       vector<FieldType>& p2);

  void multToPolynomial(vector<FieldType>& p1,
                        vector<FieldType>& p2);

  void scaleToPolynomial(FieldType c,
                         vector<FieldType>& p);

  void dividePolynomial(vector<FieldType>& p1,
                        vector<FieldType>& p2,
                        vector<FieldType>& q, // quotient
                        vector<FieldType>& r);
public:
  Interpolate();
  ~Interpolate();

  // debugging
  void printPolynomial(vector<FieldType>& p);

  // independent version: using supplied x values
  FieldType evalPolynomial(FieldType x,
                           vector<FieldType>& polynomial);

  void interpolate(vector<FieldType>& x, // input 
                   vector<FieldType>& y, // input
                   vector<FieldType>& polynomial);
};
  
// ---------------- implementations ----------------
template <class FieldType>
Interpolate<FieldType>::Interpolate(){
}

template <class FieldType>
Interpolate<FieldType>::~Interpolate(){
}


template<class FieldType> 
void Interpolate<FieldType>::
trimZeroes(vector<FieldType>& a) {
  int i = a.size();
  FieldType Zero = FieldType(0);
  while(i > 0 && a[i-1] == Zero) { i--; }
  a.resize(i);
  return;
}

template <class FieldType>
void Interpolate<FieldType>::
addToPolynomial(vector<FieldType>& p1, // input
                vector<FieldType>& p2){ // both input and output
  // store result in p2
  int p1_size = p1.size();
  int p2_size = p2.size();
  if(p1_size > p2_size){
    p2.resize( p1_size, FieldType(0) );
  }

  for(int i=0; i<p1_size; i++){
    p2[i] += p1[i];
  }
  trimZeroes(p2);
  return;
}

template <class FieldType>
void Interpolate<FieldType>::
multToPolynomial(vector<FieldType>& p1,
                 vector<FieldType>& p2){

  int p1_size = p1.size();
  int p2_size = p2.size();
  vector<FieldType> tempProduct(p1_size + p2_size - 1, FieldType(0));

  for(int i=0; i<p1_size; i++){
    if(p1[i] == FieldType(0)){
      continue;
    }
    for(int j=0; j<p2_size; j++){
      tempProduct[i+j] += p1[i] * p2[j];
    }
  }
    
  trimZeroes(tempProduct);
  p2 = tempProduct;
  return;
}

template <class FieldType>
void Interpolate<FieldType>::
scaleToPolynomial(FieldType c,
                  vector<FieldType>& p){
  int p_deg = p.size();
  for(int i=0; i<p_deg; i++){
    p[i] *= c;
  }
  return;
}

template <class FieldType>
void Interpolate<FieldType>::
dividePolynomial(vector<FieldType>& p1,
                 vector<FieldType>& p2,
                 vector<FieldType>& q, // quotient
                 vector<FieldType>& r){ // remainder
  r = p1;
  int p1Size = p1.size();
  int p2Size = p2.size();
  if(p1Size < p2Size){
    q.resize(0);
    r = p1;
    return;
  }
  if (p2Size == 0) {
    cerr << "ECC: dividing by zero (polynomail)" << endl;
    abort();
  }

  int qSize = p1Size - p2Size +1;
  q.resize(qSize);

  for(int i = p1Size - p2Size; i>=0; i--){
    FieldType topCoeff = r[ p2Size + i -1] / p2[p2Size-1];

    q[i] = topCoeff;
    vector<FieldType> xi(i+1, FieldType(0));
    xi[i] = xi[i] - topCoeff;
    vector<FieldType> p2Tmp = p2;
    
    // r -= topCoeff * p2 * x^i
    multToPolynomial(xi, p2Tmp);
    addToPolynomial(p2Tmp, r);
  }
  trimZeroes(q);
  trimZeroes(r);
  return;
}

template <class FieldType>
FieldType Interpolate<FieldType>::
evalPolynomial(FieldType x, 
               vector<FieldType>& polynomial){
  int degree = polynomial.size() - 1;
  // assert(degree >= 0);
  FieldType result = polynomial[0];
  FieldType x_value = FieldType(1);
  for(int i=0; i<degree; i++){
    x_value *= x;
    result += polynomial[i+1] * x_value;
  }
  
  return result;
}


template <class FieldType>
void Interpolate<FieldType>::
interpolate(vector<FieldType>& x, // input 
            vector<FieldType>& y, // input
            vector<FieldType>& polynomial){

  // assert(y.size() <= x.size());
  int nPoints = y.size();
  int degree = nPoints - 1;
  vector<FieldType> result(nPoints, FieldType(0));

  // ---- O(n^2) to compute all numerators ----
    vector< vector<FieldType> > numerator_before_i(nPoints);
  vector< vector<FieldType> > numerator_skip_i(nPoints);
  // fill-in numerator_before_i from left to right
  numerator_before_i[0] = vector<FieldType>(1, FieldType(1));
  for(int i=1; i<nPoints; i++){
    numerator_before_i[i] = vector<FieldType>( 2, FieldType(1) );
    numerator_before_i[i][0] = FieldType(0) - x[i-1];
    multToPolynomial( numerator_before_i[i-1], numerator_before_i[i] );
  }
  // fill-in numerator_after_i from right to left
  numerator_skip_i[nPoints-1] = vector<FieldType>(1, FieldType(1));
  for(int i=nPoints-2; i>=0; i--){
    numerator_skip_i[i] = vector<FieldType>(2, FieldType(1));
    numerator_skip_i[i][0] = FieldType(0) - x[i+1];
    multToPolynomial( numerator_skip_i[i+1], numerator_skip_i[i] );
  }
  // multiply before_i and after_i to get skip_i
  for(int i=0; i<nPoints; i++){
    multToPolynomial(numerator_before_i[i], numerator_skip_i[i]);
  }

  // ---- O(n^2) to compute all factors ----
  vector<FieldType> factor(nPoints);
  for(int i=0; i<nPoints; i++){
    FieldType denom = evalPolynomial( x[i], numerator_skip_i[i] );
    factor[i] = y[i] / denom;
  }
  
  // ---- O(n^2) add to get result ----
  for(int i=0; i<nPoints; i++){
    scaleToPolynomial( factor[i], numerator_skip_i[i] );
    addToPolynomial(numerator_skip_i[i], result);
  }

  polynomial = result;
  return;
}

template <class FieldType>
void Interpolate<FieldType>::
printPolynomial(vector<FieldType>& p){
  int deg = 0;
  for(auto e : p){
    cout << e << "*x^" << (deg++) << " ";
  }
  cout << endl;
  return;
}

#endif
