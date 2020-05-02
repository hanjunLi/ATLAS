#ifndef PROTOCOLPARTY_H_
#define PROTOCOLPARTY_H_

#include <stdlib.h>
#include "HashEncrypt.h"
#include "ProtocolTimer.h"
#include <bitset>
#include <chrono>
#include <emmintrin.h>
#include <fstream>
#include <iostream>
#include <libscapi/include/circuits/ArithmeticCircuit.hpp>
#include <libscapi/include/comm/MPCCommunication.hpp>
#include <libscapi/include/cryptoInfra/Protocol.hpp>
#include <libscapi/include/infra/Common.hpp>
#include <libscapi/include/infra/Measurement.hpp>
#include <libscapi/include/primitives/Matrix.hpp>
#include <libscapi/include/primitives/Mersenne.hpp>
#include <thread>
#include <vector>

#include "Interpolate.h"
#include "Dispute.h"
#include "Communication.h"
#include <cmath>

#define flag_print true
#define flag_print_output false

using namespace std;
using namespace std::chrono;

template <class FieldType>
class ProtocolParty : public Protocol, public HonestMajority, MultiParty {

private:
  // -- polynomial functionalities
  Interpolate<FieldType> interp; // evaluate (O(n)), interpolate polynomials (O(n^2))

  // -- global const
  int numThreads = 4;      // TODO: add as main arguments later
  int _K = 7;              // interpolation degree <=> 'shrink' factor

  // -- global variables
  int iteration;                       // current iteration number
  int currentCircuitLayer = 0;         // current circuit layer
  vector<FieldType> y_for_interpolate; // size = N

  // -- additional classes and data  
  Dispute _disp;
  HIM<FieldType> _mat_to_T, _mat_n_to_n;;
  FieldType _zero;
  int _fieldByteSize;
  vector<vector<bool>> _TMasks;

  // -- set in constructor
  int times;                // number of times to run the run function
  Measurement *timer;
  TemplateField<FieldType> *_field; // used to call some field functions
  ProtocolTimer *protocolTimer;
  // N, T - number of parties, malicious
  int _N, _T, keySize, verifyIterations, _myId;
  string s;                     // string of my_partyId
  string inputsFile, outputFile;
  ArithmeticCircuit circuit;
  // M - number of gates
  int M, numOfInputGates, numOfOutputGates, numOfMultGates;
  vector<long> myInputs;

  // -- set in initialization phase
  vector<FieldType> beta;      // a single zero element
  vector<FieldType> alpha;     // N distinct non-zero field elements
  vector<FieldType> _alpha_k;  // K distinct non-zero field elements
  vector<FieldType> _alpha_2k; // 2K-1 distinct non-zero field elements
  HIM<FieldType> matrix_for_interpolate;
  HIM<FieldType> matrix_for_t;  // interpolate first t+1 points to rest N-t-1
  HIM<FieldType> matrix_for_2t; // interpolate first 2t+1 tpoints to rest N-2t-1
  HIM<FieldType> matrix_for_k;  // interpolate first k points to rest 2k-1 - k
  VDM<FieldType> matrix_vand;
  VDMTranspose<FieldType> matrix_vand_transpose;

  // -- filled during offline-preparation phase
  vector<FieldType> _singleSharesArray; // used as pads and rand coins
  int _singleSharesOffset = 0;          // next random share to use
  vector<FieldType> _doubleSharesArray; // used for DN mult
  int _doubleSharesOffset = 0;    // next random double share to use
  vector<int> _doubleKingIds;     // for T-wise
  vector<FieldType> gateShareArr; // my share of each gate (1 share each)

  __inline__
  void intToElements(vector<int>& intVec, vector<FieldType>& elemVec) {
    int nInts = intVec.size();
    elemVec.resize(nInts);
    for (int i=0; i<nInts; i++) {
      elemVec[i] = _field->GetElement(intVec[i]+1);
    }
  }

  __inline__
  void matForRefresh(int pid, HIM<FieldType>& mat_to_T) {
    vector<int> TSet, NonTSet;
    _disp.tAndNonTSetP(pid, TSet, NonTSet);
    vector<FieldType> TElems, NonTElems;
    intToElements(TSet, TElems);
    intToElements(NonTSet, NonTElems);
    vector<FieldType> alpha_nonT = NonTElems;
    alpha_nonT.push_back(_zero);
    alpha_nonT.insert(alpha_nonT.end(), TElems.begin(), TElems.end());
    mat_to_T.allocate(_T+1, _N - _T, _field);
    mat_to_T.InitHIMVectorAndsizes(alpha_nonT, _N - _T, _T+1);
  }

  __inline__
  void make2TShares(int nShares, vector<FieldType>& secrets,
                    vector<vector<FieldType>>& fullSharesVec) {
    secrets.resize(nShares);
    fullSharesVec.clear();
    fullSharesVec.resize(_N ,vector<FieldType>(nShares));
    vector<FieldType> doubleShares(_N);
    for (int i=0; i<nShares; i++) {
      for (int j=0; j<_N; j++) {
        doubleShares[j] = _field->Random();
        fullSharesVec[j][i] = doubleShares[j];
      }
      secrets[i] = interpolate_to_zero(doubleShares);
    }
  }

public:
  bool hasOffline() { return true; }
  bool hasOnline() override { return true; }
  ~ProtocolParty();

  // -- from new_protocol_god
  void batchMultSingle(vector<FieldType>& e2Shares, vector<FieldType>& e1Shares);
  void makeRandShares(int nRands, vector<FieldType> &randShares);
  void makeRandDoubleShares(int nRands, vector<FieldType> &randDoubleShares, vector<int>& kingIds);

  // -- protocol functions
  ProtocolParty(int argc, char *argv[]);
  void readMyInputs();
  void initializationPhase();

  void run() override;
  void runOffline() override;
  bool preparationPhase();  

  void runOnline() override;
  void inputPhase();
  void inputVerification(vector<FieldType> &inputShares);

  void computationPhase();
  int processNotMult();
  int processMultiplications(int lastMultGate);
  int processMultDN(int indexInRandomArray);

  void verificationPhase();
  // compute vector products w/ DN Mult
  void interpolatePolyVec(vector<FieldType>& aShares,
                          vector<vector<FieldType>>& AShares, int groupSize);
  void evalPolyVecAt(vector<vector<FieldType>>& AShares,
                     vector<FieldType>& aShares, FieldType point);
  void evalPolyVec(vector<vector<FieldType>>& AShares,
                   vector<FieldType>& aShares, int nPoints, int offset);
  void DNMultVec(vector<FieldType>& a, vector<FieldType>& b,
                 vector<FieldType>& c, int groupSize);
  void buildPolyVec(vector<FieldType>& aShares, vector<FieldType>& bShares,
                    FieldType& cShare, vector<FieldType>& dShares, int groupSize,
                    vector<vector<FieldType>>& AShares,
                    vector<vector<FieldType>>& BShares,
                    vector<FieldType>& CShare);
  void buildPolyVecIndWorker(vector<FieldType>& aShares, vector<FieldType>& bShares,
                             vector< vector<FieldType> >& ASharesEval,
                             vector< vector<FieldType> >& BSharesEval,
                             int groupSize, int threadId);
  void buildPolyVecInd(vector<FieldType>& aShares, vector<FieldType>& bShares,
                       vector<FieldType>& dShares, int groupSize);
  void compressVerifyVec(vector<FieldType>& aShares,
                         vector<FieldType>& bShares, FieldType& cShare);
  void verifyVec(vector<FieldType>& aShares,
                 vector<FieldType>& bShares, FieldType& cShare);
  void outputPhase();

  // -- secret sharing functionalities
  void openShare(int numOfRandomShares, vector<FieldType> &Shares,
                 vector<FieldType> &secrets);
  FieldType randomCoin();

  // --  helper functions
  void getRandomShares(int numOfRandoms,
                       vector<FieldType> &randomElementsToFill);
  bool checkConsistency(vector<FieldType> &x, int d);
  // locally interpolate received sharing
  FieldType reconstructShare(vector<FieldType> &x, int d);
  FieldType interpolate_to_zero(vector<FieldType> &x);

  // -- external access
  int getN() {return _N;}
  int getT() {return _T;}
  int getMyId() {return _myId;}
  TemplateField<FieldType>* getField() {return _field;}
  bool preparationPhase(int nSingle, int nDouble);
  Communication _comm;

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
        for (int j=0; j<_N-_T-1; j++) {
          int pnt = NonTSet[j];
          nonTShares[j] = _field->Random();
          fullSharesVec[pnt][i] = nonTShares[j];
        }
      } // else: left nonT shares to zero
      _mat_to_T.MatrixMult(nonTShares, TShares);
      for (int j=0; j<_T+1; j++) {
        int pt = TSet[j];
        fullSharesVec[pt][i] = TShares[j];
      }
    }
  }
};


template <class FieldType>
void ProtocolParty<FieldType>::
makeRandShares(int nRands, vector<FieldType> &randShares) {
  int nBuckets = (nRands / (_N - _T)) + 1;
  randShares.resize(nBuckets * (_N - _T));

  vector<FieldType> secrets(nBuckets);
  for (int i=0; i<nBuckets; i++) {
    secrets[i] = _field->Random();
  }
  if(flag_print)
	  cout<<"An example of secret:"<<secrets[0]<<endl;
  vector<vector<FieldType>> randSharesAll(_N, vector<FieldType>(nBuckets));
  makeTShares(secrets, randSharesAll);

  vector<vector<byte>> sendBufs(_N, vector<byte>(nBuckets*_fieldByteSize));
  for (int i=0; i<_N; i++) {
    _field->elementVectorToByteVector(randSharesAll[i], sendBufs[i]);
  }
  vector<vector<byte>> recBufs(_N, vector<byte>(nBuckets*_fieldByteSize, 0));
  _comm.allToAll(sendBufs, recBufs);

  int count = 0;
  vector<FieldType> bufferIn(_N), bufferOut(_N);
  for (int i=0; i<nBuckets; i++) { // improve: switch loop order
    for (int j=0; j<_N; j++) {
      bufferIn[j] =
        _field->bytesToElement(recBufs[j].data() + i * _fieldByteSize);
    }
    matrix_vand_transpose.MatrixMult(bufferIn, bufferOut, _N - _T);
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
  kingIds.resize(nBuckets * (_N-_T));

  vector<FieldType> secrets;
  vector<vector<FieldType>> dSharesAll, sharesAll, randSharePairsAll;
  make2TShares(nBuckets, secrets, dSharesAll);
  makeTShares(secrets, sharesAll);
  // mix single and double shares
  randSharePairsAll.resize(_N, vector<FieldType>(nBuckets*2));
  for (int j=0; j<_N; j++) {
    for (int i = 0; i < nBuckets; i++) {
      randSharePairsAll[j][i*2] = sharesAll[j][i]; // even slots: single share
      randSharePairsAll[j][i*2+1] = dSharesAll[j][i]; // odd slots: double Share
    }
  }
  // mix single and double shares
  vector<vector<byte>> sendBufs(_N, vector<byte>(nBuckets*_fieldByteSize*2));
  for (int i=0; i<_N; i++) {
    _field->elementVectorToByteVector(randSharePairsAll[i], sendBufs[i]);
  }
  vector<vector<byte>> recBufs(_N, vector<byte>(nBuckets*_fieldByteSize*2, 0));
  _comm.allToAll(sendBufs, recBufs);
  int count = 0;
  int king = 0;
  vector<FieldType> bufferIn(_N), bufferIn2(_N), bufferOut(_N), bufferOut2(_N);
  for (int i=0; i<nBuckets; i++) {
    for (int j=0; j<_N; j++) {
      bufferIn[j] =
        _field->bytesToElement(recBufs[j].data() + 2*i*_fieldByteSize);
      bufferIn2[j] =
        _field->bytesToElement(recBufs[j].data() + (2*i+1)*_fieldByteSize);
    }
    _mat_n_to_n.MatrixMult(bufferIn, bufferOut);
    _mat_n_to_n.MatrixMult(bufferIn2, bufferOut2);
    for (int j=0; j<(_N-_T); j++) {
      randSharePairs[count*2] = bufferOut[j];
      randSharePairs[count*2+1] = bufferOut2[j];
      kingIds[count++] = (king++) % _N;
    }
  }
}

template <class FieldType>
ProtocolParty<FieldType>::ProtocolParty(int argc, char *argv[])
    : Protocol("MPCHonestMajorityNoTriples", argc, argv) {
  // -- initialize protocol parameters
  this->times = stoi(
      this->getParser().getValueByKey(arguments, "internalIterationsNumber"));
  vector<string> subTaskNames{
      "Offline",      "preparationPhase",  "Online",     "inputPhase",
      "ComputePhase", "VerificationPhase", "outputPhase"};
  this->timer = new Measurement(*this, subTaskNames);
  string fieldType = this->getParser().getValueByKey(arguments, "fieldType");
  if (fieldType.compare("ZpMersenne31") == 0) {
    _field = new TemplateField<FieldType>(2147483647);
  } else if (fieldType.compare("ZpMersenne61") == 0) {
    _field = new TemplateField<FieldType>(0);
  } else if (fieldType.compare("ZpKaratsuba") == 0) {
    _field = new TemplateField<FieldType>(0);
  } else if (fieldType.compare("GF2E") == 0) {
    _field = new TemplateField<FieldType>(8);
  } else if (fieldType.compare("Zp") == 0) {
    _field = new TemplateField<FieldType>((1UL << 61) -1);
  } else if (fieldType.compare("ZpMersenne127")==0)
  {
	  _field = new TemplateField<FieldType>(0);
  }
  string circuitFile =
    this->getParser().getValueByKey(arguments, "circuitFile");
  string outputTimerFileName =
    circuitFile + "Times" + to_string(_myId) + fieldType + ".csv";
  this->protocolTimer = new ProtocolTimer(times, outputTimerFileName);
  _N = stoi(this->getParser().getValueByKey(arguments, "partiesNumber"));
  _T = (_N - 1)/ 2;
  this->keySize = 16 / _field->getElementSizeInBytes() + 1;
  this->verifyIterations =
    (5 + _field->getElementSizeInBytes() - 1) / _field->getElementSizeInBytes();
  _myId = stoi(this->getParser().getValueByKey(arguments, "partyID"));
  this->s = to_string(_myId);
  this->inputsFile = this->getParser().getValueByKey(arguments, "inputFile");
  this->outputFile = this->getParser().getValueByKey(arguments, "outputFile");
  this->circuit.readCircuit(circuitFile.c_str());
  this->circuit.reArrangeCircuit();
  this->M = circuit.getNrOfGates();
  this->numOfInputGates = circuit.getNrOfInputGates();
  this->numOfOutputGates = circuit.getNrOfOutputGates();
  this->numOfMultGates = circuit.getNrOfMultiplicationGates();
  this->myInputs.resize(numOfInputGates);  

  _disp.reset(_N);
  string partiesFile =
    getParser().getValueByKey(arguments, "partiesFile");
  _comm.reset(_N, _myId, numThreads, partiesFile, &_disp);
  _zero = *(_field->GetZero());
  _fieldByteSize = _field->getElementSizeInBytes();
  _TMasks.resize(_N);
  for(int i=0; i<_N; i++) {
    _disp.tMaskP(i, _TMasks[i]);
  }

  // -- initialize phase
  initializationPhase();
}

template <class FieldType> void ProtocolParty<FieldType>::readMyInputs() {

  ifstream myfile;
  long input;
  int i = 0;
  myfile.open(this->inputsFile);
  do {
    myfile >> input;

    if (i >= myInputs.size()){
      cout << "inputs file " << this->inputsFile <<endl;
      cout << "have more than " << (i+1) << " inputs!" << endl;
      abort();
    }

    myInputs[i++] = input;
  } while (!(myfile.eof()));
  myfile.close();
}

template <class FieldType> void ProtocolParty<FieldType>::run() {

  for (iteration = 0; iteration < times; iteration++) {

    auto t1start = high_resolution_clock::now();
    timer->startSubTask("Offline", iteration);
    runOffline();
    timer->endSubTask("Offline", iteration);
    timer->startSubTask("Online", iteration);
    runOnline();
    timer->endSubTask("Online", iteration);
    auto t2end = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2end - t1start).count();
    protocolTimer->totalTimeArr[iteration] = duration;
  }
}

template <class FieldType> void ProtocolParty<FieldType>::runOffline() {
  auto t1 = high_resolution_clock::now();
  timer->startSubTask("preparationPhase", iteration);
  if (preparationPhase() == false) {
    if (flag_print) {
      cout << "preparationPhase failed" << '\n';
    }
    return;
  } else {
    if (flag_print) {
      cout << "finish Preparation Phase" << '\n';
    }
  }
  timer->endSubTask("preparationPhase", iteration);
  auto t2 = high_resolution_clock::now();

  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  protocolTimer->preparationPhaseArr[iteration] = duration;
}

template <class FieldType> void ProtocolParty<FieldType>::runOnline() {

  auto t1 = high_resolution_clock::now();
  timer->startSubTask("inputPhase", iteration);
  inputPhase();
  timer->endSubTask("inputPhase", iteration);
  auto t2 = high_resolution_clock::now();

  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  protocolTimer->inputPreparationArr[iteration] = duration;

  t1 = high_resolution_clock::now();
  timer->startSubTask("ComputePhase", iteration);
  computationPhase();
  timer->endSubTask("ComputePhase", iteration);
  t2 = high_resolution_clock::now();

  duration = duration_cast<milliseconds>(t2 - t1).count();
  protocolTimer->computationPhaseArr[iteration] = duration;

  t1 = high_resolution_clock::now();
  timer->startSubTask("VerificationPhase", iteration);
  verificationPhase();
  timer->endSubTask("VerificationPhase", iteration);
  t2 = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(t2 - t1).count();
  protocolTimer->verificationPhaseArr[iteration] = duration;

  t1 = high_resolution_clock::now();
  timer->startSubTask("outputPhase", iteration);
  outputPhase();
  timer->endSubTask("outputPhase", iteration);
  t2 = high_resolution_clock::now();

  duration = duration_cast<milliseconds>(t2 - t1).count();
  protocolTimer->outputPhaseArr[iteration] = duration;
}

template <class FieldType>
void ProtocolParty<FieldType>::computationPhase() {
  int countNumMult = _doubleSharesOffset;
  int numOfLayers = circuit.getLayers().size();
  for (int i = 0; i < numOfLayers - 1; i++) {
    currentCircuitLayer = i;
    // send the index of the current mult gate
    processNotMult();
    countNumMult += processMultiplications(countNumMult);
  }
  _doubleSharesOffset = countNumMult;
}

template <class FieldType>
void ProtocolParty<FieldType>::verificationPhase() {
  vector<FieldType> aShares(this->numOfMultGates);
  vector<FieldType> bShares(this->numOfMultGates);
  FieldType cShare = *(_field->GetZero());

  // -- open a random share lambda
  FieldType lambda = randomCoin();
  // -- gather all mult triples to be verified
  //    and combine them into a dot product
  int idx = 0;
  FieldType lambdaI = lambda;
  for (auto& gate : this->circuit.getGates()) {
    if(gate.gateType == MULT) {
      aShares[idx] = gateShareArr[gate.input1] * lambdaI;
      bShares[idx] = gateShareArr[gate.input2];
      cShare += gateShareArr[gate.output] * lambdaI;
      idx++;
      lambdaI *= lambda;
    }
  }
  compressVerifyVec(aShares, bShares, cShare);
}

template <class FieldType>
void ProtocolParty<FieldType>::
buildPolyVec(vector<FieldType>& aShares, vector<FieldType>& bShares,
             FieldType& cShare, vector<FieldType>& dShares, int groupSize,
             vector<vector<FieldType>>& AShares,
             vector<vector<FieldType>>& BShares, vector<FieldType>& CShare) {

  // batchDegree is _K except in the base case
  int batchDegree = (aShares.size() + groupSize - 1) / groupSize;
  AShares.resize(groupSize);
  BShares.resize(groupSize);

  // -- interploate groupSize polynomials of degree K-1
  // interpolate A[i], B[i]
  interpolatePolyVec(aShares, AShares, groupSize);
  interpolatePolyVec(bShares, BShares, groupSize);

  // evaluate A(k), ..., A(2k-2), B(k), ..., B(2k-2)
  vector<FieldType> ASharesEval;
  vector<FieldType> BSharesEval;
  evalPolyVec(AShares, ASharesEval, batchDegree-1, batchDegree);
  evalPolyVec(BShares, BSharesEval, batchDegree-1, batchDegree);

  // -- dot product batch verification
  // build dShares k .. 2k - 2 and interpolate C
  vector<FieldType> dSharesRest;
  DNMultVec(ASharesEval, BSharesEval, dSharesRest, groupSize);
  dShares.insert(dShares.end(), dSharesRest.begin(), dSharesRest.end());
  interp.interpolate(_alpha_2k, dShares, CShare);
}

template <class FieldType>
void ProtocolParty<FieldType>::
buildPolyVecIndWorker(vector<FieldType>& aShares, vector<FieldType>& bShares,
                      vector< vector<FieldType> >& ASharesEval,
                      vector< vector<FieldType> >& BSharesEval,
                      int groupSize, int threadId) {
  ZZ_p::init(ZZ((1UL << 61) -1)); // hack to init zzp

  int totalLength = aShares.size();
  vector<FieldType> ySharesA;   // tmp vector
  vector<FieldType> ySharesB;   // tmp vector

  for (int i=threadId; i<groupSize; i+= numThreads) {
    ySharesA.clear();
    ySharesB.clear();
    for (int j=i; j<totalLength; j+=groupSize) {
      ySharesA.push_back(aShares[j]);
      ySharesB.push_back(bShares[j]);
    }
    matrix_for_k.MatrixMult(ySharesA, ASharesEval[i]);
    matrix_for_k.MatrixMult(ySharesB, BSharesEval[i]);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
buildPolyVecInd(vector<FieldType>& aShares, vector<FieldType>& bShares,
                vector<FieldType>& dShares, int groupSize) {

  // batchDegree is _K except in the base case
  int batchDegree = (aShares.size() + groupSize - 1) / groupSize;

  // interpolate A[i], B[i] and evaluate at A(k), ..., A(2k-2), B(k), ..., B(2k-2)
  vector< vector<FieldType> > ASharesEval(groupSize, vector<FieldType>(_K-1));
  vector< vector<FieldType> > BSharesEval(groupSize, vector<FieldType>(_K-1));

  // Do multi-threading here
  vector<thread> threads;
  if (numThreads > 1) {
    threads.resize(numThreads-1);
    for (int t=0; t < numThreads-1; t++) {
      // bring up thread workers
      threads[t] = thread(&ProtocolParty::buildPolyVecIndWorker, this,
                          ref(aShares), ref(bShares),
                          ref(ASharesEval), ref(BSharesEval), groupSize, t+1);
    }
  }

  int totalLength = aShares.size();
  vector<FieldType> ySharesA;   // tmp vector
  vector<FieldType> ySharesB;   // tmp vector

  for (int i=0; i<groupSize; i+= numThreads) {
    ySharesA.clear();
    ySharesB.clear();
    for (int j=i; j<totalLength; j+=groupSize) {
      ySharesA.push_back(aShares[j]);
      ySharesB.push_back(bShares[j]);
    }
    matrix_for_k.MatrixMult(ySharesA, ASharesEval[i]);
    matrix_for_k.MatrixMult(ySharesB, BSharesEval[i]);
  }

  if (numThreads > 1) {
    // join thread workers
    for (int t=0; t < numThreads-1; t++) {
      threads[t].join();
    }
  }


  // vector-mult to get the rest of dShares
  vector<FieldType> ASharesEvalFlat(groupSize * (_K - 1));
  vector<FieldType> BSharesEvalFlat(groupSize * (_K - 1));
  for (int i=0; i<_K-1; i++) {
    // point i
    for (int j=0; j<groupSize; j++) {
      // poly j
      ASharesEvalFlat[i * groupSize + j] = ASharesEval[j][i];
      BSharesEvalFlat[i * groupSize + j] = BSharesEval[j][i];
    }
  }

  vector<FieldType> dSharesRest;
  DNMultVec(ASharesEvalFlat, BSharesEvalFlat, dSharesRest, groupSize);

  // append to passed-in dShares
  dShares.insert(dShares.end(), dSharesRest.begin(), dSharesRest.end());
}

// -- keep compressing the dot product to be verified until
//    it's small enough
template <class FieldType>
void ProtocolParty<FieldType>::
compressVerifyVec(vector<FieldType>& aShares, vector<FieldType>& bShares,
                  FieldType& cShare){
  int totalLength = aShares.size();
  if (totalLength <= _K) { // base case:
    verifyVec(aShares, bShares, cShare);
    return;
  } // else : recursive case:
  // -- divide into K groups
  int groupSize = (totalLength + _K -1 )/ _K;
  aShares.resize( groupSize * _K, _field->GetElement(0) );
  bShares.resize( groupSize * _K, _field->GetElement(0) );
  totalLength = groupSize * _K;

  // -- one DN mult for each group i to compute
  // build dShares 0 .. k-2
  vector<FieldType> dShares;
  DNMultVec(aShares, bShares, dShares, groupSize);

  // build dShares k-1: c - d_0 - ... - d_{k-2}
  dShares[_K-1] = cShare;
  for (int i = 0; i< _K - 1; i++) {
    dShares[_K-1] = dShares[_K-1] - dShares[i];
  }

  buildPolyVecInd(aShares, bShares, dShares, groupSize);
  vector<FieldType> aSharesNew(groupSize);
  vector<FieldType> bSharesNew(groupSize);

  FieldType lambda = randomCoin();
  HIM<FieldType> matrix_for_k_lambda;
  HIM<FieldType> matrix_for_2k_lambda;
  vector<FieldType> beta_lambda(1);
  beta_lambda[0] = lambda;
  matrix_for_k_lambda.allocate(1, _K, _field);
  matrix_for_k_lambda.InitHIMByVectors(_alpha_k, beta_lambda);
  matrix_for_2k_lambda.allocate(1, 2*_K-1, _field);
  matrix_for_2k_lambda.InitHIMByVectors(_alpha_2k, beta_lambda);

  // TODO: refactor
  vector<FieldType> ySharesA;   // tmp vector
  vector<FieldType> ySharesB;   // tmp vector
  for (int i=0; i<groupSize; i++) {
    ySharesA.clear();
    ySharesB.clear();
    for (int j=i; j<totalLength; j+=groupSize) {
      ySharesA.push_back(aShares[j]);
      ySharesB.push_back(bShares[j]);
    }
    matrix_for_k_lambda.MatrixMult(ySharesA, y_for_interpolate);
    aSharesNew[i] = y_for_interpolate[0];
    matrix_for_k_lambda.MatrixMult(ySharesB, y_for_interpolate);
    bSharesNew[i] = y_for_interpolate[0];
  }

  matrix_for_2k_lambda.MatrixMult(dShares, y_for_interpolate);
  FieldType cShareNew = y_for_interpolate[0];
  compressVerifyVec(aSharesNew, bSharesNew, cShareNew);
  return;
}

template <class FieldType>
void ProtocolParty<FieldType>::
verifyVec(vector<FieldType>& aShares, vector<FieldType>& bShares,
          FieldType& cShare){

  // -- append random shares vec(a)_N, vec(b)_N, c_N
  int vecSize = aShares.size();
  vector<FieldType> aN(vecSize);
  getRandomShares(vecSize, aN);
  vector<FieldType> bN(vecSize);
  getRandomShares(vecSize, bN);

  aShares.insert(aShares.end(), aN.begin(), aN.end());
  bShares.insert(bShares.end(), bN.begin(), bN.end());
  vector<FieldType> dShares(2);
  DNMultVec(aShares, bShares, dShares, vecSize);

  cShare += dShares[1];
  dShares[1] = cShare - dShares[0];
  // -- build polynomials, and verify by open
  // interpolate A[i], B[i], and C
  vector<vector<FieldType>> AShares;
  vector<vector<FieldType>> BShares;
  vector<FieldType> CShare;
  buildPolyVec(aShares, bShares, cShare, dShares, vecSize,
               AShares, BShares, CShare);

  // eval at random coin
  vector<FieldType> ASecrets(vecSize);
  vector<FieldType> BSecrets(vecSize);
  FieldType lambda = randomCoin();
  vector<FieldType> CSecret(1, interp.evalPolynomial(lambda, CShare));
  evalPolyVecAt(AShares, ASecrets, lambda);
  evalPolyVecAt(BShares, BSecrets, lambda);

  // open and verify
  vector<FieldType> AResults(vecSize);
  openShare(vecSize, ASecrets, AResults);
  vector<FieldType> BResults(vecSize);
  openShare(vecSize, BSecrets, BResults);
  vector<FieldType> CResult(1);
  openShare(1, CSecret, CResult);

  for (int i = 0; i < vecSize; i++) {
    CResult[0] = CResult[0] - AResults[i] * BResults[i];
  }
  if (CResult[0] != _zero) {
    cout << "verification fail: " << CResult[0] << endl;
    abort();
  }
}


/**
 * the function implements the second step of Input Phase:
 * the party broadcasts for each input gate the difference between
 * the random secret and the actual input value.
 * @param diff
 */
template <class FieldType> void ProtocolParty<FieldType>::inputPhase() {

  readMyInputs();
  int fieldByteSize = _field->getElementSizeInBytes();
  vector<FieldType> x1(_N), y1(_N);
  vector<vector<FieldType>> sendBufsElements(_N);
  vector<vector<byte>> sendBufsBytes(_N);
  vector<vector<byte>> recBufBytes(_N);
  vector<vector<FieldType>> recBufElements(_N);

  // -- prepare the shares for the input
  int index = 0;
  vector<int> sizes(_N);
  for (int k = 0; k < numOfInputGates; k++) {
    if (circuit.getGates()[k].gateType == INPUT) {
      // get the expected sizes from the other parties
      sizes[circuit.getGates()[k].party]++;

      if (circuit.getGates()[k].party == _myId) {
        auto input = myInputs[index];
        index++;
        if (flag_print) {
          cout << "input  " << input << endl;
        }
        // the value of a_0 is the input of the party.
        x1[0] = _field->GetElement(input);

        // generate random degree-T polynomial
        for (int i = 1; i < _T + 1; i++) {
          // A random field element, uniform distribution
          x1[i] = _field->Random();
        }

        matrix_vand.MatrixMult(x1, y1, _T + 1);

        // prepare shares to be sent
        for (int i = 0; i < _N; i++) {
          sendBufsElements[i].push_back(y1[i]);
        }
      }
    } // else {abort();}
  }

  // -- convert shares to bytes
  for (int i = 0; i < _N; i++) {
    sendBufsBytes[i].resize(sendBufsElements[i].size() * fieldByteSize);
    recBufBytes[i].resize(sizes[i] * fieldByteSize);
    for (int j = 0; j < sendBufsElements[i].size(); j++) {
      _field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize),
                            sendBufsElements[i][j]);
    }
  }

  for (int i=0; i<_N; i++) {
    _comm.oneToAll(recBufBytes[i], sendBufsBytes, i, false);
  }

  // -- convert received bytes to shares
  for (int i = 0; i < _N; i++) {
    recBufElements[i].resize(sizes[i]);
    for (int j = 0; j < sizes[i]; j++) {
      recBufElements[i][j] =
        _field->bytesToElement(recBufBytes[i].data() + (j * fieldByteSize));
    }
  }

  // -- store shares in input gate order
  vector<int> counters(_N);
  for (int i = 0; i < _N; i++) {
    counters[i] = 0;
  }
  vector<FieldType> inputShares(this->numOfInputGates);
  for (int k = 0; k < this->numOfInputGates; k++) {
    if (circuit.getGates()[k].gateType == INPUT) {
      auto share = recBufElements[circuit.getGates()[k].party]
                                 [counters[circuit.getGates()[k].party]];
      counters[circuit.getGates()[k].party] += 1;
      gateShareArr[circuit.getGates()[k].output] = share;
      inputShares[k] = share;
    }
  }
}


template <class FieldType>
void ProtocolParty<FieldType>::getRandomShares(
    int numOfRandoms, vector<FieldType> &randomElementsToFill) {

  randomElementsToFill.assign(_singleSharesArray.begin() + _singleSharesOffset,
                              _singleSharesArray.begin() + _singleSharesOffset +
                              numOfRandoms);

  _singleSharesOffset += numOfRandoms;
}




template <class FieldType>
void ProtocolParty<FieldType>::initializationPhase() {
  // -- allocate spaces for protocol storage
  this->y_for_interpolate.resize(_N);
  this->gateShareArr.resize(M - this->numOfOutputGates);

  // -- initialize protocol parameters
  this->beta.resize(1);
  this->beta[0] = _field->GetElement(0); // zero of the field
  this->alpha.resize(_N);
  for (int i = 0; i < _N; i++) {
    this->alpha[i] = _field->GetElement(i + 1);
  }
  _alpha_2k.resize(_K*2 - 1);
  for (int i=0; i< _K*2 - 1; i++) {
    _alpha_2k[i] = _field->GetElement(i+1);
  }
  _alpha_k.resize(_K);
  for (int i=0; i< _K; i++) {
    _alpha_k[i] = _field->GetElement(i+1);
  }

  // Interpolate first T+1 positions (deg = T)
  matrix_for_t.allocate(_N - (1+_T), _T + 1, _field);
  matrix_for_t.InitHIMVectorAndsizes(alpha, _T + 1, _N - _T - 1);

  // Interpolate first 2T+1 positions (deg = 2T)
  matrix_for_2t.allocate(_N - (1+2*_T), 2*_T + 1, _field);
  matrix_for_2t.InitHIMVectorAndsizes(alpha, 2*_T + 1, _N - 1 - 2*_T);

  // Interpolate first K positions (deg = K-1)
  matrix_for_k.allocate(2*_K-1 - _K, _K, _field);
  matrix_for_k.InitHIMVectorAndsizes(_alpha_2k, _K, 2*_K-1 - _K);

  vector<FieldType> alpha_2n(_N*2);
  for (int i = 0; i < _N; i++) {
    alpha_2n[i] = _field->GetElement(i + 1);
    alpha_2n[i+_N] = _field->GetElement(i+_N+1);
  }
  _mat_n_to_n.allocate(_N, _N, _field);
  _mat_n_to_n.InitHIMVectorAndsizes(alpha_2n, _N, _N);

  matForRefresh(_myId, _mat_to_T);

  this->matrix_for_interpolate.allocate(1, _N, _field);
  this->matrix_for_interpolate.InitHIMByVectors(alpha, beta);
  this->matrix_vand.allocate(_N, _N, _field);
  this->matrix_vand.InitVDM(); // default to alpha as defined
  this->matrix_vand_transpose.allocate(_N, _N, _field);
  this->matrix_vand_transpose.InitVDMTranspose();
}

template <class FieldType> bool ProtocolParty<FieldType>::preparationPhase() {
  // ---- # of random single shares:
  // 2. Padding vector product verification in base case
  //    [2 * vecSize] // which is at most _K
  //    in total: 2 * _K
  // 3. Random Coins
  //    [1]
  //    in total: nCompressions + 2
  int nCompressions = (int)(log(this->numOfMultGates) / log(_K) + 0.5);
  int numSingleShares = 2 * _K + nCompressions + 2;
  _singleSharesOffset = 0;
  makeRandShares(numSingleShares, _singleSharesArray);
  // cout << "generated single Shares: " << numSingleShares << endl;
  // ---- # of random double shares:
  // 1. Compute Mult gates in the circuit
  //    -- in total: numOfMultGates
  // 2. Compress Verification
  //    -- 2 * _K * (nCompressions + 1)
  int numDoubleShares =
    this->numOfMultGates + (nCompressions+1)*_K*2;
  _doubleSharesOffset = 0;
  makeRandDoubleShares(numDoubleShares, _doubleSharesArray, _doubleKingIds);
  // cout << "generated doubles: " << numDoubleShares << endl;
  return true;
}

template <class FieldType> bool ProtocolParty<FieldType>::
preparationPhase(int nSingle, int nDouble) { // for special (external) use
  _singleSharesOffset = 0;
  makeRandShares(nSingle, _singleSharesArray);
  
  _doubleSharesOffset = 0;
  makeRandDoubleShares(nDouble, _doubleSharesArray, _doubleKingIds);
  // cout << "generated doubles: " << numDoubleShares << endl;
  return true;
}


template <class FieldType>
bool ProtocolParty<FieldType>::checkConsistency(vector<FieldType> &x, int d) {
  if (d == _T) {
    vector<FieldType> y(_N - 1 - d); // the result of multiplication
    vector<FieldType> x_until_t(_T + 1);

    for (int i = 0; i < _T + 1; i++) {
      x_until_t[i] = x[i];
    }

    matrix_for_t.MatrixMult(x_until_t, y);

    // compare that the result is equal to the according positions in x
    for (int i = 0; i < _N - d - 1; i++)
    {
      if ((y[i]) != (x[d + 1 + i])) {
        return false;
      }
    }
    return true;
  } else if (d == 2 * _T) {
    vector<FieldType> y(_N - 1 - d); // the result of multiplication

    vector<FieldType> x_until_2t(2 * _T + 1);

    for (int i = 0; i < 2 * _T + 1; i++) {
      x_until_2t[i] = x[i];
    }

    matrix_for_2t.MatrixMult(x_until_2t, y);

    // compare that the result is equal to the according positions in x
    for (int i = 0; i < _N - d - 1; i++)
    {
      if ((y[i]) != (x[d + 1 + i])) {
        return false;
      }
    }
    return true;

  } else {
    vector<FieldType> alpha_until_d(d + 1);
    vector<FieldType> alpha_from_d(_N - 1 - d);
    vector<FieldType> x_until_d(d + 1);
    vector<FieldType> y(_N - 1 - d); // the result of multiplication

    for (int i = 0; i < d + 1; i++) {
      alpha_until_d[i] = alpha[i];
      x_until_d[i] = x[i];
    }
    for (int i = d + 1; i < _N; i++) {
      alpha_from_d[i - (d + 1)] = alpha[i];
    }
    // Interpolate first d+1 positions of (alpha,x)
    HIM<FieldType> matrix(_N - 1 - d, d + 1,
                          _field); // slices, only positions from 0..d
    matrix.InitHIMByVectors(alpha_until_d, alpha_from_d);
    matrix.MatrixMult(x_until_d, y);

    // compare that the result is equal to the according positions in x
    for (int i = 0; i < _N - d - 1; i++)
    {
      if (y[i] != x[d + 1 + i]) {
        return false;
      }
    }
    return true;
  }
  return true;
}


template <class FieldType>
FieldType ProtocolParty<FieldType>::reconstructShare(vector<FieldType> &x,
                                                     int d) {

  if (!checkConsistency(x, d)) {
    cout << "reconstructing " << d << "-share fail." << '\n';
    abort();
  } else
    return interpolate_to_zero(x);
}

// Interpolate polynomial at position Zero
template <class FieldType>
FieldType ProtocolParty<FieldType>::interpolate_to_zero(vector<FieldType> &x) {
  // vector<FieldType> y(N); // result of interpolate
  matrix_for_interpolate.MatrixMult(x, y_for_interpolate);
  return y_for_interpolate[0];
}


template <class FieldType> int ProtocolParty<FieldType>::processNotMult() {
  int count = 0;
  for (int k = circuit.getLayers()[currentCircuitLayer];
       k < circuit.getLayers()[currentCircuitLayer + 1]; k++) {

    // add gate
    if (circuit.getGates()[k].gateType == ADD) {
      gateShareArr[circuit.getGates()[k].output] =
          gateShareArr[circuit.getGates()[k].input1] +
          gateShareArr[circuit.getGates()[k].input2];
      count++;
    }

    else if (circuit.getGates()[k].gateType == SUB) // sub gate
    {
      gateShareArr[circuit.getGates()[k].output] =
          gateShareArr[circuit.getGates()[k].input1] -
          gateShareArr[circuit.getGates()[k].input2];
      count++;
    } else if (circuit.getGates()[k].gateType == SCALAR) {
      long scalar(circuit.getGates()[k].input2);
      FieldType e = _field->GetElement(scalar);
      gateShareArr[circuit.getGates()[k].output] =
          gateShareArr[circuit.getGates()[k].input1] * e;
      count++;
    } else if (circuit.getGates()[k].gateType == SCALAR_ADD) {
      long scalar(circuit.getGates()[k].input2);
      FieldType e = _field->GetElement(scalar);
      gateShareArr[circuit.getGates()[k].output] =
          gateShareArr[circuit.getGates()[k].input1] + e;
      count++;
    } // else ignore mult gates
  }
  return count;
}

/**
 * the Function process all multiplications which are ready.
 * @return the number of processed gates.
 */
template <class FieldType>
int ProtocolParty<FieldType>::processMultiplications(int lastMultGate) {

  return processMultDN(lastMultGate);
}

template <class FieldType>
int ProtocolParty<FieldType>::processMultDN(int indexInRandomArray) {


  int fieldByteSize = _field->getElementSizeInBytes();
  int maxNumberOfLayerMult = circuit.getLayers()[currentCircuitLayer + 1] -
                             circuit.getLayers()[currentCircuitLayer];
  vector<FieldType> xyMinusRShares(maxNumberOfLayerMult);
  vector<FieldType> xyMinusR;
  vector<byte> xyMinusRBytes;
  vector<vector<byte>> recBufsBytes(_N);
  vector<vector<byte>> sendBufsBytes(_N);
  vector<vector<FieldType>> sendBufsElements(_N);

  // -- go through current layer to collect mult gates
  int index = 0;
  for (int k = circuit.getLayers()[currentCircuitLayer];
       k < circuit.getLayers()[currentCircuitLayer + 1]; k++) {
    auto gate = circuit.getGates()[k];

    if (gate.gateType == MULT) {
      // compute the share of xy-r
      xyMinusRShares[index] =
          gateShareArr[gate.input1] * gateShareArr[gate.input2] -
          _doubleSharesArray[2 * indexInRandomArray + 1];
      indexInRandomArray++;
      index++;
    }
  }

  // -- early stop if no mult in current layer
  if (index == 0)
    return 0;

  xyMinusRShares.resize(index);
  batchMultSingle(xyMinusRShares, xyMinusR);
  // -- store mult results to gate wire array
  indexInRandomArray -= index;
  index = 0;
  for (int k = circuit.getLayers()[currentCircuitLayer];
       k < circuit.getLayers()[currentCircuitLayer + 1]; k++) {
    auto gate = circuit.getGates()[k];

    if (gate.gateType == MULT) {
      gateShareArr[gate.output] =
          _doubleSharesArray[2 * indexInRandomArray] + xyMinusR[index];
      indexInRandomArray++;
      index++;
    }
  }

  return index;
}

template <class FieldType>
void ProtocolParty<FieldType>::
interpolatePolyVec(vector<FieldType>& aShares,
                   vector<vector<FieldType>>& AShares, int groupSize) {
  int totalLength = aShares.size();
  int polyDegree = (totalLength + groupSize -1) / groupSize;
  vector<FieldType> yShares;
  for (int i = 0; i < groupSize; i++) {
    yShares.clear();
    for (int j = i; j < totalLength; j += groupSize) {
      yShares.push_back(aShares[j]);
    }
    interp.interpolate(_alpha_2k, yShares, AShares[i]);
    AShares[i].resize(polyDegree, _field->GetElement(0));
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
evalPolyVecAt(vector<vector<FieldType>>& AShares,
              vector<FieldType>& aShares, FieldType point) {
  int groupSize = AShares.size();
  aShares.resize(groupSize);
  for (int i = 0; i < groupSize; i++) {
    aShares[i] = interp.evalPolynomial(point, AShares[i]);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
evalPolyVec(vector<vector<FieldType>>& AShares,
            vector<FieldType>& aShares, int nPoints, int offset) {
  int groupSize = AShares.size();
  aShares.resize(nPoints * groupSize);
  for (int i = 0; i < nPoints; i++) {
    for (int j = 0; j < groupSize; j++) {
      aShares[i*groupSize + j] =
        interp.evalPolynomial(_alpha_2k[i+offset], AShares[j]);
    }
  }
}


template <class FieldType>
__inline__
void ProtocolParty<FieldType>::
batchMultSingle(vector<FieldType>& e2Shares, vector<FieldType>& e1Shares) {
  // compatible names from new_protocol_god
  int _doubleOffset = _doubleSharesOffset;

  auto t1 = high_resolution_clock::now();
  int nMults = e2Shares.size();
  vector<vector<FieldType>> e2SharesLoad(_N);
  vector<int> loadSizes(_N, 0);
  for (int i=0; i<nMults; i++) {
    int king = _doubleKingIds[i + _doubleOffset];
    e2SharesLoad[king].push_back(e2Shares[i]);
    loadSizes[king]++;
  }
  int maxLoad = loadSizes[0];
  for(int i=0; i<_N; i++) {
    maxLoad = max(maxLoad, loadSizes[i]);
  }
  vector<vector<byte>> recBufs(_N, vector<byte>(maxLoad*_fieldByteSize)),
    e2SharesByte(_N, vector<byte>(maxLoad*_fieldByteSize));
  for (int i=0; i<_N; i++) {
      _field->elementVectorToByteVector(e2SharesLoad[i], e2SharesByte[i]);
  }
  _comm.allToAll(e2SharesByte, recBufs);

  // -- all non-corr parties may be king
  int myLoad = e2SharesLoad[_myId].size();
  vector<FieldType> secrets(myLoad), e2ShareAll(_N);
  for (int i=0; i<myLoad; i++) {
    for (int j = 0; j < _N; j++) {
      e2ShareAll[j] =
        _field->bytesToElement(recBufs[j].data() + (i * _fieldByteSize));
    }
    secrets[i] = interpolate_to_zero(e2ShareAll);
  }

  vector<vector<FieldType>> eSharesAll(_N, vector<FieldType>(myLoad, _zero));
  makeTShares(secrets, eSharesAll, true);

  vector<vector<byte>>
    e1SharesByte(_N, vector<byte>(maxLoad * _fieldByteSize, 0)),
    sendBufs(_N, vector<byte>(maxLoad * _fieldByteSize, 0));
  for (int i=0; i < _N; i++) {
    if (_TMasks[_myId][i]) {
      _field->elementVectorToByteVector(eSharesAll[i], sendBufs[i]);
    }
  }
  _comm.allToT (sendBufs, e1SharesByte);
  // -- convert byte to elements
  vector<int> idx(_N, 0);
  e1Shares.clear();
  e1Shares.resize(nMults, _zero);
  for (int i = 0; i < nMults; i++) {
    int king = _doubleKingIds[i + _doubleOffset];
    if (_TMasks[king][_myId]) {
      e1Shares[i] = _field->bytesToElement(e1SharesByte[king].data() +
                                           (idx[king]++) * _fieldByteSize);
    }
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::
DNMultVec(vector<FieldType>& a, vector<FieldType>& b,
          vector<FieldType>& c, int groupSize) {
  // assert( a.size() == b.size() );
  int totalLength = a.size();
  int numOfMults = (totalLength + groupSize - 1) / groupSize;
  int fieldByteSize = _field->getElementSizeInBytes();
  c.resize(numOfMults);
  vector<FieldType> xyMinusR2T(numOfMults);
  vector<byte> xyMinusRSharesBytes(numOfMults * fieldByteSize);
  vector<FieldType> xyMinusRT(numOfMults);
  vector<byte> xyMinusRBytes(numOfMults * fieldByteSize);
  vector<vector<byte>> recBufsBytes;

  // -- generate the 2t-sharings for xy - r
  for (int group_i = 0; group_i < numOfMults; group_i++) {
    int group_end = (group_i + 1) * groupSize > totalLength ?
      totalLength : (group_i + 1) * groupSize;
    for (int i = group_i * groupSize; i < group_end; i++) {
      xyMinusR2T[group_i] += a[i] * b[i];
    }
    xyMinusR2T[group_i] = xyMinusR2T[group_i] -
      _doubleSharesArray[(_doubleSharesOffset + group_i)*2 + 1];
  }
  batchMultSingle(xyMinusR2T, xyMinusRT);

  // -- compute xy - r + [r]_t = t-sharing of xy
  for (int group_i = 0; group_i < numOfMults; group_i++) {
    c[group_i] =
      _doubleSharesArray[(_doubleSharesOffset + group_i)*2] + xyMinusRT[group_i];
  }

  _doubleSharesOffset += numOfMults;
}

template <class FieldType>
void ProtocolParty<FieldType>::openShare(int numOfRandomShares,
                                         vector<FieldType> &Shares,
                                         vector<FieldType> &secrets) {
  secrets.resize(numOfRandomShares);
  vector<vector<byte>> sendBufsBytes(_N);
  vector<vector<byte>> recBufsBytes(_N);

  vector<FieldType> x1(_N);
  int fieldByteSize = _field->getElementSizeInBytes();

  // resize vectors
  for (int i = 0; i < _N; i++) {
    sendBufsBytes[i].resize(numOfRandomShares * fieldByteSize);
    recBufsBytes[i].resize(numOfRandomShares * fieldByteSize);
  }

  // set the first sending data buffer
  for (int j = 0; j < numOfRandomShares; j++) {
    _field->elementToBytes(sendBufsBytes[0].data() + (j * fieldByteSize),
                          Shares[j]);
  }

  // copy the same data for all parties
  for (int i = 1; i < _N; i++) {
    sendBufsBytes[i] = sendBufsBytes[0];
  }

  _comm.allToAll(sendBufsBytes, recBufsBytes);

  // reconstruct each set of shares to get the secret

  for (int k = 0; k < numOfRandomShares; k++) {
    // get the set of shares for each element
    for (int i = 0; i < _N; i++) {
      x1[i] = _field->bytesToElement(recBufsBytes[i].data() +
                                    (k * fieldByteSize));
    }
    secrets[k] = reconstructShare(x1, _T);
  }
}

template <class FieldType>
FieldType ProtocolParty<FieldType>::randomCoin() {
  vector<FieldType> lambdaShare(1);
  vector<FieldType> lambda(1);
  getRandomShares(1, lambdaShare);
  openShare(1, lambdaShare, lambda);

  vector<vector<byte>> sendEmpty(_N), recEmpty(_N);
  _comm.allToAll(sendEmpty, recEmpty); // empty round to sync
  return lambda[0];
}

/**
 * the function Walk through the circuit and reconstruct output gates.
 * @param circuit
 * @param gateShareArr
 * @param alpha
 */
template <class FieldType> void ProtocolParty<FieldType>::outputPhase() {
  vector<FieldType> x1(_N); // vector for the shares of my outputs
  vector<vector<FieldType>> sendBufsElements(_N);
  vector<vector<byte>> sendBufsBytes(_N);
  vector<vector<byte>> recBufBytes(_N);

  FieldType num;
  ofstream myfile;
  myfile.open(this->outputFile);

  for (int k = M - numOfOutputGates; k < M; k++) {
    if (circuit.getGates()[k].gateType == OUTPUT) {
      // send to party (which need this gate) your share for this gate
      sendBufsElements[circuit.getGates()[k].party].push_back(
          gateShareArr[circuit.getGates()[k].input1]);
    }
  }

  int fieldByteSize = _field->getElementSizeInBytes();
  for (int i = 0; i < _N; i++) {
    sendBufsBytes[i].resize(sendBufsElements[i].size() * fieldByteSize);
    recBufBytes[i].resize(sendBufsElements[_myId].size() * fieldByteSize);
    for(int j=0; j<sendBufsElements[i].size();j++) {
      _field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize),
                            sendBufsElements[i][j]);
    }
  }

  for (int i=0; i<_N; i++) {
    _comm.allToOne(sendBufsBytes[i], recBufBytes, i);
  }

  int counter = 0;
  for (int k = M - numOfOutputGates; k < M; k++) {
    if (circuit.getGates()[k].gateType == OUTPUT &&
        circuit.getGates()[k].party == _myId) {
      for (int i = 0; i < _N; i++) {
        x1[i] = _field->bytesToElement(recBufBytes[i].data() +
                                      (counter * fieldByteSize));
      }

      if (flag_print_output)
        cout << "the result for " << circuit.getGates()[k].input1
             << " is : " << _field->elementToString(reconstructShare(x1, _T+1)) << '\n';
      counter++;
    }
  }

  // close output file
  myfile.close();
}

template <class FieldType> ProtocolParty<FieldType>::~ProtocolParty() {
  protocolTimer->writeToFile();
  delete protocolTimer;
  delete _field;
  delete timer;
  // delete comm;
}
#endif /* PROTOCOLPARTY_H_ */
