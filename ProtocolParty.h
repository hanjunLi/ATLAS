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
#include <libscapi/include/primitives/Prg.hpp>
#include <thread>
#include <vector>

#define flag_print false
#define flag_print_timings true
#define flag_print_output true

using namespace std;
using namespace std::chrono;

template <class FieldType>
class ProtocolParty : public Protocol, public HonestMajority, MultiParty {

private:
  // -- global const
  int numThreads = 1;           // TODO: add as main arguments later
  
  // -- global variables
  int iteration;                       // current iteration number
  int currentCircuitLayer = 0;         // current circuit layer
  int offset = 0;                      // next random double share to use
  int randomSharesOffset = 0;          // next random share to use
  vector<FieldType> y_for_interpolate; // size = N
  boost::asio::io_service io_service;

  // -- set in constructor
  int times;                // number of times to run the run function
  Measurement *timer;
  TemplateField<FieldType> *field; // used to call some field functions
  ProtocolTimer *protocolTimer;
  // N, T - number of parties, malicious
  int N, T, m_partyId;  
  string s;                     // string of my_partyId
  string inputsFile, outputFile;
  ArithmeticCircuit circuit;
  // M - number of gates
  int M, numOfInputGates, numOfOutputGates, numOfMultGates;
  vector<long> myInputs;
  vector<shared_ptr<ProtocolPartyData>> parties;

  // -- set in initialization phase
  vector<FieldType> bigR;         // a single random element
  vector<FieldType> beta;         // a single zero element
  vector<FieldType> alpha;        // N distinct non-zero field elements
  HIM<FieldType> matrix_for_interpolate;
  HIM<FieldType> matrix_for_t;
  HIM<FieldType> matrix_for_2t;  
  VDM<FieldType> matrix_vand;
  VDMTranspose<FieldType> matrix_vand_transpose;  

  // -- filled during offline-preparation phase
  vector<FieldType> randomSharesArray;  // used as seed for AES PRG
  vector<FieldType> randomTAnd2TShares; // used for DN mult
  vector<FieldType> gateShareArr; // my share of each gate (1 share each)
  vector<FieldType> gateValueArr; // the value for my input and output
  
public:
  bool hasOffline() { return true; }
  bool hasOnline() override { return true; }
  ~ProtocolParty();
  
  // -- protocol functions
  ProtocolParty(int argc, char *argv[]);
  void readMyInputs();
  void initializationPhase();
  
  void run() override;
  void runOffline() override;
  bool preparationPhase();
  void offlineDNForMultiplication(int numOfDoubleShares);
  
  void runOnline() override;
  void inputPhase();

  void computationPhase();
  int processNotMult();
  int processMultiplications(int lastMultGate);
  int processMultDN(int indexInRandomArray);
  
  void verificationPhase();
  bool comparingViews();        // TODO: most likely not needed
  void outputPhase();
  
  // -- round functions (followed by thread functions)
  void roundFunctionSync(vector<vector<byte>> &sendBufs,
                         vector<vector<byte>> &recBufs, int round);
  void exchangeData(vector<vector<byte>> &sendBufs,
                    vector<vector<byte>> &recBufs, int first, int last);
  void recToP1(vector<byte> &myShare, vector<vector<byte>> &recBufs);
  void recDataToP1(vector<vector<byte>> &recBufs, int first, int last);
  void sendFromP1(vector<byte> &sendBuf);
  void sendDataFromP1(vector<byte> &sendBuf, int first, int last);

  // -- secret sharing functionalities
  void generateRandomShares(int numOfRandoms,
                            vector<FieldType> &randomElementsToFill);
  void generateRandom2TAndTShares(int numOfRandomPairs,
                                  vector<FieldType> &randomElementsToFill);
  void batchConsistencyCheckOfShares(const vector<FieldType> &inputShares);
  void openShare(int numOfRandomShares, vector<FieldType> &Shares,
                 vector<FieldType> &secrets);

  // --  helper functions
  void getRandomShares(int numOfRandoms,
                       vector<FieldType> &randomElementsToFill);
  void getRandomSharesWithCheck(int numOfRnadoms,
                                vector<FieldType> &randomElementsToFill);
  vector<byte> generateCommonKey();
  void generatePseudoRandomElements(vector<byte> &aesKey,
                                    vector<FieldType> &randomElementsToFill,
                                    int numOfRandomElements);
  bool checkConsistency(vector<FieldType> &x, int d);
  // locally interpolate received sharing
  FieldType reconstructShare(vector<FieldType> &x, int d);
  FieldType interpolate(vector<FieldType> &x);
};

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
    this->field = new TemplateField<FieldType>(2147483647);
  } else if (fieldType.compare("ZpMersenne61") == 0) {
    this->field = new TemplateField<FieldType>(0);
  } else if (fieldType.compare("ZpKaratsuba") == 0) {
    this->field = new TemplateField<FieldType>(0);
  } else if (fieldType.compare("GF2E") == 0) {
    this->field = new TemplateField<FieldType>(8);
  } else if (fieldType.compare("Zp") == 0) {
    this->field = new TemplateField<FieldType>(2147483647);
  }
  string circuitFile =
    this->getParser().getValueByKey(arguments, "circuitFile");
  string outputTimerFileName =
    circuitFile + "Times" + to_string(m_partyId) + fieldType + ".csv";
  this->protocolTimer = new ProtocolTimer(times, outputTimerFileName);
  this->N = stoi(this->getParser().getValueByKey(arguments, "partiesNumber"));
  this->T = (this->N - 1)/ 2;
  this->m_partyId = stoi(this->getParser().getValueByKey(arguments, "partyID"));
  this->s = to_string(m_partyId);
  this->inputsFile = this->getParser().getValueByKey(arguments, "inputFile");
  this->outputFile = this->getParser().getValueByKey(arguments, "outputFile");
  this->circuit.readCircuit(circuitFile.c_str());
  this->circuit.reArrangeCircuit();
  this->M = circuit.getNrOfGates();
  this->numOfInputGates = circuit.getNrOfInputGates();
  this->numOfOutputGates = circuit.getNrOfOutputGates();
  this->numOfMultGates = circuit.getNrOfMultiplicationGates();
  this->myInputs.resize(numOfInputGates);
  readMyInputs();
  MPCCommunication comm;
  string partiesFile =
      this->getParser().getValueByKey(arguments, "partiesFile");
  this->parties =
    comm.setCommunication(this->io_service, m_partyId, N, partiesFile);
  
  // -- test communication
  string tmp = "init times";
  byte tmpBytes[20];
  for (int i = 0; i < parties.size(); i++) {
    if (parties[i]->getID() < m_partyId) {
      parties[i]->getChannel()->write(tmp);
      parties[i]->getChannel()->read(tmpBytes, tmp.size());
    } else {
      parties[i]->getChannel()->read(tmpBytes, tmp.size());
      parties[i]->getChannel()->write(tmp);
    }
  }

  // -- initialize phase
  auto t1 = high_resolution_clock::now();
  initializationPhase();
  auto t2 = high_resolution_clock::now();

  auto duration = duration_cast<milliseconds>(t2 - t1).count();
  if (flag_print_timings) {
    cout << "time in milliseconds initializationPhase: " << duration << endl;
  }
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
    cout << "time in milliseconds for protocol: " << duration << endl;
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
  if (flag_print_timings) {
    cout << "time in milliseconds preparationPhase: " << duration << endl;
  }
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
  if (flag_print_timings) {
    cout << "time in milliseconds inputPhase: " << duration << endl;
  }

  t1 = high_resolution_clock::now();
  timer->startSubTask("ComputePhase", iteration);
  computationPhase();
  timer->endSubTask("ComputePhase", iteration);
  t2 = high_resolution_clock::now();

  duration = duration_cast<milliseconds>(t2 - t1).count();
  protocolTimer->computationPhaseArr[iteration] = duration;

  if (flag_print_timings) {
    cout << "time in milliseconds computationPhase: " << duration << endl;
  }

  t1 = high_resolution_clock::now();
  timer->startSubTask("VerificationPhase", iteration);
  verificationPhase();
  timer->endSubTask("VerificationPhase", iteration);
  t2 = high_resolution_clock::now();
  duration = duration_cast<milliseconds>(t2 - t1).count();
  protocolTimer->verificationPhaseArr[iteration] = duration;

  if (flag_print_timings) {
    cout << "time in milliseconds verificationPhase: " << duration << endl;
  }

  t1 = high_resolution_clock::now();
  timer->startSubTask("outputPhase", iteration);
  outputPhase();
  timer->endSubTask("outputPhase", iteration);
  t2 = high_resolution_clock::now();

  duration = duration_cast<milliseconds>(t2 - t1).count();
  protocolTimer->outputPhaseArr[iteration] = duration;

  if (flag_print_timings) {
    cout << "time in milliseconds outputPhase: " << duration << endl;
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::computationPhase() {
  int countNumMult = this->offset / 2;
  int numOfLayers = circuit.getLayers().size();
  for (int i = 0; i < numOfLayers - 1; i++) {
    currentCircuitLayer = i;
    // send the index of the current mult gate
    processNotMult();
    countNumMult += processMultiplications(countNumMult);
  }
  this->offset = countNumMult;
}

template <class FieldType>
void ProtocolParty<FieldType>::verificationPhase() {
  vector<FieldType> aShares(this->numOfMultGates);
  vector<FieldType> bShares(this->numOfMultGates);
  vector<FieldType> cShares(this->numOfMultGates);

  // -- gather all mult triples to be verified
  int idx = 0;
  for (auto& gate : this->circuit.getGate()) {
    if(gate.gateType == MULT) {
      aShares[idx] = gateShareArr[gate.input1];
      bShares[idx] = gateShareArr[gate.input2];
      cShares[idx] = gateShareArr[gate.output];
      idx++;
    }
  }

  // -- append random shares a_N, b_N
  vector<FieldType> twoRandShares(2);
  getRandomSharesWithCheck(2, twoRandShares);
  aShares.push_back(twoRandShares[0]);
  bShares.push_back(twoRandShares[1]);

  // -- prepare HIM matrices for evaluating A(N+1, ..., 2N-1)
  // -- it's  N-1 row, N cols (N = numMultGates + 1)
  vector<FieldType> alpha_tmp(this->numMultGates+1);
  vector<FieldType> beta_tmp(this->numMultGates);
  // TODO Start from here: what batch size to use? FFT for interpolation?
}


/**
 * the function implements the second step of Input Phase:
 * the party broadcasts for each input gate the difference between
 * the random secret and the actual input value.
 * @param diff
 */
template <class FieldType> void ProtocolParty<FieldType>::inputPhase() {
  
  int fieldByteSize = field->getElementSizeInBytes();
  vector<FieldType> x1(N), y1(N);
  vector<vector<FieldType>> sendBufsElements(N);
  vector<vector<byte>> sendBufsBytes(N);
  vector<vector<byte>> recBufBytes(N);
  vector<vector<FieldType>> recBufElements(N);
  
  // -- prepare the shares for the input
  int index = 0;
  vector<int> sizes(N);
  for (int k = 0; k < numOfInputGates; k++) {
    if (circuit.getGates()[k].gateType == INPUT) {
      // get the expected sizes from the other parties
      sizes[circuit.getGates()[k].party]++;

      if (circuit.getGates()[k].party == m_partyId) {
        auto input = myInputs[index];
        index++;
        if (flag_print) {
          cout << "input  " << input << endl;
        }
        // the value of a_0 is the input of the party.
        x1[0] = field->GetElement(input);

        // generate random degree-T polynomial
        for (int i = 1; i < T + 1; i++) {
          // A random field element, uniform distribution
          x1[i] = field->Random();
        }

        matrix_vand.MatrixMult(x1, y1, T + 1);

        // prepare shares to be sent
        for (int i = 0; i < N; i++) {
          sendBufsElements[i].push_back(y1[i]);
        }
      }
    } // else {abort();}
  }

  // -- convert shares to bytes
  for (int i = 0; i < N; i++) {
    sendBufsBytes[i].resize(sendBufsElements[i].size() * fieldByteSize);
    recBufBytes[i].resize(sizes[i] * fieldByteSize);
    for (int j = 0; j < sendBufsElements[i].size(); j++) {
      field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize),
                            sendBufsElements[i][j]);
    }
  }

  roundFunctionSync(sendBufsBytes, recBufBytes, 7);

  // -- convert received bytes to shares
  for (int i = 0; i < N; i++) {
    recBufElements[i].resize(sizes[i]);    
    for (int j = 0; j < sizes[i]; j++) {
      recBufElements[i][j] =
        field->bytesToElement(recBufBytes[i].data() + (j * fieldByteSize));
    }
  }

  // -- store shares in input gate order
  vector<int> counters(N);
  for (int i = 0; i < N; i++) {
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
  inputVerification(inputShares);
}

template <class FieldType>
void ProtocolParty<FieldType>::generateRandomShares(
    int numOfRandoms, vector<FieldType> &randomElementsToFill) {
  // -- allocate communication buffers
  vector<FieldType> x1(N), y1(N), t1(N), r1(N);
  // every round N-T random shares are generated.
  int fieldByteSize = field->getElementSizeInBytes();
  int no_buckets = (numOfRandoms / (N - T)) + 1;
  randomElementsToFill.resize(no_buckets * (N - T));
  // N vectors of no_buckets size
  vector<vector<FieldType>> sendBufsElements(N);
  // N vectors of no_buckets elements in bytes
  vector<vector<byte>> sendBufsBytes(N);
  vector<vector<byte>> recBufsBytes(N);
  for (int i = 0; i < N; i++) {
    sendBufsElements[i].resize(no_buckets);
    sendBufsBytes[i].resize(no_buckets * fieldByteSize);
    recBufsBytes[i].resize(no_buckets * fieldByteSize);
  }

  // -- generate no_buckets number of random t-sharings
  for (int k = 0; k < no_buckets; k++) {
    // generate random degree-T polynomial (T+1 ... N positions are zero)
    for (int i = 0; i < T + 1; i++) {
      x1[i] = field->Random();
    }
    // evaluate using vandermonde matrix (over alpha)
    matrix_vand.MatrixMult(x1, y1, T + 1);

    // prepare shares to be sent
    for (int i = 0; i < N; i++) {
      sendBufsElements[i][k] = y1[i];
    }
  }
  // -- convert random elements to bytes
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < no_buckets; j++) {
      field->elementToBytes(sendBufsBytes[i].data() + j*fieldByteSize,
                            sendBufsElements[i][j]);
    }
  }
  // -- exchange elements all-to-all
  roundFunctionSync(sendBufsBytes, recBufsBytes, 4); // round 4
  // -- convert bytes to random elements and extract to result
  int index = 0;
  for (int k = 0; k < no_buckets; k++) {
    for (int i = 0; i < N; i++) {
      t1[i] = field->bytesToElement(recBufsBytes[i].data() +
                                    (k * fieldByteSize));
    }
    matrix_vand_transpose.MatrixMult(t1, r1, N - T);
    // copy the resulting vector to the array of randoms
    for (int i = 0; i < N - T; i++) {
      randomElementsToFill[index++] = r1[i];
    }
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::getRandomShares(
    int numOfRandoms, vector<FieldType> &randomElementsToFill) {

  randomElementsToFill.assign(randomSharesArray.begin() + randomSharesOffset,
                              randomSharesArray.begin() + randomSharesOffset +
                              numOfRandoms);

  this->randomSharesOffset += numOfRandoms;
}

template <class FieldType>
void ProtocolParty<FieldType>::getRandomSharesWithCheck(
    int numOfRandoms, vector<FieldType> &randomElementsToFill) {

  getRandomShares(numOfRandoms, randomElementsToFill);

  batchConsistencyCheckOfShares(randomElementsToFill);
}

template <class FieldType>
void ProtocolParty<FieldType>::generateRandom2TAndTShares(
    int numOfRandomPairs, vector<FieldType> &randomElementsToFill) {
  // -- allocate communication buffers
  vector<FieldType> x1(N), y1(N), x2(N), y2(N), t1(N), r1(N), t2(N), r2(N);
  int fieldByteSize = field->getElementSizeInBytes();
  int no_buckets = (numOfRandomPairs / (N - T)) + 1;
  randomElementsToFill.resize(no_buckets * (N - T) * 2);
  vector<vector<FieldType>> sendBufsElements(N);
  vector<vector<byte>> sendBufsBytes(N);
  vector<vector<byte>> recBufsBytes(N);
  // a copy of the t-shares
  vector<FieldType> randomElementsOnlyTshares(no_buckets * (N - T));
  for (int i = 0; i < N; i++) {
    sendBufsElements[i].resize(no_buckets * 2);
    sendBufsBytes[i].resize(no_buckets * 2 * fieldByteSize);
    recBufsBytes[i].resize(no_buckets * 2 * fieldByteSize);
  }

  // -- generate no_buckets number of random t-sharings, and 2t-sharings
  for (int k = 0; k < no_buckets; k++) {
    // generate random degree-T polynomial
    for (int i = 0; i < T + 1; i++) {
      x1[i] = field->Random();
    }
    matrix_vand.MatrixMult(x1, y1, T + 1);
    x2[0] = x1[0];
    // generate random degree-2T polynomial
    for (int i = 1; i < 2 * T + 1; i++) {
      x2[i] = field->Random();
    }
    matrix_vand.MatrixMult(x2, y2, 2 * T + 1);

    // prepare shares to be sent
    for (int i = 0; i < N; i++) {
      // cout << "y1[ " <<i<< "]" <<y1[i] << endl;
      sendBufsElements[i][2 * k] = y1[i];
      sendBufsElements[i][2 * k + 1] = y2[i];
    }
  }

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < no_buckets * 2; j++) {
      field->elementToBytes(sendBufsBytes[i].data() + j*fieldByteSize,
                            sendBufsElements[i][j]);
    }    
  }
  // -- exchange elements all-to-all
  roundFunctionSync(sendBufsBytes, recBufsBytes, 5); // round 5
  // -- convert bytes to random elements and extract to result
  int index = 0;
  for (int k = 0; k < no_buckets; k++) {
    for (int i = 0; i < N; i++) {
      t1[i] = field->bytesToElement(recBufsBytes[i].data() +
                                    (2 * k * fieldByteSize));
      t2[i] = field->bytesToElement(recBufsBytes[i].data() +
                                    ((2 * k + 1) * fieldByteSize));
    }
    matrix_vand_transpose.MatrixMult(t1, r1, N - T);
    matrix_vand_transpose.MatrixMult(t2, r2, N - T);
    // copy the resulting vector to the array of randoms
    for (int i = 0; i < (N - T); i++) {
      randomElementsToFill[index * 2] = r1[i];
      // also copy to a separate vector for checking consistency
      randomElementsOnlyTshares[index] = r1[i];
      randomElementsToFill[index * 2 + 1] = r2[i];
      index++;
    }
  }

  // -- check validity of the t-shares
  // 2t-shares do not have to be checked
  batchConsistencyCheckOfShares(randomElementsOnlyTshares);  
}

template <class FieldType>
void ProtocolParty<FieldType>::batchConsistencyCheckOfShares(
    const vector<FieldType> &inputShares) { 
  // first generate the common aes key
  auto key = generateCommonKey();
  // calc the number of times we need to run the verification -- ceiling
  int iterations =
      (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();

  vector<FieldType> randomElements(inputShares.size() * iterations);
  generatePseudoRandomElements(key, randomElements, inputShares.size());

  for (int j = 0; j < iterations; j++) {
    vector<FieldType> r(1); // vector holding the random shares generated
    vector<FieldType> v(1);
    vector<FieldType> secret(1);

    getRandomShares(1, r);

    for (int i = 0; i < inputShares.size(); i++)
      v[0] += randomElements[i + j * inputShares.size()] * inputShares[i];

    v[0] += r[0];

    // if all the the parties share lie on the same polynomial this will not
    // abort
    openShare(1, v, secret);
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::initializationPhase() {
  // -- allocate spaces for protocol storage
  this->y_for_interpolate.resize(N);
  this->gateShareArr.resize(M - this->numOfOutputGates);

  // -- initialize protocol parameters
  this->bigR.resize(1);
  this->beta.resize(1);
  this->beta[0] = field->GetElement(0); // zero of the field
  this->alpha.resize(N);
  for (int i = 0; i < N; i++) {
    this->alpha[i] = field->GetElement(i + 1);
  }

  // Interpolate first d+1 positions of (alpha,x)
  matrix_for_t.allocate(N - 1 - T, T + 1, field); 
  // slices, only positions from 0..d
  matrix_for_t.InitHIMVectorAndsizes(alpha, T + 1, N - T - 1);

  // Interpolate first d+1 positions of (alpha,x)
  matrix_for_2t.allocate(N - 1 - 2*T, 2*T + 1, field);
  // slices, only positions from 0..d
  matrix_for_2t.InitHIMVectorAndsizes(alpha, 2*T + 1, N - 1 - 2*T);
  
  this->matrix_for_interpolate.allocate(1, N, field);
  this->matrix_for_interpolate.InitHIMByVectors(alpha, beta);
  this->matrix_vand.allocate(N, N, field);
  this->matrix_vand.InitVDM(); // default to alpha as defined
  this->matrix_vand_transpose.allocate(N, N, field);
  this->matrix_vand_transpose.InitVDMTranspose();
}

template <class FieldType> bool ProtocolParty<FieldType>::preparationPhase() {
  int iterations =
      (5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();
  int keysize = 16 / field->getElementSizeInBytes() + 1;
  // TODO: this is loose
  int numOfRandomShares = 6 * keysize + 3 * iterations + 2;
  this->randomSharesArray.resize(numOfRandomShares);

  // -- generate random shares for the AES key
  this->randomSharesOffset = 0;
  generateRandomShares(numOfRandomShares, this->randomSharesArray);
  
  // -- generate t-shareing and 2t-sharings for DN mult to use
  this->offset = 0;
  int numDoubleShares = this->numOfMultGates*2 + 1;
  offlineDNForMultiplication(numDoubleShares);
  return true;
}

template <class FieldType>
bool ProtocolParty<FieldType>::checkConsistency(vector<FieldType> &x, int d) {
  if (d == T) {
    vector<FieldType> y(N - 1 - d); // the result of multiplication
    vector<FieldType> x_until_t(T + 1);

    for (int i = 0; i < T + 1; i++) {
      x_until_t[i] = x[i];
    }

    matrix_for_t.MatrixMult(x_until_t, y);

    // compare that the result is equal to the according positions in x
    for (int i = 0; i < N - d - 1; i++)
    {
      if ((y[i]) != (x[d + 1 + i])) {
        return false;
      }
    }
    return true;
  } else if (d == 2 * T) {
    vector<FieldType> y(N - 1 - d); // the result of multiplication

    vector<FieldType> x_until_2t(2 * T + 1);

    for (int i = 0; i < 2 * T + 1; i++) {
      x_until_2t[i] = x[i];
    }

    matrix_for_2t.MatrixMult(x_until_2t, y);

    // compare that the result is equal to the according positions in x
    for (int i = 0; i < N - d - 1; i++)
    {
      if ((y[i]) != (x[d + 1 + i])) {
        return false;
      }
    }
    return true;

  } else {
    vector<FieldType> alpha_until_d(d + 1);
    vector<FieldType> alpha_from_d(N - 1 - d);
    vector<FieldType> x_until_d(d + 1);
    vector<FieldType> y(N - 1 - d); // the result of multiplication

    for (int i = 0; i < d + 1; i++) {
      alpha_until_d[i] = alpha[i];
      x_until_d[i] = x[i];
    }
    for (int i = d + 1; i < N; i++) {
      alpha_from_d[i - (d + 1)] = alpha[i];
    }
    // Interpolate first d+1 positions of (alpha,x)
    HIM<FieldType> matrix(N - 1 - d, d + 1,
                          field); // slices, only positions from 0..d
    matrix.InitHIMByVectors(alpha_until_d, alpha_from_d);
    matrix.MatrixMult(x_until_d, y);

    // compare that the result is equal to the according positions in x
    for (int i = 0; i < N - d - 1; i++)
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
    return interpolate(x);  
}

// Interpolate polynomial at position Zero
template <class FieldType>
FieldType ProtocolParty<FieldType>::interpolate(vector<FieldType> &x) {
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
      FieldType e = field->GetElement(scalar);
      gateShareArr[circuit.getGates()[k].output] =
          gateShareArr[circuit.getGates()[k].input1] * e;
      count++;
    } else if (circuit.getGates()[k].gateType == SCALAR_ADD) {
      long scalar(circuit.getGates()[k].input2);
      FieldType e = field->GetElement(scalar);
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

  
  int fieldByteSize = field->getElementSizeInBytes();
  int maxNumberOfLayerMult = circuit.getLayers()[currentCircuitLayer + 1] -
                             circuit.getLayers()[currentCircuitLayer];
  vector<FieldType> xyMinusRShares(maxNumberOfLayerMult);
  vector<FieldType> xyMinusR;
  vector<byte> xyMinusRBytes;
  vector<vector<byte>> recBufsBytes(N);
  vector<vector<byte>> sendBufsBytes(N);
  vector<vector<FieldType>> sendBufsElements(N);

  // -- go through current layer to collect mult gates
  int index = 0;
  for (int k = circuit.getLayers()[currentCircuitLayer];
       k < circuit.getLayers()[currentCircuitLayer + 1]; k++) {
    auto gate = circuit.getGates()[k];

    if (gate.gateType == MULT) {
      // compute the share of xy-r
      xyMinusRShares[index] =
          gateShareArr[gate.input1] * gateShareArr[gate.input2] -
          randomTAnd2TShares[2 * indexInRandomArray + 1];
      indexInRandomArray++;
      index++;
    }
  }

  // -- early stop if no mult in current layer
  if (index == 0)
    return 0;

  // -- spread reconstruction work to all parties
  // set the acctual number of mult gate proccessed in this layer
  int acctualNumOfMultGates = index;
  int numOfElementsForParties = acctualNumOfMultGates / N;
  int indexForDecreasingSize =
      acctualNumOfMultGates - numOfElementsForParties * N;

  int counter = 0;
  int currentNumOfElements;
  for (int i = 0; i < N; i++) {
    currentNumOfElements = numOfElementsForParties;
    if (i < indexForDecreasingSize)
      currentNumOfElements++;

    // fill the send buf according to the number of elements to send to each
    // party
    sendBufsElements[i].resize(currentNumOfElements);
    sendBufsBytes[i].resize(currentNumOfElements * fieldByteSize);
    for (int j = 0; j < currentNumOfElements; j++) {
      sendBufsElements[i][j] = xyMinusRShares[counter];
      field->elementToBytes(sendBufsBytes[i].data() + (j*fieldByteSize),
                            sendBufsElements[i][j]);
      counter++;
    }
  }

  // -- resize the recbuf array
  int myNumOfElementsToExpect = numOfElementsForParties;
  if (m_partyId < indexForDecreasingSize) {
    myNumOfElementsToExpect = numOfElementsForParties + 1;
  }
  for (int i = 0; i < N; i++) {
    recBufsBytes[i].resize(myNumOfElementsToExpect * fieldByteSize);
  }

  roundFunctionSync(sendBufsBytes, recBufsBytes, 8);

  // -- convert received bytes to shares, and reconstruct
  xyMinusR.resize(myNumOfElementsToExpect);
  xyMinusRBytes.resize(myNumOfElementsToExpect * fieldByteSize);
  vector<FieldType> xyMinurAllShares(N);
  for (int k = 0; k < myNumOfElementsToExpect; k++) {
    for (int i = 0; i < N; i++) {
      xyMinurAllShares[i] =
          field->bytesToElement(recBufsBytes[i].data() + (k * fieldByteSize));
    }
    xyMinusR[k] = interpolate(xyMinurAllShares);
    field->elementToBytes(xyMinusRBytes.data() + (k * fieldByteSize),
                          xyMinusR[k]);
  }

  // -- send reconstruct results to all, and receive from all
  for (int i = 0; i < N; i++) {
    sendBufsBytes[i] = xyMinusRBytes;
  }
  for (int i = 0; i < N; i++) {
    currentNumOfElements = numOfElementsForParties;
    if (i < indexForDecreasingSize)
      currentNumOfElements++;

    recBufsBytes[i].resize(currentNumOfElements * fieldByteSize);
  }
  roundFunctionSync(sendBufsBytes, recBufsBytes, 9);

  // -- convert byte to elements, reusing xyMinusR buffer
  xyMinusR.resize(acctualNumOfMultGates);
  counter = 0;
  for (int i = 0; i < N; i++) {
    currentNumOfElements = numOfElementsForParties;
    if (i < indexForDecreasingSize)
      currentNumOfElements++;

    for (int j = 0; j < currentNumOfElements; j++) {
      xyMinusR[counter] =
          field->bytesToElement(recBufsBytes[i].data() + (j * fieldByteSize));
      counter++;
    }
  }

  // -- store mult results to gate wire array
  indexInRandomArray -= index;
  index = 0;
  for (int k = circuit.getLayers()[currentCircuitLayer];
       k < circuit.getLayers()[currentCircuitLayer + 1]; k++) {
    auto gate = circuit.getGates()[k];

    if (gate.gateType == MULT) {
      gateShareArr[gate.output] =
          randomTAnd2TShares[2 * indexInRandomArray] + xyMinusR[index];
      indexInRandomArray++;
      index++;
    }
  }

  return index;
}

// template <class FieldType>
// void ProtocolParty<FieldType>::DNHonestMultiplication(
//     FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfMults) {

//   int fieldByteSize = field->getElementSizeInBytes();
//   // hold both in the same vector to send in one batch
//   vector<FieldType> xyMinusRShares(numOfMults);
//   vector<byte> xyMinusRSharesBytes(numOfMults * fieldByteSize);
//   vector<FieldType> xyMinusR(numOfMults);
//   vector<byte> xyMinusRBytes(numOfMults * fieldByteSize);
//   vector<vector<byte>> recBufsBytes;

//   // -- generate the 2t-sharings for xy - r
//   for (int k = 0; k < numOfMults; k++) {
//     xyMinusRShares[k] = a[k] * b[k] - randomTAnd2TShares[offset + 2 * k + 1];
//   }

//   for(int j=0; j<numOfMults;j++) {
//     field->elementToBytes(xyMinusRSharesBytes.data() +
//                           (j * fieldByteSize), xyMinusRShares[j]);
//   }

//   // -- gather from all to p0
//   if (m_partyId == 0) {
//     // p0 receive the shares from all the other parties
//     recBufsBytes.resize(N);
//     for (int i = 0; i < N; i++) {
//       recBufsBytes[i].resize(numOfMults * fieldByteSize);
//     }
//     recToP1(xyMinusRSharesBytes, recBufsBytes);
//   } else {
//     // send the shares to p0
//     parties[0]->getChannel()->write(xyMinusRSharesBytes.data(),
//                                     xyMinusRSharesBytes.size());
//   }

//   // -- p0 reconstruct the shares and send to all
//   if (m_partyId == 0) {
//     vector<FieldType> xyMinurAllShares(N);

//     for (int k = 0; k < numOfMults; k++) {
//       for (int i = 0; i < N; i++) {
//         xyMinurAllShares[i] =
//             field->bytesToElement(recBufsBytes[i].data() + (k * fieldByteSize));
//       }
//       // reconstruct the shares by p0
//       xyMinusR[k] = interpolate(xyMinurAllShares);
//       // convert to bytes
//       field->elementToBytes(xyMinusRBytes.data() + (k * fieldByteSize),
//                             xyMinusR[k]);
//     }

//     // send the reconstructed vector to all the other parties
//     sendFromP1(xyMinusRBytes);
//   } else { 
//     // each party get the xy-r reconstruced vector from p0
//     parties[0]->getChannel()->read(xyMinusRBytes.data(), xyMinusRBytes.size());
//   }

//   // -- other party convert by to element
//   if (m_partyId != 0) {
//     for (int i = 0; i < numOfMults; i++) {
//       xyMinusR[i] = field->bytesToElement(xyMinusRBytes.data() +
//                                           (i * fieldByteSize));
//     }
//   }

//   // -- add the xy-r bytes to h (accumulated view)
//   h.insert(h.end(), xyMinusRBytes.begin(), xyMinusRBytes.end());
  
//   // -- compute xy - r + [r]_t = t-sharing of xy
//   for (int k = 0; k < numOfMults; k++) {
//     cToFill[k] = randomTAnd2TShares[offset + 2 * k] + xyMinusR[k];
//   }
//   offset += numOfMults * 2;
// }

template <class FieldType>
void ProtocolParty<FieldType>::
offlineDNForMultiplication(int numOfDoubleShares) {
  generateRandom2TAndTShares(numOfDoubleShares,
                             this->randomTAnd2TShares);
}

template <class FieldType>
void ProtocolParty<FieldType>::openShare(int numOfRandomShares,
                                         vector<FieldType> &Shares,
                                         vector<FieldType> &secrets) {

  vector<vector<byte>> sendBufsBytes(N);
  vector<vector<byte>> recBufsBytes(N);

  vector<FieldType> x1(N);
  int fieldByteSize = field->getElementSizeInBytes();

  // resize vectors
  for (int i = 0; i < N; i++) {
    sendBufsBytes[i].resize(numOfRandomShares * fieldByteSize);
    recBufsBytes[i].resize(numOfRandomShares * fieldByteSize);
  }

  // set the first sending data buffer
  for (int j = 0; j < numOfRandomShares; j++) {
    field->elementToBytes(sendBufsBytes[0].data() + (j * fieldByteSize),
                          Shares[j]);
  }

  // copy the same data for all parties
  for (int i = 1; i < N; i++) {
    sendBufsBytes[i] = sendBufsBytes[0];
  }

  roundFunctionSync(sendBufsBytes, recBufsBytes, 6);

  // reconstruct each set of shares to get the secret

  for (int k = 0; k < numOfRandomShares; k++) {
    // get the set of shares for each element
    for (int i = 0; i < N; i++) {
      x1[i] = field->bytesToElement(recBufsBytes[i].data() +
                                    (k * fieldByteSize));
    }
    secrets[k] = reconstructShare(x1, T);
  }
}

/**
 * the function Walk through the circuit and reconstruct output gates.
 * @param circuit
 * @param gateShareArr
 * @param alpha
 */
template <class FieldType> void ProtocolParty<FieldType>::outputPhase() {
  vector<FieldType> x1(N); // vector for the shares of my outputs
  vector<vector<FieldType>> sendBufsElements(N);
  vector<vector<byte>> sendBufsBytes(N);
  vector<vector<byte>> recBufBytes(N);

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

  int fieldByteSize = field->getElementSizeInBytes();
  for (int i = 0; i < N; i++) {
    sendBufsBytes[i].resize(sendBufsElements[i].size() * fieldByteSize);
    recBufBytes[i].resize(sendBufsElements[m_partyId].size() * fieldByteSize);
    for(int j=0; j<sendBufsElements[i].size();j++) {
      field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize),
                            sendBufsElements[i][j]);
    }
  }

  roundFunctionSync(sendBufsBytes, recBufBytes, 11);

  int counter = 0;
  for (int k = M - numOfOutputGates; k < M; k++) {
    if (circuit.getGates()[k].gateType == OUTPUT &&
        circuit.getGates()[k].party == m_partyId) {
      for (int i = 0; i < N; i++) {
        x1[i] = field->bytesToElement(recBufBytes[i].data() +
                                      (counter * fieldByteSize));
      }

      if (flag_print_output)
        cout << "the result for " << circuit.getGates()[k].input1
             << " is : " << field->elementToString(interpolate(x1)) << '\n';
      counter++;
    }
  }

  // close output file
  myfile.close();
}

template <class FieldType>
void ProtocolParty<FieldType>::roundFunctionSync(vector<vector<byte>> &sendBufs,
                                                 vector<vector<byte>> &recBufs,
                                                 int round) {
  // -- assign works to threads
  int numPartiesForEachThread;
  if (parties.size() <= numThreads) {
    numThreads = parties.size();
    numPartiesForEachThread = 1;
  } else {
    numPartiesForEachThread = (parties.size() + numThreads - 1) / numThreads;
  }
  
  // -- recieve the data using threads
  // direct send to myself
  recBufs[m_partyId] = move(sendBufs[m_partyId]);
  vector<thread> threads(numThreads);
  for (int t = 0; t < numThreads; t++) {
    if ((t + 1) * numPartiesForEachThread <= parties.size()) {
      threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs),
                          ref(recBufs), t * numPartiesForEachThread,
                          (t + 1) * numPartiesForEachThread);
    } else {
      threads[t] =
          thread(&ProtocolParty::exchangeData, this, ref(sendBufs),
                 ref(recBufs), t * numPartiesForEachThread, parties.size());
    }
  }
  for (int t = 0; t < numThreads; t++) {
    threads[t].join();
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::exchangeData(vector<vector<byte>> &sendBufs,
                                            vector<vector<byte>> &recBufs,
                                            int first, int last) {

  // cout<<"in exchangeData";
  for (int i = first; i < last; i++) {

    if ((m_partyId) < parties[i]->getID()) {

      // send shares to my input bits
      parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(),
                                      sendBufs[parties[i]->getID()].size());
      // cout<<"write the data:: my Id = " << m_partyId - 1<< "other ID = "<<
      // parties[i]->getID() <<endl;

      // receive shares from the other party and set them in the shares array
      parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(),
                                     recBufs[parties[i]->getID()].size());
      // cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<<
      // parties[i]->getID()<<endl;

    } else {

      // receive shares from the other party and set them in the shares array
      parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(),
                                     recBufs[parties[i]->getID()].size());
      // cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<<
      // parties[i]->getID()<<endl;

      // send shares to my input bits
      parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(),
                                      sendBufs[parties[i]->getID()].size());
      // cout<<"write the data:: my Id = " << m_partyId-1<< "other ID = "<<
      // parties[i]->getID() <<endl;
    }
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::recToP1(
    vector<byte> &myShare, vector<vector<byte>> &recBufs) {
  
  int numPartiesForEachThread;
  if (parties.size() <= numThreads) {
    numThreads = parties.size();
    numPartiesForEachThread = 1;
  } else {
    numPartiesForEachThread = (parties.size() + numThreads - 1) / numThreads;
  }

  recBufs[m_partyId] = myShare;
  // recieve the data using threads
  vector<thread> threads(numThreads);
  for (int t = 0; t < numThreads; t++) {
    if ((t + 1) * numPartiesForEachThread <= parties.size()) {
      threads[t] = thread(&ProtocolParty::recDataToP1, this, ref(recBufs),
                          t * numPartiesForEachThread,
                          (t + 1) * numPartiesForEachThread);
    } else {
      threads[t] = thread(&ProtocolParty::recDataToP1, this, ref(recBufs),
                          t * numPartiesForEachThread, parties.size());
    }
  }
  for (int t = 0; t < numThreads; t++) {
    threads[t].join();
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::recDataToP1(vector<vector<byte>> &recBufs,
                                           int first, int last) {
  for (int i = first; i < last; i++) {
    parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(),
                                   recBufs[parties[i]->getID()].size());
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::sendFromP1(vector<byte> &sendBuf) {
  
  int numPartiesForEachThread;
  if (parties.size() <= numThreads) {
    numThreads = parties.size();
    numPartiesForEachThread = 1;
  } else {
    numPartiesForEachThread = (parties.size() + numThreads - 1) / numThreads;
  }

  // send data using threads
  vector<thread> threads(numThreads);
  for (int t = 0; t < numThreads; t++) {
    if ((t + 1) * numPartiesForEachThread <= parties.size()) {
      threads[t] = thread(&ProtocolParty::sendDataFromP1, this, ref(sendBuf),
                          t * numPartiesForEachThread,
                          (t + 1) * numPartiesForEachThread);
    } else {
      threads[t] = thread(&ProtocolParty::sendDataFromP1, this, ref(sendBuf),
                          t * numPartiesForEachThread, parties.size());
    }
  }
  for (int t = 0; t < numThreads; t++) {
    threads[t].join();
  }
}

template <class FieldType>
void ProtocolParty<FieldType>::sendDataFromP1(vector<byte> &sendBuf, int first,
                                              int last) {
  for (int i = first; i < last; i++) {
    parties[i]->getChannel()->write(sendBuf.data(), sendBuf.size());
  }
}

template <class FieldType> ProtocolParty<FieldType>::~ProtocolParty() {
  protocolTimer->writeToFile();
  delete protocolTimer;
  delete field;
  delete timer;
  // delete comm;
}

#endif /* PROTOCOLPARTY_H_ */
