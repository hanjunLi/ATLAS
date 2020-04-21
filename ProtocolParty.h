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


#include "Interpolate.h"
#include <cmath>

#define flag_print false
#define flag_print_output false

using namespace std;
using namespace std::chrono;

template<class FieldType> class CompareGate;

template <class FieldType>
class ProtocolParty : public Protocol, public HonestMajority, MultiParty {

	public:
		// -- polynomial functionalities
		Interpolate<FieldType> interp; // evaluate (O(n)), interpolate polynomials (O(n^2))
		vector<FieldType> gateShareArr; // my share of each gate (1 share each)
		int numOfCompareGates;
		int N, T;
		VDM<FieldType> matrix_vand;
		//inputs for helen circuit
		vector<vector<FieldType> > _Ai;
		vector<FieldType> _bi;
	private:
		// -- global const
		int numThreads = 2;      // TODO: add as main arguments later
		int _K = 4;              // interpolation degree <=> 'shrink' factor
		// -- global variables
		int iteration;                       // current iteration number
		int currentCircuitLayer = 0;         // current circuit layer
		vector<FieldType> y_for_interpolate; // size = N
		boost::asio::io_service io_service;

		// -- set in constructor
		int times;                // number of times to run the run function
		Measurement *timer;
		TemplateField<FieldType> *field; // used to call some field functions
		ProtocolTimer *protocolTimer;
		CompareGate<FieldType> *comper;
		// N, T - number of parties, malicious
		int keySize, verifyIterations, m_partyId, eleSize;
		string s;                     // string of my_partyId
		string inputsFile, outputFile;
		ArithmeticCircuit circuit;
		// M - number of gates
		int M, numOfInputGates, numOfOutputGates, numOfMultGates;// numOfCompareGates;
		vector<long> myInputs;
		vector<shared_ptr<ProtocolPartyData>> parties;

		FieldType INV2;
		// -- set in initialization phase
		vector<FieldType> beta;      // a single zero element
		vector<FieldType> alpha;     // N distinct non-zero field elements
		vector<FieldType> _alpha_k;  // K distinct non-zero field elements
		vector<FieldType> _alpha_2k; // 2K-1 distinct non-zero field elements
		HIM<FieldType> matrix_for_interpolate;
		HIM<FieldType> matrix_for_t;  // interpolate first t+1 points to rest N-t-1
		HIM<FieldType> matrix_for_2t; // interpolate first 2t+1 tpoints to rest N-2t-1
		HIM<FieldType> matrix_for_k;  // interpolate first k points to rest 2k-1 - k
		VDMTranspose<FieldType> matrix_vand_transpose;  
		
		// -- filled during offline-preparation phase
		vector<FieldType> _singleSharesArray; // used as seed for AES PRG
		int _singleSharesOffset = 0;          // next random share to use
		vector<FieldType> _doubleSharesArray; // used for DN mult
		int _doubleSharesOffset = 0;    // next random double share to use
		vector<FieldType> gateValueArr; // the value for my input and output
		//shares of [a]_b and [a]_p, will be initialized in preparation phase
		
	public:
		bool hasOffline() { return true; }
		bool hasOnline() override { return true; }
		~ProtocolParty();

		// -- protocol functions
		ProtocolParty(int argc, char *argv[]);
		void readMyInputs();
	
		//reading for helen circuit
		void readLassoInputs();


		void initializationPhase();

		void run() override;
		void runOffline() override;
		bool preparationPhase();
		void offlineDNForMultiplication(int numOfDoubleShares);

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
		FieldType randomCoin();
		
		int processComp();
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
		FieldType interpolate_to_zero(vector<FieldType> &x);
};


template<class FieldType>
class CompareGate 
{
private:
	ProtocolParty<FieldType> *helper;
	int eleSize;
	vector<vector<FieldType> > OrVector; //the share of coefficients of fan-in Or function
	vector<FieldType> _bitSharesValue;
	vector<vector<FieldType> > _bitSharesBits;
	vector<FieldType> _zeroShares;
	int _zeroShareOffset = 0;
	int _bitShareOffset = 0;
	//for safety, we use 12 bit floats
	//fixed-point floats, 4 decimals
	int _k = 24;
	int _m = 8;
	int _kappa = 8;
	int m_partyID;
	TemplateField<FieldType> *field; // used to call some field functions
public:
	CompareGate(ProtocolParty<FieldType> *ptr,int siz,int m_partyID,TemplateField<FieldType> *field);
	~CompareGate();

	// Compgate functions
	void getRandomBitShare(int num,vector<FieldType> &res, vector<vector<FieldType> > &bits);
	//void genRandom01(int num,vector<vector<FieldType>> &bin,vector<FieldType> &val);
	void compRandom(vector<FieldType> &a,vector<FieldType> &b,vector<FieldType> &res);
	//void compRandom(vector<FieldType> &a, vector<FieldType> &b, vector<FieldType> &dest);
	// given known element c and binary share a, store (c<a) to dest
	void compGiven(vector<vector<FieldType> > &a, vector<FieldType> &c, vector<FieldType> &dest);
	void compGivenNoPrefixOr(vector<vector<FieldType> > &a, vector<FieldType> &c, vector<FieldType> &dest);
	void computeAnswer(vector<vector<FieldType> > &b, vector<vector<FieldType> > &c,vector<FieldType> &res);
	void compP(vector<vector<FieldType> > &a, vector<FieldType> &dest);
	// compare the function (a < p/2)
	void compHalfP(vector<FieldType> &a, vector<FieldType> &dest);
	// compute square inverse of a given element
	FieldType compSqrInverse(FieldType &a);
	//save the prefix or of array a to array res
	void compPrefixOr(vector<vector<FieldType> > &a, vector<vector<FieldType> > &res);
	void compFanInOr(int num, vector<vector<FieldType> > &a, vector<FieldType> &res);
	int generateBitShares(int numOfBitShares);
	void getRandomInvPairs(int num, vector<FieldType> &b, vector<FieldType> &invb, vector<FieldType> &tmpc);
	void getBitDecomp(FieldType &e, vector<FieldType> &dest);
	///////////////////////////////////////
	//only used for multiplication prod. Par: m and k
	void TruncPR(vector<FieldType> &a,vector<FieldType> &res);//vector<int> &k,vector<int> &m);
	void doubleVecMult(vector<FieldType> &a,vector<FieldType> &b,vector<FieldType> &res);
	void SoftThres(vector<FieldType> &thres, vector<FieldType> &a, vector<FieldType> &res);
	void runLasso(int iter,FieldType lambda, FieldType rho, vector<vector<FieldType> > & Ai, vector<FieldType> &bi, vector<FieldType> &res);
	//compute 1/a[i] under double
	void doubleInverse(FieldType a,FieldType &res);
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
			this->field = new TemplateField<FieldType>((1ll<<61)-1);
		}
		//save the value of inversion of 2
		//this->INV2 = this->field->inv(2);
		string circuitFile =
			this->getParser().getValueByKey(arguments, "circuitFile");
		string outputTimerFileName =
			circuitFile + "Times" + to_string(m_partyId) + fieldType + ".csv";
		this->protocolTimer = new ProtocolTimer(times, outputTimerFileName);
		this->eleSize = field->getElementSizeInBytes() * 8;
		if(flag_print)
			cout<<"flag size:"<<eleSize<<endl;
		this->N = stoi(this->getParser().getValueByKey(arguments, "partiesNumber"));
		this->T = (this->N - 1)/ 2;
		this->keySize = 16 / field->getElementSizeInBytes() + 1;
		this->verifyIterations =
			(5 + field->getElementSizeInBytes() - 1) / field->getElementSizeInBytes();
		this->m_partyId = stoi(this->getParser().getValueByKey(arguments, "partyID"));
		this->s = to_string(m_partyId);
		this->inputsFile = this->getParser().getValueByKey(arguments, "inputFile");
		this->outputFile = this->getParser().getValueByKey(arguments, "outputFile");
		this->circuit.readCircuit(circuitFile.c_str());
		this->circuit.reArrangeCircuit();
		this->M = circuit.getNrOfGates();
		this->numOfInputGates = circuit.getNrOfInputGates();
		this->numOfOutputGates = circuit.getNrOfOutputGates();
		//should we +1 to avoid log(0)?
		this->numOfMultGates = circuit.getNrOfMultiplicationGates() + 1;
		this->numOfCompareGates = circuit.getNrOfCompareGates();
		this->myInputs.resize(numOfInputGates);
		cout<<"# of inputs:"<<numOfInputGates<<endl;
		readLassoInputs();
		//readMyInputs();
		MPCCommunication comm;
		string partiesFile =
			this->getParser().getValueByKey(arguments, "partiesFile");
		if(flag_print)
			cout<<"Ready for Communication"<<endl;
		this->parties =
			comm.setCommunication(this->io_service, m_partyId, N, partiesFile);
		if(flag_print)
			cout<<"Set Communication"<<endl;
		// -- test communication
		string tmp = "init times";
		byte tmpBytes[20];
		for (int i = 0; i < parties.size(); i++) {
			if(flag_print)
				cout<<"Comm attempt #"<<i<<endl;
			if (parties[i]->getID() < m_partyId) {
				parties[i]->getChannel()->write(tmp);
				parties[i]->getChannel()->read(tmpBytes, tmp.size());
			} else {
				parties[i]->getChannel()->read(tmpBytes, tmp.size());
				parties[i]->getChannel()->write(tmp);
			}
		}

		// -- initialize phase
		initializationPhase();
	}

template <class FieldType> void ProtocolParty<FieldType>::readMyInputs() {

	ifstream myfile;
	long input;
	int i = 0;
	myfile.open(this->inputsFile);
	while(myfile>>input)
	{
		//myfile >> input;
		//cout <<"input read "<<i<<":"<<input<<endl;
		if (i >= myInputs.size()){
			cout << "inputs file " << this->inputsFile <<endl;
			cout << "have more than " << (i+1) << " inputs!" << endl;
			abort();
		}

		myInputs[i++] = input;
	} //while (!(myfile.eof()));
	if(flag_print)
		cout<<"Read Done!"<<endl;
	myfile.close();
}
template <class FieldType> void ProtocolParty<FieldType>::readLassoInputs()
{
	if(flag_print)
		cout<<"Reading lasso"<<endl;
	ifstream myfile;
	myfile.open(this->inputsFile);
	int n; myfile>>n;
	if(flag_print)
		cout<<"Read "<<n<<endl;
	//read Ai
	_Ai.resize(n);
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
		{
			long tmp;
			myfile>>tmp;
			if(flag_print)
				cout<<"Read "<<tmp<<","<<FieldType(tmp)<<endl;
			_Ai[i].push_back(FieldType(tmp));
		}
	//read Bi
	for(int i=0; i<n; i++)
	{
		long tmp;
		myfile>>tmp;
		if(flag_print)
			cout<<"Read "<<tmp<<endl;
		_bi.push_back(FieldType(tmp));
	}
	myfile.close();
	if(flag_print)
		cout<<"Reading lasso ended"<<endl;
}

template <class FieldType> void ProtocolParty<FieldType>::run() {
	int tottme = 0;
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
		tottme += duration;
	}
	cout<<"Total time:"<<tottme<<endl;
	cout<<"Gates:"<<numOfMultGates + numOfCompareGates<<endl;
	cout<<"Ave time:"<<(double)tottme / times / (numOfMultGates + numOfCompareGates);
}

template <class FieldType> void ProtocolParty<FieldType>::runOffline() {
	if(flag_print)
		cout<<"runOffline()"<<endl;
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
	//ONLY FOR TEST PURPOSE
	//verificationPhase();
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
	vector<FieldType> res;
	if(flag_print)
		cout<<"running lasso"<<endl;
	comper->runLasso(10,FieldType(1ll<<(7)),FieldType(1ll<<(7)),_Ai,_bi,res);
	//if(flag_print)
	{
		cout<<"Lasso finished"<<endl;
	cout<<"Used single:"<<_singleSharesOffset<<"/"<<_singleSharesArray.size()<<","<<(double)_singleSharesOffset / _singleSharesArray.size()<<",left:"<<_singleSharesArray.size() - _singleSharesOffset << endl;
		cout<<"Used double:"<<2 * _doubleSharesOffset<<"/"<<_doubleSharesArray.size()<<","<<(double)2*_doubleSharesOffset / _doubleSharesArray.size()<<",left:"<<-2*_doubleSharesOffset + _doubleSharesArray.size()<<endl;
	}
	abort();
	int countNumMult = _doubleSharesOffset;
	int numOfLayers = circuit.getLayers().size();
	for (int i = 0; i < numOfLayers - 1; i++) {
		if(flag_print)
			cout<<"Now at layer "<<i<<"/"<<numOfLayers<<endl;
		currentCircuitLayer = i;
		// send the index of the current mult gate
		processNotMult();
		processComp();
		if(flag_print)
			cout<<countNumMult<<endl;
		//TODO: just for test
		countNumMult += processMultiplications(countNumMult);
		if(flag_print)
			cout<<"layer finished"<<endl;

	}
	_doubleSharesOffset = countNumMult;
}

template <class FieldType>
void ProtocolParty<FieldType>::verificationPhase() {
	if(flag_print)
		cout<<"Now at Verification"<<endl;
	vector<FieldType> aShares(this->numOfMultGates);
	vector<FieldType> bShares(this->numOfMultGates);
	FieldType cShare = *(field->GetZero());

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
	if(flag_print)
		cout<<"Ready for Verify"<<endl;
	compressVerifyVec(aShares, bShares, cShare);
	if(flag_print)
		cout<<"Verify end"<<endl;
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
		if(flag_print)
			cout<<"Base case"<<endl;
		verifyVec(aShares, bShares, cShare);
		return;
	} // else : recursive case:
	// -- divide into K groups
	int groupSize = (totalLength + _K -1 )/ _K;

	// cout << "---- cur groupSize: " << groupSize << endl;

	// -- one DN mult for each group i to compute   
	// build dShares 0 .. k-2
	vector<FieldType> dShares;
	DNMultVec(aShares, bShares, dShares, groupSize);
	if(flag_print)
		cout<<"1st DNMult end"<<endl;
	// build dShares k-1: c - d_0 - ... - d_{k-2}
	dShares[_K-1] = cShare;
	for (int i = 0; i< _K - 1; i++) {
		dShares[_K-1] = dShares[_K-1] - dShares[i];
	}
	aShares.resize( groupSize * _K, field->GetElement(0) );
	bShares.resize( groupSize * _K, field->GetElement(0) );
	buildPolyVecInd(aShares, bShares, dShares, groupSize);
	if(flag_print)
		cout<<"Ind built"<<endl;
	vector<FieldType> aSharesNew(groupSize);
	vector<FieldType> bSharesNew(groupSize);
	if(flag_print)
		cout<<"At middle"<<endl;
	FieldType lambda = randomCoin();
	if(flag_print)
		cout<<"Get randomcoin"<<endl;
	HIM<FieldType> matrix_for_k_lambda;
	HIM<FieldType> matrix_for_2k_lambda;
	vector<FieldType> beta_lambda(1);
	beta_lambda[0] = lambda;
	matrix_for_k_lambda.allocate(1, _K, field);
	matrix_for_k_lambda.InitHIMByVectors(_alpha_k, beta_lambda);
	matrix_for_2k_lambda.allocate(1, 2*_K-1, field);
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
	getRandomSharesWithCheck(vecSize, aN);
	vector<FieldType> bN(vecSize);
	getRandomSharesWithCheck(vecSize, bN);

	aShares.insert(aShares.end(), aN.begin(), aN.end());
	bShares.insert(bShares.end(), bN.begin(), bN.end());
	vector<FieldType> dShares(2);
	DNMultVec(aShares, bShares, dShares, vecSize);

	cShare += dShares[1];

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
	if (CResult[0] != this->field->GetElement(0)) {
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

	int fieldByteSize = field->getElementSizeInBytes();
	vector<FieldType> x1(N), y1(N);
	vector<vector<FieldType>> sendBufsElements(N);
	vector<vector<byte>> sendBufsBytes(N);
	vector<vector<byte>> recBufBytes(N);
	vector<vector<FieldType>> recBufElements(N);

	// -- prepare the shares for the input
	// also need to prepare shares for P and P/2
	int index = 0;
	vector<int> sizes(N);
	if(flag_print)
		cout<<"#input:"<<numOfInputGates<<",#output:"<<numOfOutputGates<<endl;
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
	/*
	   if(this->m_partyId == 0) //send info about P and P/2
	   {
	//first generate shares. For convenience, suppose P = 2^k-1
	for(int i=0; i<2*fieldByteSize-1; i++)
	{
	x1[0] = FieldType(1);
	for (int j=1; j<T+1; j++)
	x1[i] = field->Random();
	matrix_vand.MatrixMult(x1,y1,T+1);
	for(int j=0; j<N; j++)
	sendBufsElements[j].push_back(y1[j]);
	}
	x1[0] = FieldType(0);
	for(int j=1; j<T+1; j++)
	x1[j] = field->Random();
	matrix_vand.MatrixMult(x1,y1,T+1);
	for(int j=0; j<N; j++)
	sendBufsElements[j].push_back(y1[j]);
	}
	else
	sizes[0]+=2*fieldByteSize;
	 */
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
	vector<FieldType> inputShares(this->numOfInputGates + 2*fieldByteSize);
	for (int k = 0; k < this->numOfInputGates; k++) {
		if (circuit.getGates()[k].gateType == INPUT) {
			auto share = recBufElements[circuit.getGates()[k].party]
				[counters[circuit.getGates()[k].party]];
			counters[circuit.getGates()[k].party] += 1;
			gateShareArr[circuit.getGates()[k].output] = share;
			inputShares[k] = share;
		}
	}
	/*
	//the next bytes are for P and P/2
	_PBits.resize(N);
	_HalfPBits.resize(N);
	for(int i=0; i<fieldByteSize; i++)
	{
	auto share = recBufElements[0][counters[0]];
	counters[0]++;
	_PBits[i] = share;
	inputShares[this->numOfInputGates + i] = share;
	}
	for(int i=0; i<fieldByteSize; i++)
	{
	auto share = recBufElements[0][counters[0]];
	counters[0]++;
	_HalfPBits[i] = share;
	inputShares[this->numOfInputGates + fieldByteSize + i]=share;
	}
	 */
	inputVerification(inputShares);
}

template <class FieldType>
void ProtocolParty<FieldType>::inputVerification(
		vector<FieldType> &inputShares) {
	batchConsistencyCheckOfShares(inputShares);
}

template <class FieldType>
void ProtocolParty<FieldType>::generateRandomShares(
		int numOfRandoms, vector<FieldType> &randomElementsToFill) {
	// -- allocate communication buffers
	vector<FieldType> x1(N), y1(N), t1(N), r1(N);
	// every round N-T random shares are generated.
	int fieldByteSize = field->getElementSizeInBytes();
	int no_buckets = (numOfRandoms / (N - T)) + 1;
	if(flag_print)
		cout<<"Resizing "<<no_buckets<<","<<N<<","<<T<<","<<no_buckets * (N-T) <<endl;
	randomElementsToFill.resize(no_buckets * (N - T));
	// N vectors of no_buckets size
	vector<vector<FieldType>> sendBufsElements(N);
	// N vectors of no_buckets elements in bytes
	vector<vector<byte>> sendBufsBytes(N);
	vector<vector<byte>> recBufsBytes(N);
	if(flag_print)
		cout<<"Resizing"<<endl;
	for (int i = 0; i < N; i++) {
		sendBufsElements[i].resize(no_buckets);
		sendBufsBytes[i].resize(no_buckets * fieldByteSize);
		recBufsBytes[i].resize(no_buckets * fieldByteSize);
	}
	if(flag_print)
		cout<<"Resized"<<endl;

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

//modified, in case that not enough shares
template <class FieldType>
void ProtocolParty<FieldType>::getRandomShares(
		int numOfRandoms, vector<FieldType> &randomElementsToFill) {
	if(randomElementsToFill.size()<numOfRandoms)
		randomElementsToFill.resize(numOfRandoms);
	if(_singleSharesOffset + numOfRandoms > _singleSharesArray.size())
	{
		//if(flag_print)
		cout<<"Not enough single elements!"<<endl;
		//int lft = _singleShareArray.size() - _singleSharesOffset;
		//randomElementsToFill.assign(_singleSharesArray.begin() + _singleSharesOffset,
		//                         _singleSharesArray.begin() + _singleSharesOffset +
		//                         lft);
		//numOfRandoms -= lft;
		generateRandomShares(numOfRandoms,_singleSharesArray);
		if(flag_print)
			cout<<"Genereated"<<endl;
		_singleSharesOffset = 0;
	}
	// since assign will reset the size of array, we use iteration instead
	for(int i=0; i<numOfRandoms; i++)
		randomElementsToFill[i] = _singleSharesArray[_singleSharesOffset++];
	//  randomElementsToFill.assign(_singleSharesArray.begin() + _singleSharesOffset,
	//                            _singleSharesArray.begin() + _singleSharesOffset +
	//                          numOfRandoms);

	// _singleSharesOffset += numOfRandoms;
}

template <class FieldType>
void ProtocolParty<FieldType>::getRandomSharesWithCheck(
		int numOfRandoms, vector<FieldType> &randomElementsToFill) {

	getRandomShares(numOfRandoms, randomElementsToFill);

	batchConsistencyCheckOfShares(randomElementsToFill);
}

template <class FieldType>
vector<byte> ProtocolParty<FieldType>::generateCommonKey() {
	// -- calc the number of elements needed for 128 bit AES key
	int fieldByteSize = field->getElementSizeInBytes();

	// -- allocate buffer for AES key
	vector<FieldType> randomSharesArray(keySize);
	vector<FieldType> aesArray(keySize);
	vector<byte> aesKey(keySize * fieldByteSize);

	// -- generate random shares for the AES key
	getRandomShares(keySize, randomSharesArray);
	openShare(keySize, randomSharesArray, aesArray);

	// -- turn the aes array into bytes to get the common aes key.
	for (int i = 0; i < keySize; i++) {
		for (int j = 0; j < keySize; j++) {
			field->elementToBytes(aesKey.data() + (j * fieldByteSize), aesArray[j]);
		}
	}
	aesKey.resize(16);
	return aesKey;
}

template <class FieldType>
void ProtocolParty<FieldType>::generatePseudoRandomElements(
		vector<byte> &aesKey, vector<FieldType> &randomElementsToFill,
		int numOfRandomElements) {

	int fieldSize = field->getElementSizeInBytes();
	int fieldSizeBits = field->getElementSizeInBits();
	bool isLongRandoms;
	int size;
	if (fieldSize > 4) {
		isLongRandoms = true;
		size = 8;
	} else {

		isLongRandoms = false;
		size = 4;
	}

	if (flag_print) {
		cout << "size is " << size << " for party : " << m_partyId << endl;
	}

	PrgFromOpenSSLAES prg((numOfRandomElements * size / 16) + 1);
	SecretKey sk(aesKey, "aes");
	prg.setKey(sk);

	for (int i = 0; i < numOfRandomElements; i++) {

		if (isLongRandoms)
			randomElementsToFill[i] = field->GetElement(
					((unsigned long)prg.getRandom64()) >> (64 - fieldSizeBits));
		else
			randomElementsToFill[i] = field->GetElement(prg.getRandom32());
	}
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
	vector<FieldType> randomElements(inputShares.size() * verifyIterations);
	generatePseudoRandomElements(key, randomElements, inputShares.size());

	for (int j = 0; j < verifyIterations; j++) {
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
	if(flag_print)
		cout<<"Entering Self Init"<<endl;
	// -- allocate spaces for protocol storage
	this->y_for_interpolate.resize(N);
	this->gateShareArr.resize(M - this->numOfOutputGates);

	// -- initialize protocol parameters
	this->beta.resize(1);
	this->beta[0] = field->GetElement(0); // zero of the field
	this->alpha.resize(N);
	for (int i = 0; i < N; i++) {
		this->alpha[i] = field->GetElement(i + 1);
	}
	_alpha_2k.resize(_K*2 - 1);
	for (int i=0; i< _K*2 - 1; i++) {
		_alpha_2k[i] = field->GetElement(i+1);
	}
	_alpha_k.resize(_K);
	for (int i=0; i< _K; i++) {
		_alpha_k[i] = field->GetElement(i+1);
	}

	// Interpolate first T+1 positions (deg = T)
	matrix_for_t.allocate(N - (1+T), T + 1, field); 
	matrix_for_t.InitHIMVectorAndsizes(alpha, T + 1, N - T - 1);

	// Interpolate first 2T+1 positions (deg = 2T)
	matrix_for_2t.allocate(N - (1+2*T), 2*T + 1, field);
	matrix_for_2t.InitHIMVectorAndsizes(alpha, 2*T + 1, N - 1 - 2*T);

	// Interpolate first K positions (deg = K-1)
	matrix_for_k.allocate(2*_K-1 - _K, _K, field);
	matrix_for_k.InitHIMVectorAndsizes(_alpha_2k, _K, 2*_K-1 - _K);

	this->matrix_for_interpolate.allocate(1, N, field);
	this->matrix_for_interpolate.InitHIMByVectors(alpha, beta);
	this->matrix_vand.allocate(N, N, field);
	this->matrix_vand.InitVDM(); // default to alpha as defined
	this->matrix_vand_transpose.allocate(N, N, field);
	this->matrix_vand_transpose.InitVDMTranspose();
}

template <class FieldType> bool ProtocolParty<FieldType>::preparationPhase() {
	// ---- # of random single shares:
	// 1. Generating AES keys and verification
	//    [keySize + verifyIterations]
	//    -- inputVerification() -> 1 time
	//    -- generateRandom2TAndTShares() -> 1 time
	//    -- verificationPhase() -> 2 times in base case
	//    in total: 4 * (keySize + verifyIterations)
	// 2. Padding vector product verification in base case
	//    [2 * vecSize] // which is at most _K
	//    in total: 2 * _K
	// 3. Random Coins
	//    [1]
	//    in total: nCompressions + 2
	// 4. Random number with bit representation
	//	3 * fieldSize (+ possibility may fail?
	if(flag_print)
		cout<<"entering preparation phase"<<endl;
	//cout<<"eleSize:"<<eleSize<<endl;
	//cout<<"#comp:"<<numOfCompareGates<<endl;
	int nCompressions = (int)(log(this->numOfMultGates) / log(_K) + 0.5);
	int numSingleShares =
		4 * (keySize + verifyIterations) + 2 * _K + nCompressions + 2
		+ 60 * eleSize *numOfCompareGates;
	//cout<<"should have num:"<<numSingleShares<<endl;
	numSingleShares += 100000;
	_singleSharesOffset = 0;
	if(flag_print)
		cout<<"Generating single Share"<<endl;
	generateRandomShares(numSingleShares, _singleSharesArray);
	if(flag_print)
		cout << "generated single Shares: " << numSingleShares << endl;

	// ---- # of random double shares:
	// 1. Compute Mult gates in the circuit
	//    -- in total: numOfMultGates + 279 * numOfComp(?)
	// 2. Compress Verification
	//    -- 2 * _K * (nCompressions + 1)
	int numDoubleShares =
		this->numOfMultGates + (nCompressions+1)*_K*2 + 100 * eleSize * this->numOfCompareGates;
	//used for Lasso?
	numDoubleShares += 100000;
	_doubleSharesOffset = 0;
	offlineDNForMultiplication(numDoubleShares);

	if(flag_print)
		cout << "generated doubles: " << numDoubleShares << endl;
	//a brute force way for generating Or polys are applied.
	//can be implemented in a more efficient way.

	//put it here to let bit share to be generated
	this->comper = new CompareGate<FieldType>(this,eleSize,m_partyId,field);
	//vector<FieldType> tmp,tmpnum;
	//int tot = eleSize * (N+2);
	/*OrVector.resize(tot);
	  tmp.resize(1);
	  tmpnum.resize(1);
	  tmp[0] = FieldType(0);
	  tmpnum[0] = FieldType(1);
	  if(flag_print)
	  cout<<"tot:"<<tot<<" Intering 0"<<endl;
	  interp.interpolate(tmpnum, tmp, OrVector[0]);
	  for(int i=1; i<4; i++)
	  {
	  if(flag_print)
	  cout<<"Intering "<<i<<endl;
	  tmp.push_back(FieldType(1));
	  tmpnum.push_back(FieldType(i+1));
	  interp.interpolate(tmpnum, tmp, OrVector[i]);
	  if(flag_print)
	  for(int j=0; j<=i; j++)
	  cout<<OrVector[i][j]<<","<<endl;
	  }*/
	return true;
}

//template <>
//ZpMersenneIntElement ProtocolParty<ZpMersenneIntElement>::compSqrInverse(ZpMersenneIntElement &num)


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
		/*else if (circuit.getGates()[k].gateType == COMPARE)
		  {
		//first test if we can work with random 01
		if(flag_print)
		{
		cout<<"find a compare, dealing"<<endl;
		vector<FieldType> d0,d1;
		openShare(2,gateShareArr,d1);
		cout<<"input:"<<d1[0]<<","<<d1[1]<<endl;
		//cout<<circuit.getGates()[k].input1<<","<<circuit.getGates()[k].input2<<endl;
		}
		compRandom(circuit.getGates()[k]);
		}*/
	}
	return count;
}

// notice we only deal with 1 compare gate one time
//template <class FieldType>
//int ProtocolParty<FieldType>::processComparsion(int 
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
	if(flag_print)
		cout<<"limit this layer:"<<maxNumberOfLayerMult<<endl;
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
			if(flag_print)
				cout<<"Find a mult Gate"<<endl;
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
		xyMinusR[k] = interpolate_to_zero(xyMinurAllShares);
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
		AShares[i].resize(polyDegree, this->field->GetElement(0));
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
void ProtocolParty<FieldType>::
DNMultVec(vector<FieldType>& a, vector<FieldType>& b,
		vector<FieldType>& c, int groupSize) {
	// assert( a.size() == b.size() );
	if(flag_print)
		cout<<"Entering DNMult, sizes:"<<a.size()<<","<<b.size()<<endl;
	int totalLength = a.size();
	int numOfMults = (totalLength + groupSize - 1) / groupSize;
	int fieldByteSize = field->getElementSizeInBytes();
	c.resize(numOfMults);
	vector<FieldType> xyMinusRShares(numOfMults);
	vector<byte> xyMinusRSharesBytes(numOfMults * fieldByteSize);
	vector<FieldType> xyMinusR(numOfMults);
	vector<byte> xyMinusRBytes(numOfMults * fieldByteSize);
	vector<vector<byte>> recBufsBytes;
	if(flag_print)
		cout<<"# of mult:"<<numOfMults<<endl;
	// -- generate the 2t-sharings for xy - r
	for (int group_i = 0; group_i < numOfMults; group_i++) {
		int group_end = (group_i + 1) * groupSize > totalLength ?
			totalLength : (group_i + 1) * groupSize;
		for (int i = group_i * groupSize; i < group_end; i++) {
			xyMinusRShares[group_i] += a[i] * b[i];
		}
		xyMinusRShares[group_i] = xyMinusRShares[group_i] -
			_doubleSharesArray[(_doubleSharesOffset + group_i)*2 + 1];
	}
	for (int group_i = 0; group_i < numOfMults; group_i++) {
		field->elementToBytes(xyMinusRSharesBytes.data() +
				(group_i * fieldByteSize), xyMinusRShares[group_i]);
	}
	if(flag_print)
		cout<<"generated 2t-sharing"<<endl;
	// -- gather from all to p0
	if (this->m_partyId == 0) {
		// p0 receive the shares from all the other parties
		recBufsBytes.resize(N);
		for (int i = 0; i < N; i++) {
			recBufsBytes[i].resize(numOfMults * fieldByteSize);
		}
		recToP1(xyMinusRSharesBytes, recBufsBytes);
	} else {
		// send the shares to p0
		this->parties[0]->getChannel()->write(xyMinusRSharesBytes.data(),
				xyMinusRSharesBytes.size());
	}
	if(flag_print)
		cout<<"p0 received"<<endl;
	// -- p0 reconstruct the shares and send to all
	if (this->m_partyId == 0) {
		vector<FieldType> xyMinurAllShares(N);

		for (int group_i = 0; group_i < numOfMults; group_i++) {
			for (int i = 0; i < N; i++) {
				xyMinurAllShares[i] =
					field->bytesToElement(recBufsBytes[i].data() +
							(group_i * fieldByteSize));
			}
			// reconstruct the shares by p0
			xyMinusR[group_i] = interpolate_to_zero(xyMinurAllShares);
			// convert to bytes
			field->elementToBytes(xyMinusRBytes.data() + (group_i * fieldByteSize),
					xyMinusR[group_i]);
		}

		// send the reconstructed vector to all the other parties
		sendFromP1(xyMinusRBytes);
	} else { 
		// each party get the xy-r reconstruced vector from p0
		this->parties[0]->getChannel()->read(xyMinusRBytes.data(),
				xyMinusRBytes.size());
	}
	//if(flag_print)
	//	cout<<"P0 reconstructed"<<endl;
	// -- other party convert by to element
	if (this->m_partyId != 0) {
		for (int group_i = 0; group_i < numOfMults; group_i++) {
			xyMinusR[group_i] = field->bytesToElement(xyMinusRBytes.data() +
					(group_i * fieldByteSize));
		}
	}

	// -- compute xy - r + [r]_t = t-sharing of xy
	for (int group_i = 0; group_i < numOfMults; group_i++) {
		c[group_i] =
			_doubleSharesArray[(_doubleSharesOffset + group_i)*2] + xyMinusR[group_i];
	}
	if(flag_print)
		cout<<"DN ended"<<endl;
	_doubleSharesOffset += numOfMults;
}

template <class FieldType>
void ProtocolParty<FieldType>::
offlineDNForMultiplication(int numOfDoubleShares) {
	generateRandom2TAndTShares(numOfDoubleShares,
			_doubleSharesArray);
}

template <class FieldType>
void ProtocolParty<FieldType>::openShare(int numOfRandomShares,
		vector<FieldType> &Shares,
		vector<FieldType> &secrets) {
	//	if(flag_print)
	//		cout<<"Entering openshare, num:"<<numOfRandomShares<<","<<Shares.size()<<endl;
	if(secrets.size()<Shares.size())
		secrets.resize(Shares.size());
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
	// if(flag_print)
	//	 cout<<"Waiting Sync"<<endl;
	roundFunctionSync(sendBufsBytes, recBufsBytes, 6);
	// if(flag_print)
	//	 cout<<"Sync ended"<<endl;

	// reconstruct each set of shares to get the secret

	for (int k = 0; k < numOfRandomShares; k++) {
		// get the set of shares for each element
		for (int i = 0; i < N; i++) {
			x1[i] = field->bytesToElement(recBufsBytes[i].data() +
					(k * fieldByteSize));
		}
		//if(flag_print)
		//	cout<<"Reconstructing #"<<k<<endl;
		secrets[k] = reconstructShare(x1, T);
	}
}

template <class FieldType>
FieldType ProtocolParty<FieldType>::randomCoin() {
	vector<FieldType> lambdaShare(1);
	vector<FieldType> lambda(1);
	getRandomShares(1, lambdaShare);
	openShare(1, lambdaShare, lambda);
	return lambda[0];
}

/**
 * the function Walk through the circuit and reconstruct output gates.
 * @param circuit
 * @param gateShareArr
 * @param alpha
 */
template <class FieldType> void ProtocolParty<FieldType>::outputPhase() {
	if(flag_print)
		cout<<"Entering Output"<<endl;
	vector<FieldType> x1(N); // vector for the shares of my outputs
	vector<vector<FieldType>> sendBufsElements(N);
	vector<vector<byte>> sendBufsBytes(N);
	vector<vector<byte>> recBufBytes(N);

	FieldType num;
	ofstream myfile;
	myfile.open(this->outputFile);
	if(flag_print)
		cout<<"#tot gate:"<<M<<endl;
	for (int k = M - numOfOutputGates; k < M; k++) {
		if (circuit.getGates()[k].gateType == OUTPUT) {
			if(flag_print)
				cout<<"FOUND OUTPUT GATE"<<endl;
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

			if (flag_print)
				cout << "the result for " << circuit.getGates()[k].input1
					<< " is : " << field->elementToString(interpolate_to_zero(x1)) << '\n';
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
	vector<thread> threads(numThreads-1);
	for (int t = 0; t < numThreads-1; t++) {
		threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs),
				ref(recBufs), t * numPartiesForEachThread,
				(t + 1) * numPartiesForEachThread);
	}
	exchangeData(sendBufs, recBufs,
			(numThreads-1) * numPartiesForEachThread, parties.size());
	for (int t = 0; t < numThreads-1; t++) {
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
	vector<thread> threads(numThreads-1);
	for (int t = 0; t < numThreads-1; t++) {
		threads[t] = thread(&ProtocolParty::recDataToP1, this, ref(recBufs),
				t * numPartiesForEachThread,
				(t + 1) * numPartiesForEachThread);
	}
	recDataToP1(recBufs, (numThreads-1) * numPartiesForEachThread,
			parties.size());

	for (int t = 0; t < numThreads-1; t++) {
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
	vector<thread> threads(numThreads-1);
	for (int t = 0; t < numThreads-1; t++) {
		threads[t] = thread(&ProtocolParty::sendDataFromP1, this, ref(sendBuf),
				t * numPartiesForEachThread,
				(t + 1) * numPartiesForEachThread);
	}

	sendDataFromP1(sendBuf, (numThreads-1) * numPartiesForEachThread, parties.size());

	for (int t = 0; t < numThreads-1; t++) {
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


	template<class FieldType>
int ProtocolParty<FieldType>::processComp()
{
	//step 1: go through gates
	int cnt = 0;
	vector<FieldType> a,b,res;
	for(int k = circuit.getLayers()[currentCircuitLayer];
			k<circuit.getLayers()[currentCircuitLayer + 1]; k++)
	{
		auto gate = circuit.getGates()[k];
		if(gate.gateType == COMPARE)
		{
			cnt++;
			a.push_back(gateShareArr[gate.input1]);
			b.push_back(gateShareArr[gate.input2]);	
			//tmp.push_back(gate);
		}
	}
	comper->compRandom(a,b,res);
	cnt = 0;
	for(int k = circuit.getLayers()[currentCircuitLayer];
			k<circuit.getLayers()[currentCircuitLayer + 1]; k++)
	{
		auto gate = circuit.getGates()[k];
		if(gate.gateType == COMPARE)
			gateShareArr[gate.output] = res[cnt++];
	}
	if(flag_print)
	{
		cout<<"compRandom ended"<<endl;
		cout<<"Used single:"<<_singleSharesOffset<<"/"<<_singleSharesArray.size()<<","<<(double)_singleSharesOffset / _singleSharesArray.size()<<",left:"<<_singleSharesArray.size() - _singleSharesOffset << endl;
		cout<<"Used double:"<<2 * _doubleSharesOffset<<"/"<<_doubleSharesArray.size()<<","<<(double)2*_doubleSharesOffset / _doubleSharesArray.size()<<",left:"<<-2*_doubleSharesOffset + _doubleSharesArray.size()<<endl;
	}
	return cnt;
}

template <class FieldType> ProtocolParty<FieldType>::~ProtocolParty() {
	protocolTimer->writeToFile();
	delete protocolTimer;
	delete field;
	delete timer;
	// delete comm;
}

/////////////////////////////////////////////
//START OF COMPAREGATE
	template<class FieldType>
CompareGate<FieldType>::CompareGate(ProtocolParty<FieldType> *ptr,int siz,int id,TemplateField<FieldType> *f)
{
	m_partyID = id;
	helper = ptr;
	_bitShareOffset=0;
	eleSize = siz;
	field = f;
	int numBitShares = helper->numOfCompareGates + 100;
	if(flag_print)
		cout<<"generating bit shares:"<<numBitShares<<endl;
	generateBitShares(numBitShares);
}

	template<class FieldType>
CompareGate<FieldType>::~CompareGate()
{
}

	template <class FieldType>
FieldType CompareGate<FieldType>::compSqrInverse(FieldType &num)
{
}

	template<>
inline ZZ_p CompareGate<ZZ_p>::compSqrInverse(ZZ_p &num)
{
	//the mod is 2147483647
	//we don't need to init ZZ (?)
	//if(flag_print)
	//	cout<<"Computing Sqr Inv of"<<num<<":"<<endl;;ZZ md = _ZZPRIME; //rep(_PRIME);
	ZZ md((1ll<<61)-1); //rep(_PRIME);
	//	cout<<"Prime:"<<md<<endl;
	ZZ tmp = (md+1)/4;
	//if(flag_print)
	//	cout<<"Times:"<<tmp<<endl;
	ZZ_p x = power(num,tmp); //qpow(num,tmp); //power(num,tmp);
	//if(flag_print)
	//	cout<<"Res:"<<x<<endl<<"Check:"<<x * num<<endl;
	return x;
}

//return value: succesfully generated shares
	template <class FieldType>
int CompareGate<FieldType>::generateBitShares(int num)
{
	int suc_cnt=0;
	_bitShareOffset = 0;
	//we pick out 4*required bits
	int tot = 4 * eleSize * num;
	if(flag_print)
		cout<<"Required #bitshare:"<<tot<<endl;
	vector<FieldType> tempShares,resShares,secrets;
	_bitSharesValue.resize(tot);
	helper->getRandomShares(tot, tempShares);
	//if(flag_print)
	//	cout<<"Generating random shares, #"<<tot<<endl;
	vector<FieldType> PlainText;
	//if(flag_print)
	//	helper->openShare(tot,tempShares,PlainText);
	//now do multiplications and rule out 0 shares
	helper->DNMultVec(tempShares,tempShares,resShares,1);
	//if(flag_print)
	//	cout<<"DNMult finished"<<endl;
	helper->openShare(tot,resShares,secrets);
	if(flag_print)
	{
		cout<<"Share Opened in GenB(),"<<tot<<endl;
		cout<<secrets[0]<<","<<secrets[1]<<endl;
	}
	//for each opened, check if is 0, if not we generate the share
	for(int i=0; i<tot; i++)
		if(secrets[i] != FieldType(0)) //valid
		{
			//if(flag_print)
			//	cout<<"Trying to calc inv for "<<secrets[i]<<endl;
			FieldType cur = compSqrInverse(secrets[i]);
			//if(flag_print)
			//	cout<<"find invs:"<<cur<<","<<PlainText[i]<<endl;
			//cout<<"Inv2:"<<FieldType(1) / FieldType(2)<<endl;
			cur = ( (FieldType(1) / FieldType(cur)) * tempShares[i]+ FieldType(1))/ FieldType(2);
			//if(flag_print)
			//	cout<<"found share:"<<cur<<endl;
			/*if(flag_print)
			  {
			  vector<FieldType> tp(1),res;
			  tp[0] = cur;
			  helper->openShare(1,tp,res);
			  cout<<"Generated a "<<res[0]<<endl;
			  }*/
			_bitSharesValue[_bitShareOffset] = cur;
			suc_cnt++;
			_bitShareOffset++;
			if(flag_print)
			{
				vector<FieldType> d0(1),d1;
				d0[0] = cur;
				helper->openShare(1,d0,d1);
				//cout<<"generate a bit "<<d1[0]<<endl;
				if(d1[0]!=FieldType(0) && d1[0]!=FieldType(1))
				{
					cout<<"invalid generated bit!"<<endl;
					abort();
				}
			}
			//if(flag_print)
			//	cout<<"Suc count:"<<suc_cnt<<endl;
		}
	//else if(flag_print)
	//	cout<<i<<" failed"<<endl;
	if(flag_print)
		cout<<"Bit Generate done"<<endl;
	_bitShareOffset = 0;
	return suc_cnt;
}

	template<class FieldType>
void CompareGate<FieldType>::compRandom(vector<FieldType> &a, vector<FieldType> &b, vector<FieldType> &res)
{
	//here a,b,c means the *share* [a],[b],[c] that are hold by this party
	int cnt = a.size();
	if(flag_print)
		cout<<"Entering compRandom with size "<<cnt<<endl;
	vector<FieldType> c,w,x,y,tmp,tmp2,tmp3;
	for(int l=0; l<cnt; l++)
	{
		//auto gate = gates[l];
		//a.push_back(helper->gateShareArr[gate.input1]);
		//b.push_back(helper->gateShareArr[gate.input2]);
		c.push_back(a[l] - b[l]);
	}
	//FieldType c = a - b;
	//same notation as the paper, share of [a<p/2],[b<p/2],[a-b<p/2]
	//FieldType w,x,y,tmp,tmp2,tmp3;
	//if(flag_print)
	//	cout<<"Entering comphalf"<<endl;
	/*if(flag_print)
	  {
	  cout<<"checking input:"<<endl;
	  vector<FieldType> d0(3),d1;
	  d0[0] = a; d0[1] = b; d0[2]=c;
	  helper->openShare(3,d0,d1);
	  cout<<d1[0]<<","<<d1[1]<<","<<d1[2]<<endl;
	  }*/
	compHalfP(a,w);
	//	cout<<"Entering comphalf 2"<<endl;
	compHalfP(b,x);
	compHalfP(c,y);
	/*if(flag_print)
	  {
	  cout<<"checking compHalfP results:"<<endl;
	  vector<FieldType> d0(3),d1;
	  d0[0]=w; d0[1]=x; d0[2]=y;
	  helper->openShare(3,d0,d1);
	  cout<<d1[0]<<","<<d1[1]<<","<<d1[2]<<endl;
	  }*/
	//    cout<<"Exiting compHalf"<<endl;
	//two innovations since we need to compute xy and w*(x+y-2xy)
	//first round: compute xy
	helper->DNMultVec(x,y,tmp,1);
	for(int l=0; l<cnt; l++)
		tmp2.push_back(x[l] + y[l] - FieldType(2) * tmp[l]);
	helper->DNMultVec(w,tmp2,tmp3,1);
	if(res.size()<cnt)
		res.resize(cnt);
	for(int l=0; l<cnt; l++)
	{
		//auto gate = gates[l];
		//helper->gateShareArr[gate.output] = tmp3[l] + FieldType(1) - y[l] - x[l] + tmp[l];
		res[l] = tmp3[l] + FieldType(1) - y[l] - x[l] + tmp[l];
	}
}
// given known element c and binary share a, store (c<a) to dest
	template<class FieldType>
void CompareGate<FieldType>::compGiven(vector<vector<FieldType> > &a, vector<FieldType> &c, vector<FieldType> &dest)
{
}
	template<>
inline void CompareGate<ZZ_p>::compGiven(vector<vector<ZZ_p> > &a, vector<ZZ_p> &c, vector<ZZ_p> &dest)
{
	if(a.size()!=c.size())
	{
		cout<<"INVALID PARM for compGiven"<<endl;
		abort();
	}
	vector<ZZ_p> d0,d1;
	int tot = a.size();
	if(flag_print)
	{
		cout<<"Entering compGiven"<<endl;
	}
	vector<vector<ZZ_p> > cBits;
	cBits.resize(tot);
	for(int j=0; j<tot; j++)
	{
		ZZ tmp = rep(c[j]);
		cBits[j].resize(eleSize);
		for(int i=eleSize-1; i>=0; i--)
		{
			cBits[j][i] = tmp % 2;
			tmp /=2;
		}
		/*if(flag_print)
		  {
		  for(int i=0; i<eleSize; i++)
		  cout<<cBits[i]<<",";
		  cout<<endl;
		  for(int i=0; i<eleSize; i++)
		  cout<<d1[i]<<",";
		  cout<<endl;
		  }*/
	}
	if(flag_print)
		cout<<"compGiven::xorSum"<<endl;
	//step 1: compute bitwise xor
	vector<vector<ZZ_p> > xorSum;
	xorSum.resize(tot);
	for(int j=0; j<tot; j++)
	{
		xorSum[j].resize(eleSize);
		for(int i=0; i<eleSize; i++)
		{
			xorSum[j][i] = a[j][i] + cBits[j][i];
			if(cBits[j][i]==1)
				xorSum[j][i] -= 2 * a[j][i];
		}
	}
	/*if(flag_print)
	  {
	  cout<<"Xor sum:"<<endl;
	  helper->openShare(eleSize,xorSum,d1);
	  for(int i=0; i<eleSize; i++)
	  cout<<d1[i]<<",";
	  cout<<endl;
	  }*/
	//step 2: compute prefix or
	vector<vector<ZZ_p> > prefixOr,E;
	if(flag_print)
		cout<<"Entering PrefixOr"<<endl;
	compPrefixOr(xorSum,prefixOr);
	/*if(flag_print)
	  {
	  cout<<"Exiting compGiven::prefixor"<<endl;
	  helper->openShare(eleSize,prefixOr,d1);
	  for(int i=0; i<eleSize; i++)
	  cout<<d1[i]<<",";
	  cout<<endl;
	  }*/
	E.resize(tot);
	dest.resize(tot);
	for(int j=0; j<tot; j++)
	{
		E[j].resize(eleSize);
		E[j][0] = prefixOr[j][0];
		for(int i=1; i<eleSize; i++)
			E[j][i] = prefixOr[j][i] - prefixOr[j][i-1];
	}
	for(int j=0; j<tot; j++)
	{
		dest[j] = 0;
		for(int i=0; i<eleSize; i++)
			if(cBits[j][i]==0)
				dest[j] = dest[j] + E[j][i];
	}
	/*if(flag_print)
	  {
	  vector<ZZ_p> d0(1),d1;
	  d0[0]=dest;
	  helper->openShare(1,d0,d1);
	  if(d1[0]!=0 && d1[0]!=1)
	  {
	  cout<<"Invalid compGiven:"<<d1[0]<<endl;
	  abort();
	  }
	  cout<<"Exiting CompGiven"<<endl;
	  }*/
}
	template<class FieldType>
void CompareGate<FieldType>::compGivenNoPrefixOr(vector<vector<FieldType> > &a, vector<FieldType> &c, vector<FieldType> &dest)
{
}
	template<>
inline void CompareGate<ZZ_p>::compGivenNoPrefixOr(vector<vector<ZZ_p> > &x, vector<ZZ_p> &r, vector<ZZ_p> &dest)
{
	if(x.size()!=r.size())
	{
		cout<<"INVALID PARM for compGiven"<<endl;
		abort();
	}
	vector<ZZ_p> d0,d1;
	int tot = x.size();
	if(flag_print)
	{
		cout<<"Entering compGiven"<<endl;
	}

	vector<vector<ZZ_p> > rBits;
	rBits.resize(tot);
	for(int j=0; j<tot; j++)
	{
		ZZ tmp = rep(r[j]);
		rBits[j].resize(eleSize);
		for(int i=eleSize-1; i>=0; i--)
		{
			rBits[j][i] = tmp % 2;
			tmp /=2;
		}
	}
	if(flag_print)
		cout<<"compGivenNew:step 1"<<endl;
	//step 1: compute b and c
	vector<vector<ZZ_p> > b,c;
	b.resize(tot);
	c.resize(tot);
	for(int j=0; j<tot; j++)
	{
		b[j].resize(eleSize);
		c[j].resize(eleSize);
		for(int i=0; i<eleSize; i++)
		{
			b[j][i] = 1 - rBits[j][i] - x[j][i] + 2 * rBits[j][i] * x[j][i];
			if(rBits[j][i] == 0) c[j][i] = x[j][i];
			else c[j][i] = 0;
		}
	}
	//step 2: compute r
	computeAnswer(b,c,dest);
}
	template<class FieldType>
void CompareGate<FieldType>::computeAnswer(vector<vector<FieldType> > &b,vector<vector<FieldType> > &c, vector<FieldType> &res)
	//{
	//}
	//tempate<>
	//inline void CompareGate<ZZ_p>::computeAnswer(vector<vector<ZZ_p> > &b,vector<vector<ZZ_p> > &c, vector<ZZ_p> &res)
{
	//step 1: compute C and B (by dividing half)
	int tot = b.size();
	if(b.size()!=c.size())
	{
		cout<<"Invalid Param in compAnswer"<<endl;
		abort();
	}
	if(flag_print)
		cout<<"Entering compAnswer with "<<b.size()<<endl;
	if(res.size()<b.size())
		res.resize(b.size());
	vector<FieldType> veca,vecb,vecres,tmpans;
	int _cnt=0;
	//all length should be even, since we have length 32 here
	for(int l=0; l<tot; l++)
	{
		if(flag_print)
			cout<<"No "<<l<<", length:"<<b[l].size()<<endl;
		if(b[l].size()==1) //the answer is c_1
			res[l] = c[l][0];
		else if(b[l].size()>1)
		{
			for(int i=0; i<b[l].size(); i+=2)
			{
				veca.push_back(b[l][i]);
				vecb.push_back(c[l][i+1]);
				veca.push_back(b[l][i]);
				vecb.push_back(b[l][i+1]);
			}
		}
		else //size = 0 ?
		{
			cout<<"0-sized array detected!"<<endl;
			abort();
		}
	}
	helper->DNMultVec(veca,vecb,vecres,1);
	//restore back to get B and C
	vector<vector<ZZ_p> > B,C;
	for(int l=0; l<tot; l++)
		if(b[l].size()>1)
		{
			B.resize(B.size()+1);
			C.resize(C.size()+1);
			for(int i=0; i<b[l].size(); i+=2)
			{
				B[B.size()-1].push_back(vecres[_cnt+1]);
				C[C.size()-1].push_back(vecres[_cnt] + c[l][i]);
				_cnt+=2;
			}
		}
	//call sub to get answer
	if(B.size()>0)
		computeAnswer(B,C,tmpans);
	//restore answer
	for(int l=0,_t=0; l<tot; l++)
		if(b[l].size()>1)
		{
			res[l] = tmpans[_t];
			_t++;
		}
}
	template<class FieldType>
void CompareGate<FieldType>::compP(vector<vector<FieldType> > &a, vector<FieldType> &dest)
{
}
	template<>
inline void CompareGate<ZZ_p>::compP(vector<vector<ZZ_p> > &a, vector<ZZ_p> &dest)
{
	if(flag_print)
		cout<<"Entering compP"<<endl;
	vector<ZZ_p> cBits;
	cBits.resize(eleSize);
	ZZ tmp((1ll<<61)-1);
	for(int i=eleSize-1; i>=0; i--)
	{
		cBits[i] = tmp % 2;
		tmp /=2;
	}
	//step 1: compute bitwise xor
	int _cnt = a.size();
	if(dest.size() < _cnt)
		dest.resize(_cnt);
	vector<vector<ZZ_p> > xorSum(_cnt);
	for(int l=0; l<_cnt; l++)
	{
		for(int i=0; i<eleSize; i++)
		{
			xorSum[l].resize(eleSize);
			xorSum[l][i] = a[l][i] + cBits[i];
			if(cBits[i]==1)
				xorSum[l][i] -= 2 * a[l][i];
		}
	}
	//step 2: compute prefix or
	vector<vector<ZZ_p> > prefixOr;
	vector<ZZ_p> E;
	if(flag_print)
		cout<<"Entering PrefixOr"<<endl;
	compPrefixOr(xorSum,prefixOr);
	if(flag_print)
		cout<<"Exiting compP::prefixor"<<endl;
	/*if(flag_print)
	  {
	  cout<<"Checking validation of prefixor"<<endl;
	  vector<ZZ_p> tri;
	  helper->openShare(1,prefixOr[0],tri);
	  cout<<"Passed"<<endl;
	  }*/
	E.resize(eleSize);
	for(int l=0; l<_cnt; l++)
	{
		E[0] = prefixOr[l][0];
		for(int i=1; i<eleSize; i++)
			E[i] = prefixOr[l][i] - prefixOr[l][i-1];
		dest[l] = 0;
		for(int i=0; i<eleSize; i++)
			if(cBits[i]==0)
				dest[l] = dest[l] + E[i];
	}
	/*if(flag_print)
	  {
	  cout<<"Checking dest"<<endl;
	  vector<ZZ_p> t1(1),t2(1);
	  t1[0] = dest;
	  helper->openShare(1,t1,t2);
	  if(t2[0]!=0 && t2[0]!=1)
	  {
	  cout<<"Invalid Comp:"<<t2[0]<<endl;
	  abort();
	  }
	  cout<<"Exiting compP"<<endl;
	  }*/
	//step 1: compute bitwise xor
}
//compare bit sharing (a<b)
// bits are given from MSB to LSB

// compare the function (x < p/2)
	template<class FieldType>
void CompareGate<FieldType>::compHalfP(vector<FieldType> &x,vector<FieldType> &dest)
{
}

//there are several pairs to be compared [ ? < p/2]
	template<>
inline void CompareGate<ZZ_p>::compHalfP(vector<ZZ_p> &X, vector<ZZ_p> &dest)
{
	if(flag_print)
		cout<<"Now at compHalfP"<<endl;
	//we want to compute ([2x]_0)
	int tot = X.size();
	//step 1: generate [r]_B and [r]_p
	vector<ZZ_p> r;
	vector<vector<ZZ_p> > rBits;
	vector<ZZ_p> d0,d1;
	if(flag_print)
		cout<<"Getting random share r"<<endl;
	getRandomBitShare(tot,r,rBits);
	if(flag_print)
		cout<<"compHalfP::Successfully get shares"<<endl;
	if(flag_print)
	{
		cout<<"checking valid of bit share"<<endl;
		vector<ZZ_p> tmp;
		vector<ZZ> tmpres;
		for(int i=0; i<tot; i++)
		{
			helper->openShare(eleSize,rBits[i],tmp);
			ZZ _t(0);
			for(int j=0; j<eleSize; j++)
			{
				cout<<"bit "<<j<<":"<<tmp[j]<<endl;
				_t = _t * 2 + conv<ZZ>(tmp[j]);
			}
			if(_t > (1ll<<61))
			{
				cout<<"too large generated!"<<endl;
				abort();
			}
			tmpres.push_back(_t);
		}
		helper->openShare(tot,r,tmp);
		for(int i=0; i<tot; i++)
			if((tmp[i]-conv<ZZ_p>(tmpres[i])) != ZZ_p(0))
			{
				cout<<"find difference:"<<tmp[i]<<","<<tmpres[i]<<endl;
			}
	}
	//step 2: compute [c] = [x] + [r], and reveal c
	vector<ZZ_p> resvec,tmpvec;
	for(int i=0; i<tot; i++)
	{
		ZZ_p res = 2 * X[i] + r[i];
		resvec.push_back(res);
		/*ZZ_p res = x + r[0], c;
		  vector<ZZ_p> resvec(1),tmpvec(1);
		//resvec.resize(1); 
		resvec[0]=res;*/
	}
	helper->openShare(tot,resvec,tmpvec);
	if(flag_print)
	{
		cout<<"compHalfP::Successfully opened shares"<<endl;
		cout<<"value of c:"<<tmpvec[0]<<endl;
	}
	//c = tmpvec[0];
	//step 3: compute [c <_b r]
	vector<ZZ_p> comp_res;
	//compGiven(rBits,tmpvec,comp_res);
	compGivenNoPrefixOr(rBits,tmpvec,comp_res);
	if(flag_print)
		cout<<"compGiven finished"<<endl;
	if(comp_res.size()!=tot)
	{
		cout<<"Invalid compGiven length"<<endl;
		abort();
	}
	//step 4: compute the overall result
	vector<ZZ_p> rescr0;
	for(int i=0; i<tot; i++)
	{
		if(rep(tmpvec[i]) % 2 ==0)
			rescr0.push_back(rBits[i][this->eleSize-1]);
		else rescr0.push_back(1 - rBits[i][this->eleSize-1]);
	}
	vector<ZZ_p> tmpres;
	helper->DNMultVec(rescr0,comp_res,tmpres,1);
	if(flag_print)
		cout<<"compHalfP::Single Mult finished"<<endl;

	if(dest.size()<tot)
		dest.resize(tot);
	for(int i=0; i<tot; i++)
	{
		//ZZ_p tmpr = 2 * tmpres[i];
		dest[i] = 1 - (comp_res[i] + rescr0[i] - 2 * tmpres[i]);
		//dest = 1 - dest;
	}
	if(flag_print)
	{
		cout<<"compHalfP, result checking:"<<endl;
		vector<ZZ_p> d1;
		helper->openShare(tot,dest,d1);
		cout<<d1[0]<<endl;
		if(d1[0] != 0 && d1[0] != 1)
		{
			cout<<"Invalid CompHalfP Result"<<endl;
			abort();
		}
	}

}
	template<class FieldType>
void CompareGate<FieldType>::getRandomBitShare(int num,vector<FieldType> &res,vector<vector<FieldType> > &bits)
{
	if(flag_print)
		cout<<"Now at GetRandomBit with "<<num<<endl;
	//for safety we generate 2 * num
	bits.resize(num);
	vector<ZZ_p> tmp;
	tmp.resize(num);
	res.resize(num);
	int suc;
	for(suc = 0; suc < num; suc++)
	{
		bits[suc].resize(eleSize);
		res[suc]=0;
		//generate from MSB to LSB
		//first 3 bits: 0
		if(_bitSharesValue.size() <= _bitShareOffset + eleSize)
		{
			cout<<"Not enough bit share!"<<endl;
			generateBitShares(num+5-suc);
		}
		for(int j=0; j<3; j++)
		{
			bits[suc][j] = 0;
			//bits[suc][j] = _zeroShares[_zeroShareOffset++];
			//res[suc] = res[suc] * FieldType(2) + bits[suc][j];
			if(_zeroShareOffset >= _zeroShares.size())
			{
				cout<<"not enough zero shares"<<endl;
				abort();
			}
		}
		for(int j=3; j<eleSize && _bitShareOffset < _bitSharesValue.size(); j++)
		{
			bits[suc][j] = _bitSharesValue[_bitShareOffset++];
			res[suc] = res[suc] * FieldType(2) + bits[suc][j];
		}
		if(_bitSharesValue.size() <= _bitShareOffset) //not enough?
		{
			cout<<"This can't happen!"<<endl;
			abort();
			if(flag_print)
				cout<<"Not enough bit share! Regenerating....."<<endl;
			generateBitShares(num+5-suc);
		}
	}
	if(flag_print)
	{
		cout<<"gen bit done, checking"<<endl;
		vector<FieldType> tmp;
		helper->openShare(num,res,tmp);
		for(int i=0; i<num; i++)
			cout<<"generated a number "<<tmp[i]<<endl;
		for(int i=0; i<num; i++)
		{
			vector<FieldType> tmp2;
			helper->openShare(eleSize,bits[i],tmp2);
			FieldType cur(0);
			for(int j=0; j<eleSize; j++)
				cur = cur * FieldType(2) + tmp2[j];
			if(cur != tmp[i])
				cout<<"difference in share:"<<cur<<","<<tmp[i]<<endl;
		}
	}
	

	/*compP(bits,chk_suc);
	  if(flag_print)
	  cout<<"genBit:compP ended"<<endl;
	  if(flag_print)
	  {
	  cout<<"Checking consistency of bit share"<<endl;
	  vector<FieldType> r1(1),r2;
	  r1[0]=res[suc];
	  helper->openShare(1,r1,r2);
	  cout<<"Check passed"<<endl;
	  }
	//TODO: should we consider multi-thread here?
	//no need to compare ??
	//IF NO NEED TO COMPARE, JUST COMMENT THE FOLLOWING LINES BY /*
	//if(flag_print)
	//	cout<<"Exiting CompP"<<endl;
	vector<FieldType> vc2;
	helper->openShare(num,chk_suc,vc2);*/
	//TODO: can we use '=' like this?
	//notice: compGiven gives [p < c], here we need to change direction
	/*
	int tcnt=0;
	for(int l=0; l<num && tcnt<num; l++)
		//if(vc2[l]==FieldType(0)) //valid share 
	{
		//copy
		for(int j=0; j<eleSize; j++)
			bits[tcnt][j] = bits[l][j];
		res[tcnt] = res[l];
		tcnt++;
	}
	if(tcnt<num) 
	{
		if(flag_print)
			cout<<"First try failed!"<<endl;
		vector<vector<FieldType> > _bits;
		vector<FieldType> _res;
		getRandomBitShare(num - tcnt, _res, _bits);
		for(int l=0,pos=tcnt; l<num-tcnt; l++,pos++)
		{
			res[pos] = _res[l];
			for(int j=0; j<eleSize; j++)
				bits[pos][j] = _bits[l][j];
		}
	}
	else if(flag_print)
		cout<<"First try passed!"<<endl;
*/
	if(flag_print)
	{
		vector<FieldType> tmp;
		helper->openShare(eleSize,bits[0],tmp);
		cout<<"Check passed"<<endl;
	}
	//res.resize(num);
	//bits.resize(num);
}

	template<class FieldType>
void CompareGate<FieldType>::compPrefixOr(vector<vector<FieldType> > &a, vector<vector<FieldType> > &res)
{
	//for convenience we make a copy, as stated in the paper 
	int _cnt = a.size();
	vector<vector<FieldType> > A;
	int bas = 0;
	int bastot = 0;
	int bassiz = 0;
	vector<int> tot,groupSize,groupNum;
	for(int l=0; l<_cnt; l++)
	{
		tot.push_back(a[l].size());
		groupSize.push_back(int(sqrt(tot[l]))+1);
		//if(flag_print)
		//	cout<<"Entering PrefixOr with length "<<tot<<","<<groupSize<<endl;
		//except the last group, each group has groupSize elements
		groupNum.push_back((tot[l]-1) / groupSize[l] + 1);
		A.resize(A.size() + groupNum[l]);
		//vector<vector<FieldType> > A;

		for(int i=0,tmp=0; i<groupNum[l]; i++)
		{
			A[i+bas].resize(groupSize[l]);
			for(int j=0; j<groupSize[l] && tmp < tot[l]; j++,tmp++)
			{
				A[i+bas][j] = a[l][tmp];
				if(tmp == tot[l])
					A[i+bas].resize(j+1);
			}
		}
		bas += groupNum[l];
	}
	//step 1: compute blockwise fanin-or
	//each group has lambda+1 elements, and there is lambda groups
	vector<FieldType> X,Y,F;
	if(flag_print)
		cout<<"Entering FaninOr"<<endl;
	compFanInOr(bas,A,X);
	/*if(flag_print)
	  {
	  cout<<"Checking 1st FanInOr"<<endl;
	  vector<FieldType> t0;
	  helper->openShare(groupNum,X,t0);
	  for(int i=0; i<groupNum; i++)
	  cout<<t0[i]<<",";
	  cout<<"Passed"<<endl;
	  }*/
	vector<vector<FieldType> > tmpX;
	bas = 0;
	for(int l=0; l<_cnt; l++)
	{
		tmpX.resize(tmpX.size() + groupNum[l]);
		for(int i=0; i<groupNum[l]; i++)
		{
			tmpX[bas + i].resize(i+1);
			for(int j=0; j<=i; j++)
				tmpX[bas + i][j] = X[j];
		}
		bas += groupNum[l];
	}
	if(flag_print)
		cout<<"Entering 2nd FaninOr"<<endl;
	compFanInOr(bas,tmpX,Y);
	/*if(flag_print)
	  {
	  cout<<"Checking 2nd FaninOr"<<endl;
	  vector<FieldType> t0;
	  helper->openShare(groupNum,Y,t0);
	  for(int i=0; i<groupNum; i++)
	  cout<<t0[i]<<",";
	  cout<<"Passed"<<endl;
	  }*/
	bas = 0;
	for(int l=0; l<_cnt; l++)
	{
		F.resize(F.size() + groupNum[l]);
		F[bas + 0] = X[bas + 0];
		for(int i=1; i<groupNum[l]; i++)
			F[bas + i] = Y[bas + i] - Y[bas + i-1];
		bas += groupNum[l];
	}

	bas = 0;
	bastot = 0;
	vector<FieldType> tmpVec1,tmpVec2,G,C,B,S;
	for(int l=0; l<_cnt; l++)
	{
		tmpVec1.resize(tmpVec1.size() + tot[l]);
		tmpVec2.resize(tmpVec2.size() + tot[l]);
		for(int i=0,cnt=0; i<groupNum[l]; i++)
			for(int j=0; j<groupSize[l] && cnt<tot[l]; j++,cnt++)
			{
				tmpVec1[bastot + cnt] = F[bas + i];
				tmpVec2[bastot + cnt] = A[bas + i][j];
			}
		bastot += tot[l];
		bas += groupNum[l];
	}
	if(flag_print)
		cout<<"Entering helper->DNMultVec after 2nd FaninOr"<<endl;
	helper->DNMultVec(tmpVec1,tmpVec2,G,1);
	bas = 0;
	bastot = 0;
	bassiz = 0;
	for(int l=0; l<_cnt; l++)
	{
		C.resize(C.size() + groupSize[l]);
		for(int i=0; i<groupSize[l]; i++)
		{
			C[i+bassiz]=0;
			for(int j=0; j<groupNum[l]; j++)
			{
				//calculate the correct label
				int pos = j * groupSize[l] + i;
				if(pos < tot[l])
					C[i+bassiz] += G[pos+bastot];
			}
		}
		bastot += tot[l];
		bassiz += groupSize[l];
	}

	bassiz = 0;
	tmpX.resize(0);
	for(int l=0; l<_cnt; l++)
	{
		tmpX.resize(tmpX.size() + groupSize[l]);
		for(int i=0; i<groupSize[l]; i++)
		{
			tmpX[i+bassiz].resize(i+1);
			for(int j=0; j<=i; j++)
				tmpX[i+bassiz][j] = C[j];
		}
		bassiz += groupSize[l];
	}
	if(flag_print)
		cout<<"Entering 3rd FanInOr"<<endl;
	compFanInOr(bassiz,tmpX,B);
	/*if(flag_print)
	  {
	  cout<<"Checking 3rd FaninOr"<<endl;
	  vector<FieldType> t0;
	  helper->openShare(groupSize,B,t0);
	  for(int i=0; i<groupSize; i++)
	  cout<<t0[i]<<",";
	  cout<<"Passed"<<endl;
	  }*/
	bastot = 0;
	bassiz = 0;
	for(int l=0; l<_cnt; l++)
	{
		for(int i=0; i<groupNum[l]; i++)
			for(int j=0; j<groupSize[l]; j++)
			{
				int pos = i*groupSize[l] + j;
				if(pos<tot[l])
					tmpVec2[pos + bastot] = B[j + bassiz];
			}
		bastot += tot[l];
		bassiz += groupSize[l];
	}
	helper->DNMultVec(tmpVec1,tmpVec2,S,1);

	bastot = 0;
	bas = 0;
	if(res.size()<_cnt)
		res.resize(_cnt);
	for(int l=0; l<_cnt; l++)
	{
		res[l].resize(tot[l]);
		for(int i=0; i<groupNum[l]; i++)
			for(int j=0; j<groupSize[l]; j++)
			{
				int pos = i*groupSize[l] + j;
				if(pos < tot[l])
					res[l][pos] = S[pos + bastot] + Y[i + bas] - F[i + bas];
			}
		bas += groupNum[l];
		bastot += tot[l];
	}
}
	template<class FieldType>
void CompareGate<FieldType>::compFanInOr(int num,vector<vector<FieldType> >&a, vector<FieldType> &res)
{
	if(num!=a.size())
	{
		cout<<"Invalid parameters for COMPFANINOR"<<endl;
		abort();
	}
	//notice all steps should be done in parallel
	//step 1: compute A
	vector<FieldType> A,Avec;
	A.resize(num);
	int tot = 0;
	for(int i=0; i<num; i++)
	{
		tot += a[i].size();
		A[i]=1;
		for(int j=0; j<a[i].size(); j++)
			A[i] = A[i] + a[i][j];
	}
	vector<FieldType> debug0,debug1;
	/*if(flag_print)
	  {
	  cout<<"revealing a"<<endl;
	  for(int i=0; i<num; i++)
	  {
	  cout<<i<<":"<<endl;
	  helper->openShare(a[i].size(),a[i],debug0);
	  for(int j=0; j<a[i].size(); j++)
	  cout<<debug0[j]<<",";
	  cout<<endl;
	  }
	  cout<<"checking A"<<endl;
	  helper->openShare(num,A,debug0);
	  for(int i=0; i<num; i++)
	  cout<<debug0[i]<<",";
	  cout<<endl;
	  }*/

	//step 2: generate pairs of inversion
	vector<FieldType> b,invb,tmpc,C,Cres;
	if(flag_print)
		cout<<"Getting Inv Pairs"<<endl;
	getRandomInvPairs(tot,b,invb,tmpc);
	if(flag_print)
		cout<<"Inv Pair generated,"<<tot<<","<<b.size()<<","<<invb.size()<<","<<tmpc.size()<<endl;
	Avec.resize(tot);
	//step 3: finish multiplications on C
	for(int i=0,cnt=0; i<num; i++)
	{
		//if(flag_print)
		//	cout<<i<<","<<cnt<<","<<a[i].size()<<endl;
		if(i!=0) //tiny modification on array tmpc
			tmpc[cnt] = tmpc[cnt] * invb[cnt-1];
		for(int j=0; j<a[i].size(); j++,cnt++)
			Avec[cnt] = A[i];
	}
	//if(flag_print)
	//	cout<<"helper->DNMultVec"<<endl;
	helper->DNMultVec(Avec,tmpc,C,1);
	//if(flag_print)
	//	cout<<"fanInOr:helper->openShare(C):"<<C.size()<<endl;
	helper->openShare(tot,C,Cres);
	//if(flag_print)
	//	cout<<"Share(C) opened"<<endl;
	//step 4: retrieve the polynomials and do calculations
	res.resize(num);
	for(int i=0,cnt=0; i<num; i++)
	{
		if(OrVector[a[i].size()].size()==0) //not initied ?
		{
			if(flag_print)
				cout<<"middle Intering "<<a[i].size()<<endl;
			vector<FieldType> tmp,tmpnum; //(a[i].size()+1),tmpnum(a[i].size()+1);
			tmp.resize(a[i].size()+1);
			tmpnum.resize(a[i].size()+1);

			tmp[0] = FieldType(0);
			tmpnum[0] = FieldType(1);
			//if(flag_print)
			//  cout<<"tot:"<<tot<<" Intering 0"<<endl;
			//interp.interpolate(tmpnum, tmp, OrVector[]);
			for(int j=1; j<=a[i].size(); j++)
			{

				tmp[j]=FieldType(1);
				tmpnum[j]=FieldType(j+1);
			}
			//if(flag_print)
			//	cout<<"Ready for Interpolate"<<endl;
			helper->interp.interpolate(tmpnum, tmp, OrVector[a[i].size()]);
			/*if(flag_print)
			  {
			  for(int j=0; j<a[i].size()+1; j++)
			  cout<<OrVector[a[i].size()][j]<<","<<endl;
			  }
			  if(flag_print)
			  {
			//interp.printPolynomial(OrVector[a[i].size()]);
			cout<<"Interpolate Done"<<endl;
			}*/
		}
		//if(flag_print)
		//	cout<<a[i].size()<<","<<C.size()<<","<<OrVector[a[i].size()].size()<<","<<b.size()<<endl;
		res[i] = OrVector[a[i].size()][0];
		FieldType Cprod = FieldType(1);
		for(int j=0; j<a[i].size(); j++,cnt++)
		{
			Cprod = Cprod * Cres[cnt];
			/*if(flag_print)
			  {
			  debug0.resize(1);
			  debug0[0]=Cprod * b[cnt];
			  helper->openShare(1,debug0,debug1);
			  cout<<"value of A^"<<j+1<<":"<< debug1[0]<<endl;
			  }*/
			res[i] = res[i] + OrVector[a[i].size()][j+1] * Cprod * b[cnt];
		}
		/*if(flag_print)
		  {
		  cout<<"Checking res array "<<i<<endl;
		  vector<FieldType> t0(1),t1;
		  t0[0]=res[i];
		  helper->openShare(1,t0,t1);
		  if(t1[0]!=0 && t1[0]!=1)
		  {
		  cout<<"Invalid FanInOr value"<<endl;
		  abort();
		  }
		  cout<<"Passed, or value:"<<t1[0]<<endl;
		  }*/
	}
	if(flag_print)
		cout<<"Gen Inv finished"<<endl;
}
	template<class FieldType>
void CompareGate<FieldType>::getRandomInvPairs(int num, vector<FieldType> &b, vector<FieldType> &invb, vector<FieldType> &ctmp)
{
	if(flag_print)
		cout<<"Entering genInv:"<<num;
	vector<FieldType> bp,B,Bres;
	bp.resize(2*num-1);
	b.resize(2*num-1);
	if(flag_print)
		cout<<b.size()<<" "<<bp.size();
	helper->getRandomShares(num,b);
	helper->getRandomShares(num,bp);
	if(flag_print)
		cout<<b.size()<<" "<<bp.size();
	for(int i=num; i<2*num-1; i++)
	{
		b[i] = b[i-num];
		bp[i] = bp[i-num+1];
	}
	if(flag_print)
		cout<<"genInv:DNMult,size:"<<bp.size()<<","<<b.size()<<endl;
	helper->DNMultVec(b,bp,B,1);
	if(flag_print)
		cout<<"genInv:helper->openShare,size:"<<B.size()<<endl;
	helper->openShare(2*num-1,B,Bres);
	if(flag_print)
		cout<<"genInv:DNMult finished"<<endl;
	invb.resize(num);
	ctmp.resize(num);
	for(int i=0; i<num; i++)
		invb[i] = FieldType(1) / Bres[i] * bp[i];
	vector<FieldType> tr,tr2;
	if(flag_print)
	{
		cout<<"genInv:validating"<<endl;
		helper->openShare(num,invb,tr);
		helper->openShare(num,b,tr2);
		for(int i=0; i<num; i++)
			if(tr[i] * tr2[i] != FieldType(1))
				cout<<"Not inv for "<<tr[i]<<" and "<<tr2[i]<<endl;
		cout<<"genInv:validat ended"<<endl;
	}
	ctmp[0] = invb[0];
	//b_i * b_{i+1}^{-1} = b_i * B_{i+1}^{-1} * b'_{i+1}
	for(int i=1; i<num; i++)
		ctmp[i] = FieldType(1) / Bres[i] * Bres[i+num-1];
	if(flag_print)
	{
		vector<ZZ_p> ttr;
		cout<<"cTmp:validating"<<endl;
		helper->openShare(num,ctmp,ttr);
		for(int i=1; i<num; i++)
			if(ttr[i] != tr2[i-1] * tr[i])
				cout<<"Not CTEMP for "<<tr2[i-1]<<" and "<<tr[i]<<endl;
		cout<<"passed"<<endl;
	}
}

	template<class FieldType>
void CompareGate<FieldType>::getBitDecomp(FieldType &a, vector<FieldType> &dest)
{
}

	template<>
inline void CompareGate<ZZ_p>::getBitDecomp(ZZ_p &a, vector<ZZ_p> &dest)
{
	/*	//step 1: generate bit
		vector<ZZ_p> rBits,tmp,Ctmp,Cres;
		ZZ_p p = _ZZPRIME;
		ZZ_p r,c;
		getRandomBitShare(1,r,rBits);
		if(dest.size()<eleSize)
		dest.resize(eleSize);
	//step 2: reveal c
	c = a-r;
	Ctmp[0]=c;
	helper->openShare(1,Ctmp,Cres);
	c = Cres[0];
	if(c==0)
	{
	//dest = r, by coincidence
	for(int i=0; i<eleSize; i++)
	dest[i] = rBits[i];
	return;
	}
	//step 3: compute q
	ZZ_p t0 = p - c - 1,q;
	compGiven(rBits,t0,q);
	//now q is the share of compGiven result
	//step 4: restore bit decomp of g = 2^l + c - qp
	vector<ZZ_p> f1Bits,f2Bits,gBits;*/
}

template<class FieldType>
void CompareGate<FieldType>::TruncPR(vector<FieldType> &a,vector<FieldType> &res)
{
}

template<>
inline void CompareGate<ZZ_p>::TruncPR(vector<ZZ_p> &a,vector<ZZ_p> &res)
{
	
	int tot = a.size();
	//step 1: generate r and r'
	int m = _m; // # of decimal digits
	int k = _k; // # of range
	int kappa = _kappa;
	ZZ_p khalf(1ll<<(k-1));
	ZZ_p twom(1ll<<_m);
	vector<ZZ_p> r,r1,r2,r2open;
	vector<vector<ZZ_p> > rBits;
	if(flag_print)
	{
		cout<<"TruncPR: checking valid"<<endl;
		helper->openShare(tot,a,r);
		cout<<"passed"<<endl;
	}
	getRandomBitShare(tot,r,rBits);
	if(flag_print)
		cout<<"TruncPR: get bit share done"<<endl;
	//step 2: add everything up and open them
	//notice the overall process will [not] exceed 2^31
	for(int i=0; i<tot; i++)
	{
		//first use the lower m digits to restore r1
		r1.push_back(ZZ_p(0));
		for(int j=eleSize - m; j < eleSize; j++)
			r1[i] = r1[i] * 2 + rBits[i][j];
		//now restore r2
		r2.push_back(ZZ_p(0));
		for(int j = eleSize - k - kappa; j < eleSize; j++)
			r2[i] = r2[i] * 2 + rBits[i][j];
		r2[i] = r2[i] + a[i] + khalf;
		//notice: r2 should be smaller than about 2^24, overwhelmed by a[i]
	}
	helper->openShare(tot,r2,r2open);
	if(flag_print)
		cout<<"TruncPR: r2 passed"<<endl;
	if(res.size()<tot)
		res.resize(tot);
	for(int i=0; i<tot; i++)
	{
		r2open[i] = ZZ_p(rep(r2open[i]) % (1ll<<m));
		res[i] = (a[i] - r2open[i] + r1[i]) / twom;
	}
	/*
	int tot = a.size();
	if(res.size()<tot)
		res.resize(tot);
	for(int i=0; i<tot; i++)
		res[i]=a[i];
		*/
}
template<class FieldType> 
void CompareGate<FieldType>::doubleVecMult(vector<FieldType> &a,vector<FieldType> &b,vector<FieldType> &res)
{
	vector<FieldType> tmp;
	//step 1: normal product
	helper->DNMultVec(a,b,tmp,1);
	//step 2: do trunc
	TruncPR(tmp,res);
}

template<class FieldType>
void CompareGate<FieldType>::SoftThres(vector<FieldType> &thres, vector<FieldType> &a, vector<FieldType> &res)
{
}

template<>
inline void CompareGate<ZZ_p>::SoftThres(vector<ZZ_p> &thres,vector<ZZ_p> &a,vector<ZZ_p> &res)
{
	//thres is [open value]
	int tot = a.size();
	vector<ZZ_p> b,thresPos,thresNeg;	
	ZZ_p pad(1ll<<(_k-1));
	//thres = thres + (1ll << (_k-1)); //automaically moded by q
	for(int i=0; i<tot; i++)
	{
		b.push_back(a[i]+pad);
		thresPos.push_back(thres[i] + pad);
		thresNeg.push_back(pad - thres[i]);
	}
	vector<ZZ_p> res1,res2,tmp1,tmp2,tmp3,add1,add2,add3;
	//test 1: a[i] ?> thres
	//compGiven(a,thresPos,res1);
	//compGiven(a,thresNeg,res2);
	//for(int i=0; i<dim; i++)
	//	res2[i] = 1-res2[i];
	compRandom(b, thresPos,res1);
	compRandom(thresNeg, b,res2);
	if(flag_print)
	{
		cout<<"SoftThres: checking result of comprandom"<<endl;
		vector<ZZ_p> _t0,_t1;
		helper->openShare(res1.size(),res1,_t0);
		helper->openShare(res2.size(),res2,_t1);
		cout<<"comprandom passed"<<endl;
		for(int i=0; i<tot; i++)
			cout<<"result of "<<i<<":"<<_t0[i]<<","<<_t1[i]<<endl;
	}
	//notice res1 and res2 are 0 ot 1 (without 2^-k)
	helper->DNMultVec(res1,res2,tmp3,1);
	if(flag_print)
		cout<<"end of mult1"<<endl;
	helper->DNMultVec(tmp3, a, add3, 1); 
	if(flag_print)
		cout<<"end of mult2"<<endl;
	for(int i=0; i<tot; i++)
	{
		tmp1.push_back(thres[i] * (1-res1[i]));
		tmp2.push_back(thres[i] * (1-res2[i]));
	}
	//TruncPR(tmp1,add1);
	//TruncPR(tmp2,add2);
	if(res.size()<tot)
		res.resize(tot);
	for(int i=0; i<tot; i++)
		res[i]= a[i] + tmp2[i] - add3[i] -tmp1[i];
}
template<class FieldType>
void CompareGate<FieldType>::doubleInverse(FieldType a,FieldType &res)
{
}

template<>
inline void CompareGate<ZZ_p>::doubleInverse(ZZ_p a,ZZ_p &res)
{
	ZZ bas(1ll<<_m);
	ZZ r = bas / rep(a);
	conv(res,r); //res = r;
	if(flag_print)
	{
		cout<<bas<<","<<a<<","<<res<<endl;
	}
}

	template<class FieldType>
void CompareGate<FieldType>::runLasso(int iter,FieldType lambda, FieldType rho, vector<vector<FieldType> > & Ai, vector<FieldType> &bi, vector<FieldType> &res)
{
}
template<>
inline void CompareGate<ZZ_p>::runLasso(int iter,ZZ_p lambda, ZZ_p rho, vector<vector<ZZ_p> > & Ai, vector<ZZ_p> &bi, vector<ZZ_p> &res)
{
	int dim = Ai.size(); //Ai: dim * dim matrix, bi: dim * 1 vector
	int N = helper->N;
	int T = helper->T;
	int fieldByteSize = eleSize / 8; //field->getElementSizeInBytes();	
	ZZ_p invN,invrho;
	ZZ_p ZN(N);
	doubleInverse(ZN,invN);
	doubleInverse(rho,invrho);
	if(flag_print)
		cout<<"InvN:"<<N<<","<<invN<<endl;

	
	//each party are holding #dim shares, shareOfW[i][j] means the share of the w_i[j]
	//shareofAi[i][j] means the share of Ai[j] 
	vector<vector<ZZ_p> > shareOfW(N),shareOfU(N),shareOfB(N);
	vector<vector<vector<ZZ_p> > > shareOfA(N);
	vector<ZZ_p> shareOfZ;
	// -- prepare the shares for the input
	// also need to prepare shares for P and P/2
	int index = 0;
	//vector<int> sizes(N);
	//before start: send shares of Ai and bi to other parties.
	//Since Ai and bi are plaintext, we need to generate corresponding shares
	vector<ZZ_p> x1(N), y1(N);
	vector<vector<ZZ_p>> sendBufsElements(N);
	vector<vector<byte>> sendBufsBytes(N);
	vector<vector<byte>> recBufBytes(N);
	vector<vector<ZZ_p>> recBufElements(N); 
	vector<int> sizes(N);
	//step 1: send and receive shares of Ai, bi, w, z and u 

	// Ai
	if(flag_print)
		cout<<"Lasso:step 1:"<<endl;
	for(int l1=0; l1<dim; l1++)
		for(int l2=0; l2<dim; l2++)
		{
			// get the expected sizes from the other parties

			// the value of a_0 is the input of the party.
			x1[0] = Ai[l1][l2]; //field->GetElement(Ai[l1][l2]);

			// generate random degree-T polynomial
			for (int i = 1; i < T + 1; i++) {
				// A random field element, uniform distribution
				x1[i] = field->Random();
			}

			helper->matrix_vand.MatrixMult(x1, y1, T + 1);

			// prepare shares to be sent
			for (int i = 0; i < N; i++) {
				sendBufsElements[i].push_back(y1[i]);
			}
		} 
	//bi
	for(int l1=0; l1<dim; l1++)
	{
		// get the expected sizes from the other parties

		// the value of a_0 is the input of the party.
		x1[0] = bi[l1]; //field->GetElement(bi[l1]);

		// generate random degree-T polynomial
		for (int i = 1; i < T + 1; i++) {
			// A random field element, uniform distribution
			x1[i] = field->Random();
		}

		helper->matrix_vand.MatrixMult(x1, y1, T + 1);

		// prepare shares to be sent
		for (int i = 0; i < N; i++) {
			sendBufsElements[i].push_back(y1[i]);
		}
	}
	//w,u (initied to 0)
	for(int l1=0; l1<2*dim; l1++)
	{
		// get the expected sizes from the other parties

		// the value of a_0 is the input of the party.
		x1[0] = ZZ_p(0); //field->GetElement(0);

		// generate random degree-T polynomial
		for (int i = 1; i < T + 1; i++) {
			// A random field element, uniform distribution
			x1[i] = field->Random();
		}

		helper->matrix_vand.MatrixMult(x1, y1, T + 1);

		// prepare shares to be sent
		for (int i = 0; i < N; i++) {
			sendBufsElements[i].push_back(y1[i]);
		}
	}
	int zero_cnt = 50000;
	if(m_partyID==0) //z
	{
		if(flag_print)
			cout<<"find ID 0"<<endl;
		for(int l1=0; l1<dim + zero_cnt; l1++)
		{
			// get the expected sizes from the other parties

			// the value of a_0 is the input of the party.
			x1[0] = ZZ_p(0); //field->GetElement(0);

			// generate random degree-T polynomial
			for (int i = 1; i < T + 1; i++) {
				// A random field element, uniform distribution
				x1[i] = field->Random();
			}

			helper->matrix_vand.MatrixMult(x1, y1, T + 1);

			// prepare shares to be sent
			for (int i = 0; i < N; i++) {
				sendBufsElements[i].push_back(y1[i]);
			}
		}
	}
	for(int i=1; i<N; i++)
		sizes[i] = dim*(dim+3);
	sizes[0] = dim*(dim+4) + zero_cnt;
	// -- convert shares to bytes
	for (int i = 0; i < N; i++) {
		sendBufsBytes[i].resize(sendBufsElements[i].size() * fieldByteSize);
		recBufBytes[i].resize(sizes[i] * fieldByteSize);
		for (int j = 0; j < sendBufsElements[i].size(); j++) {
			field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize),
					sendBufsElements[i][j]);
		}
	}

	if(flag_print)
		cout<<"before round"<<endl;
	helper->roundFunctionSync(sendBufsBytes, recBufBytes, 1);
	if(flag_print)
		cout<<"after round"<<endl;
	// -- convert received bytes to shares
	for (int i = 0; i < N; i++) {
		recBufElements[i].resize(sizes[i]);    
		for (int j = 0; j < sizes[i]; j++) {
			recBufElements[i][j] =
				field->bytesToElement(recBufBytes[i].data() + (j * fieldByteSize));
		}
	}	
	if(flag_print)
		cout<<"received"<<endl;
	//receive shares of other parties
	vector<ZZ_p> _chk;
	for(int i=0; i<N; i++) //each
	{
		int _tot = 0;
		shareOfA[i].resize(dim);
		//Ai
		for(int l1=0; l1<dim; l1++)
		{
			shareOfA[i][l1].resize(dim);
			for(int l2=0; l2<dim; l2++)
				shareOfA[i][l1][l2] = recBufElements[i][_tot++];
			if(flag_print)
			{
				cout<<"opening A"<<endl;
				helper->openShare(dim,shareOfA[i][l1],_chk);
			}
		}
		//Bi
		shareOfB[i].resize(dim);
		for(int l=0; l<dim; l++)
			shareOfB[i][l] = recBufElements[i][_tot++];
		//Wi
		shareOfW[i].resize(dim);
		for(int l=0; l<dim; l++)
			shareOfW[i][l] = recBufElements[i][_tot++];
		//Ui
		shareOfU[i].resize(dim);
		for(int l=0; l<dim; l++)
			shareOfU[i][l] = recBufElements[i][_tot++];
		//Z
		if(i==0)
		{
			shareOfZ.resize(dim);
			for(int l=0; l<dim; l++)
				shareOfZ[l] = recBufElements[i][_tot++];
			_zeroShares.resize(zero_cnt);
			for(int l=0; l<zero_cnt; l++)
				_zeroShares[l] = recBufElements[i][_tot++];
			if(flag_print)
			{
				cout<<"opening zero:"<<endl;
				helper->openShare(zero_cnt,_zeroShares,_chk);
			}
		}
	}


//the following are double number testing
	if(flag_print)
	{
		vector<ZZ_p> a1,a2,a3;
		/*a1[0] = 16;
		a1[1] = 256;
		a1[2] = 255;
		a1[3] = -256;
		a1[4] = 15;
		a1[5] = -324;*/
		TruncPR(shareOfB[0],a2);
		helper->openShare(dim,a2,a3);
		for(int i=0; i<dim; i++)
			cout<<bi[i]<<"->"<<a3[i]<<endl;
	}

	//testing softThres
	if(flag_print)
	{
		vector<ZZ_p> _res,_ans;
		helper->openShare(dim,shareOfB[0],_ans);
		cout<<"Thres:"<<endl;
		for(int i=0; i<dim; i++)
			cout<<_ans[i]<<endl;
		helper->openShare(dim,shareOfA[0][0],_ans);
		cout<<"A:"<<endl;
		for(int i=0; i<dim; i++)
			cout<<_ans[i]<<endl;
		SoftThres(bi,shareOfA[0][0],_res);
		cout<<"checking valid of softThres"<<endl;
		helper->openShare(dim,_res,_ans);
		cout<<"passed, checking softThres"<<endl;
		for(int i=0; i<dim; i++)
			cout<<Ai[0][i]<<"->"<<_ans[i]<<endl;
	}
	//start iteration.
	for(int _t=0; _t<iter; _t++)
	{
		if(flag_print)
			cout<<"Iteration "<<_t<<endl;
		//step 4(a)
		for(int i=0; i<N; i++) //repeat with each w[i]
		{
			//compute vector z - u[i]
			vector<ZZ_p> tmp,tmp1,tmp2;
			for(int j=0; j<dim; j++)
			{
				tmp1.push_back(shareOfZ[j] - shareOfU[i][j]);
				tmp2.push_back(rho);
			}
			doubleVecMult(tmp1,tmp2,tmp);
			//for(int j=0; j<dim; j++) //TODO: can we mult directly?
			//	tmp.push_back(shareOfB[i][j] + rho * (shareOfZ[j] - shareOfU[i][j]));
			//do matrix multiplication
			for(int j=0; j<dim; j++) //fill W[i][j]
			{
				//TODO: this is computed one by one. Can we batch?
				vector<ZZ_p> _t0;
				doubleVecMult(shareOfA[i][j],tmp,_t0);
				shareOfW[i][j] = 0;
				for(int l=0; l<dim; l++)
					shareOfW[i][j] = shareOfW[i][j] + _t0[l];
			}
		}
		if(flag_print)
			cout<<"4(b)"<<endl;
		//step 4(b)
		vector<ZZ_p> _t0,_t00;
		for(int j=0; j<dim; j++)
		{
			_t0.push_back(ZZ_p(0));
			for(int i=0; i<N; i++)
				_t0[j] = _t0[j] + shareOfW[i][j] + shareOfU[i][j];
			_t0[j] = _t0[j] * invN; //TODO: trunc?
		}
		TruncPR(_t0,_t00);
		if(flag_print)
		{
			cout<<"checking t0"<<endl;
			helper->openShare(dim,_t00,_t0);
			cout<<"t0 check passed"<<endl;
		}
		//compute the trunc of lambda * invN * invrho
		vector<ZZ_p> _t1,_t11,_t12;
		ZZ_p _tmp = lambda * invN;
		_tmp = conv<ZZ_p>(conv<ZZ>(_tmp) / (1ll<<_m));
		_tmp = _tmp * invrho;
		_tmp = conv<ZZ_p>(conv<ZZ>(_tmp) / (1ll<<_m));
		for(int j=0; j<dim; j++)
			_t12.push_back(_tmp);
		SoftThres(_t12, _t00, shareOfZ);
		//step 4(c)
		for(int i=0; i<N; i++)
		{
			for(int j=0; j<dim; j++)
				shareOfU[i][j] = shareOfU[i][j] + shareOfW[i][j] - shareOfZ[j];
		}
	}//end of iter
	if(res.size()<dim)
		res.resize(dim);
	for(int i=0; i<dim; i++)
		res[i] = shareOfZ[i];
	if(flag_print)
		cout<<"Used zero:"<<_zeroShareOffset<<"/"<<zero_cnt<<endl;
}


#endif /* PROTOCOLPARTY_H_ */
