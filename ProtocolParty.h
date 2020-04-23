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
// #include "Comparison.h"
#include <cmath>

#define flag_print true
#define flag_print_output false

using namespace std;
using namespace std::chrono;

// template<class FieldType> class CompareGate;

template <class FieldType>
class ProtocolParty : public Protocol, public HonestMajority, MultiParty {

	public:
		// -- polynomial functionalities
		Interpolate<FieldType> interp; // evaluate (O(n)), interpolate polynomials (O(n^2))
		vector<FieldType> gateShareArr; // my share of each gate (1 share each)
		// int numOfCompareGates;
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
		// CompareGate<FieldType> *comper;
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
		bool preparationPhase(int _cnt = 0);
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
		
		// int processComp();
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

                // -- external access
                int getN() {return N;}
                int getT() {return T;}
                int getMyId() {return m_partyId;}
                TemplateField<FieldType>* getField() {return field;}
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
		// this->numOfCompareGates = circuit.getNrOfCompareGates();
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
	// cout<<"Gates:"<<numOfMultGates + numOfCompareGates<<endl;
	// cout<<"Ave time:"<<(double)tottme / times / (numOfMultGates + numOfCompareGates);
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
	// comper->runLasso(10,FieldType(1ll<<(7)),FieldType(1ll<<(7)),_Ai,_bi,res);
	// //if(flag_print)
	// {
	// 	cout<<"Lasso finished"<<endl;
	// cout<<"Used single:"<<_singleSharesOffset<<"/"<<_singleSharesArray.size()<<","<<(double)_singleSharesOffset / _singleSharesArray.size()<<",left:"<<_singleSharesArray.size() - _singleSharesOffset << endl;
	// 	cout<<"Used double:"<<2 * _doubleSharesOffset<<"/"<<_doubleSharesArray.size()<<","<<(double)2*_doubleSharesOffset / _doubleSharesArray.size()<<",left:"<<-2*_doubleSharesOffset + _doubleSharesArray.size()<<endl;
	// }
	// abort();
	int countNumMult = _doubleSharesOffset;
	int numOfLayers = circuit.getLayers().size();
	for (int i = 0; i < numOfLayers - 1; i++) {
		if(flag_print)
			cout<<"Now at layer "<<i<<"/"<<numOfLayers<<endl;
		currentCircuitLayer = i;
		// send the index of the current mult gate
		processNotMult();
		// processComp();
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
		if(flag_print)
			cout<<"Base case"<<endl;
		verifyVec(aShares, bShares, cShare);
		return;
	} // else : recursive case:
	// -- divide into K groups
	int groupSize = (totalLength + _K -1 )/ _K;
        aShares.resize( groupSize * _K, field->GetElement(0) );
        bShares.resize( groupSize * _K, field->GetElement(0) );
        totalLength = groupSize * _K;
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

template <class FieldType> bool ProtocolParty<FieldType>::preparationPhase(int tcnt) {
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
                4 * (keySize + verifyIterations) + 2 * _K + nCompressions + 2;
		// + 60 * eleSize *numOfCompareGates;
	//cout<<"should have num:"<<numSingleShares<<endl;
	numSingleShares += tcnt;
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
                this->numOfMultGates + (nCompressions+1)*_K*2;
        // + 100 * eleSize * this->numOfCompareGates;
	//used for Lasso?
	numDoubleShares += tcnt;
	_doubleSharesOffset = 0;
	offlineDNForMultiplication(numDoubleShares);

	if(flag_print)
		cout << "generated doubles: " << numDoubleShares << endl;
	//a brute force way for generating Or polys are applied.
	//can be implemented in a more efficient way.

	//put it here to let bit share to be generated
	// this->comper = new CompareGate<FieldType>(this,eleSize,m_partyId,field);
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


// 	template<class FieldType>
// int ProtocolParty<FieldType>::processComp()
// {
// 	//step 1: go through gates
// 	int cnt = 0;
// 	vector<FieldType> a,b,res;
// 	for(int k = circuit.getLayers()[currentCircuitLayer];
// 			k<circuit.getLayers()[currentCircuitLayer + 1]; k++)
// 	{
// 		auto gate = circuit.getGates()[k];
// 		if(gate.gateType == COMPARE)
// 		{
// 			cnt++;
// 			a.push_back(gateShareArr[gate.input1]);
// 			b.push_back(gateShareArr[gate.input2]);	
// 			//tmp.push_back(gate);
// 		}
// 	}
// 	comper->compRandom(a,b,res);
// 	cnt = 0;
// 	for(int k = circuit.getLayers()[currentCircuitLayer];
// 			k<circuit.getLayers()[currentCircuitLayer + 1]; k++)
// 	{
// 		auto gate = circuit.getGates()[k];
// 		if(gate.gateType == COMPARE)
// 			gateShareArr[gate.output] = res[cnt++];
// 	}
// 	if(flag_print)
// 	{
// 		cout<<"compRandom ended"<<endl;
// 		cout<<"Used single:"<<_singleSharesOffset<<"/"<<_singleSharesArray.size()<<","<<(double)_singleSharesOffset / _singleSharesArray.size()<<",left:"<<_singleSharesArray.size() - _singleSharesOffset << endl;
// 		cout<<"Used double:"<<2 * _doubleSharesOffset<<"/"<<_doubleSharesArray.size()<<","<<(double)2*_doubleSharesOffset / _doubleSharesArray.size()<<",left:"<<-2*_doubleSharesOffset + _doubleSharesArray.size()<<endl;
// 	}
// 	return cnt;
// }

template <class FieldType> ProtocolParty<FieldType>::~ProtocolParty() {
	protocolTimer->writeToFile();
	delete protocolTimer;
	delete field;
	delete timer;
	// delete comm;
}


#endif /* PROTOCOLPARTY_H_ */
