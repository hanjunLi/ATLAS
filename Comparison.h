// Xiaoqi duan <duanxq17@mails.tsinghua.edu.cn>

#ifndef COMPARISION_H_
#define COMPARISION_H_

#include <libscapi/include/infra/Measurement.hpp>
#include "ProtocolParty.h"
#include "Interpolate.h"
#include <vector>

template<class FieldType>
class CompareGate 
{
private:
	ProtocolParty<FieldType> *helper;
        Interpolate<FieldType> interp;
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
	
	//compute 1/a[i] under double
	void doubleInverse(FieldType a,FieldType &res);

        // TODO:
        void readLassoInput(vector<vector<FieldType> > & Ai, vector<FieldType> &bi);
        void runLasso(int iter,FieldType lambda, FieldType rho, vector<vector<FieldType> > & Ai, vector<FieldType> &bi, vector<FieldType> &res);
};

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
        // TODO: generate bit sharings in runLasso()
	// int numBitShares = helper->numOfCompareGates + 100;
	// if(flag_print)
	// 	cout<<"generating bit shares:"<<numBitShares<<endl;
	// generateBitShares(numBitShares);
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
			interp.interpolate(tmpnum, tmp, OrVector[a[i].size()]);
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
	int N = helper->getN();
	int T = helper->getT();
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
	// for(int l1=0; l1<dim; l1++)
	// 	for(int l2=0; l2<dim; l2++)
	// 	{
	// 		// get the expected sizes from the other parties

	// 		// the value of a_0 is the input of the party.
	// 		x1[0] = Ai[l1][l2]; //field->GetElement(Ai[l1][l2]);

	// 		// generate random degree-T polynomial
	// 		for (int i = 1; i < T + 1; i++) {
	// 			// A random field element, uniform distribution
	// 			x1[i] = field->Random();
	// 		}

	// 		helper->matrix_vand.MatrixMult(x1, y1, T + 1);

	// 		// prepare shares to be sent
	// 		for (int i = 0; i < N; i++) {
	// 			sendBufsElements[i].push_back(y1[i]);
	// 		}
	// 	} 
	// //bi
	// for(int l1=0; l1<dim; l1++)
	// {
	// 	// get the expected sizes from the other parties

	// 	// the value of a_0 is the input of the party.
	// 	x1[0] = bi[l1]; //field->GetElement(bi[l1]);

	// 	// generate random degree-T polynomial
	// 	for (int i = 1; i < T + 1; i++) {
	// 		// A random field element, uniform distribution
	// 		x1[i] = field->Random();
	// 	}

	// 	helper->matrix_vand.MatrixMult(x1, y1, T + 1);

	// 	// prepare shares to be sent
	// 	for (int i = 0; i < N; i++) {
	// 		sendBufsElements[i].push_back(y1[i]);
	// 	}
	// }
	// //w,u (initied to 0)
	// for(int l1=0; l1<2*dim; l1++)
	// {
	// 	// get the expected sizes from the other parties

	// 	// the value of a_0 is the input of the party.
	// 	x1[0] = ZZ_p(0); //field->GetElement(0);

	// 	// generate random degree-T polynomial
	// 	for (int i = 1; i < T + 1; i++) {
	// 		// A random field element, uniform distribution
	// 		x1[i] = field->Random();
	// 	}

	// 	helper->matrix_vand.MatrixMult(x1, y1, T + 1);

	// 	// prepare shares to be sent
	// 	for (int i = 0; i < N; i++) {
	// 		sendBufsElements[i].push_back(y1[i]);
	// 	}
	// }
	int zero_cnt = 50000;
	// if(m_partyID==0) //z
	// {
	// 	if(flag_print)
	// 		cout<<"find ID 0"<<endl;
	// 	for(int l1=0; l1<dim + zero_cnt; l1++)
	// 	{
	// 		// get the expected sizes from the other parties

	// 		// the value of a_0 is the input of the party.
	// 		x1[0] = ZZ_p(0); //field->GetElement(0);

	// 		// generate random degree-T polynomial
	// 		for (int i = 1; i < T + 1; i++) {
	// 			// A random field element, uniform distribution
	// 			x1[i] = field->Random();
	// 		}

	// 		helper->matrix_vand.MatrixMult(x1, y1, T + 1);

	// 		// prepare shares to be sent
	// 		for (int i = 0; i < N; i++) {
	// 			sendBufsElements[i].push_back(y1[i]);
	// 		}
	// 	}
	// }
	// for(int i=1; i<N; i++)
	// 	sizes[i] = dim*(dim+3);
	// sizes[0] = dim*(dim+4) + zero_cnt;
	// // -- convert shares to bytes
	// for (int i = 0; i < N; i++) {
	// 	sendBufsBytes[i].resize(sendBufsElements[i].size() * fieldByteSize);
	// 	recBufBytes[i].resize(sizes[i] * fieldByteSize);
	// 	for (int j = 0; j < sendBufsElements[i].size(); j++) {
	// 		field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize),
	// 				sendBufsElements[i][j]);
	// 	}
	// }

	// if(flag_print)
	// 	cout<<"before round"<<endl;
	// helper->roundFunctionSync(sendBufsBytes, recBufBytes, 1);
	// if(flag_print)
	// 	cout<<"after round"<<endl;
	// // -- convert received bytes to shares
	// for (int i = 0; i < N; i++) {
	// 	recBufElements[i].resize(sizes[i]);    
	// 	for (int j = 0; j < sizes[i]; j++) {
	// 		recBufElements[i][j] =
	// 			field->bytesToElement(recBufBytes[i].data() + (j * fieldByteSize));
	// 	}
	// }	
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
#endif /* COMPARISON_H_ */
