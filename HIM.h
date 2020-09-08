// Hanjun Li <lihanjun1212@gmail.com>
// slightly modifeid from libscapi/include/primitives/Matrix.hpp

#ifndef HIM_H_
#define HIM_H_

#include <iostream>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/ZZ_p.h>
#include <NTL/GF2XFactoring.h>
#include <iostream>
#include <vector>
#include <array>
#include <libscapi/include/primitives/Mersenne.hpp>

#include "Share.h"

// TODO: use only HIM in protocol, then remove VDM.

using namespace std;
using namespace NTL;

template <typename FieldType>
class HIMp {
private:
    int m_n,m_m;
    FieldType** m_matrix;
    TemplateField<FieldType> *field;
public:

    /**
     * This method allocate m-by-n matrix.
     * m rows, n columns.
     */
    HIMp(int m, int n, TemplateField<FieldType> *field);

    HIMp();

    /**
     * This method is a construction of a hyper-invertible m-by-n matrix M over a finite field F with |F| ≥ 2n.
     * Let α1,...,αn , β1,...,βm denote fixed distinct elements in F according the vectors alpha and beta,
     * and consider the function f:Fn → Fm,
     * mapping (x1,...,xn) to (y1,...,ym) such that the points (β1,y1),...,(βm,ym) lie on the polynomial g(·)
     * of degree n−1 defined by the points (α1,x1),...,(αn,xn).
     * Due to the linearity of Lagrange interpolation, f is linear and can be expressed as a matrix:
     * M = {λi,j} j=1,...n i=1,...,m
     * where λ i,j = {multiplication}k=1,..n (βi−αk)/(αj−αk)
     */
    FieldType** InitHIMByVectors(vector<FieldType> &alpha, vector<FieldType> &beta);

    FieldType** InitHIMVectorAndsizes(vector<FieldType> &alpha, int n, int m);

    /**
     * This method create vectors alpha and beta,
     * and init the matrix by the method InitHIMByVectors(alpha, beta).
     */
    FieldType** InitHIM();

    /**
     * This method print the matrix
     */
    void Print();

    /**
     * matrix/vector multiplication.
     * The result is the answer vector.
     */
    void MatrixMult(std::vector<FieldType> &vector, std::vector<FieldType> &answer);

    void allocate(int m, int n, TemplateField<FieldType> *field);

    virtual ~HIMp();

  // added functions
  int getNRows() {return m_m;}
  int getNCols() {return m_n;}
  void getRow(int idx, vector<FieldType>& row) {
    row.resize(m_n);
    for (int i=0; i<m_n; i++) {
      row[i] = m_matrix[idx][i];
    }
  }
  void MatrixMultShares(vector<Share<FieldType>>& inVec,
                        vector<Share<FieldType>>& outVec) {
    assert(inVec.size() == m_n);
    outVec.clear(); outVec.resize(m_m);
    Share<FieldType> temp1;
    for(int i = 0; i < m_m; i++)
    {
        for(int j=0; j < m_n; j++)
        {
          temp1 = inVec[j];
          temp1.multByConst(m_matrix[i][j]);
          outVec[i] += temp1;
        }
    }    
  }
};



template <typename FieldType>
HIMp<FieldType>::HIMp(){}

template <typename FieldType>
HIMp<FieldType>::HIMp(int m, int n, TemplateField<FieldType> *field) {
    // m rows, n columns
    this->m_m = m;
    this->m_n = n;
    this->field = field;
    this->m_matrix = new FieldType*[m_m];

    for (int i = 0; i < m_m; i++)
    {
        m_matrix[i] = new FieldType[m_n];
    }
}

template <typename FieldType>
FieldType** HIMp<FieldType>::InitHIMByVectors(vector<FieldType> &alpha, vector<FieldType> &beta)
{
    FieldType lambda;


    int m = beta.size();
    int n = alpha.size();
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            // lambda = 1
            lambda = *(field->GetOne());

            // compute value for matrix[i,j]
            for (int k = 0; k < n; k++)
            {
                if (k == j)
                {
                    continue;
                }

                lambda *= ((beta[i]) - (alpha[k])) / ((alpha[j]) - (alpha[k]));
            }

            // set the matrix
            (m_matrix[i][j]) = lambda;
        }
    }
    return m_matrix;
}


template <typename FieldType>
FieldType** HIMp<FieldType>::InitHIMVectorAndsizes(vector<FieldType> &alpha, int n, int m)
{
    FieldType lambda;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            // lambda = 1
            lambda = *(field->GetOne());

            // compute value for matrix[i,j]
            for (int k = 0; k < n; k++)
            {
                if (k == j)
                {
                    continue;
                }

                lambda *= ((alpha[n+i]) - (alpha[k])) / ((alpha[j]) - (alpha[k]));
            }

            // set the matrix
            (m_matrix[i][j]) = lambda;
        }
    }
    return m_matrix;
}



template <typename FieldType>
void HIMp<FieldType>::allocate(int m, int n, TemplateField<FieldType> *field)
{
    // m rows, n columns
    this->m_m = m;
    this->m_n = n;
    this->field = field;
    this->m_matrix = new FieldType*[m_m];
    for (int i = 0; i < m_m; i++)
    {
        m_matrix[i] = new FieldType[m_n];
    }
}

template <typename FieldType>
FieldType** HIMp<FieldType>::InitHIM()
{
    int i;
    vector<FieldType> alpha(m_n);
    vector<FieldType> beta(m_m);

    // check if valid
    if (256 <= m_m+m_n)
    {
        cout << "error";
    }

    // Let alpha_j and beta_i be arbitrary field elements
    for (i = 0; i < m_n; i++)
    {
        alpha[i] = field->GetElement(i);
    }

    for (i = 0; i < m_m; i++)
    {
        beta[i] = field->GetElement(m_n+i);
    }

    return(InitHIMByVectors(alpha,beta));
}

template <typename FieldType>
void HIMp<FieldType>::Print()
{
    for (int i = 0; i < m_m; i++) {
        for (int j = 0; j < m_n; j++) {
            cout << (m_matrix[i][j]) << " ";
        }

        cout << " " << '\n';
    }

}

template <typename FieldType>
void HIMp<FieldType>::MatrixMult(std::vector<FieldType> &vector, std::vector<FieldType> &answer)
{
    FieldType temp1;
    for(int i = 0; i < m_m; i++)
    {
        // answer[i] = 0
        answer[i] = *(field->GetZero());

        for(int j=0; j < m_n; j++)
        {
            temp1 = m_matrix[i][j] * vector[j];
            //answer[i] = answer[i] + temp1;
            answer[i] += temp1;
        }
    }
}

template <typename FieldType>
HIMp<FieldType>::~HIMp() {
    for (int i = 0; i < m_m; i++)
    {
        delete[] m_matrix[i];
    }
    delete[] m_matrix;
}

template<typename FieldType>
class VDMt {
private:
    int m_n,m_m;
    FieldType** m_matrix;
    TemplateField<FieldType> *field;
public:
    VDMt(int n, int m, TemplateField<FieldType> *field);
    VDMt() {};
    ~VDMt();
    void InitVDMTranspose();
    void Print();
    void MatrixMult(std::vector<FieldType> &vector, std::vector<FieldType> &answer, int length);

    void allocate(int n, int m, TemplateField<FieldType> *field);

  // added functions
  int getNRows() {return m_m;}
  int getNCols() {return m_n;}
  void getRow(int idx, vector<FieldType>& row) {
    row.resize(m_n);
    for (int i=0; i<m_n; i++) {
      row[i] = m_matrix[idx][i];
    }
  }
};


template<typename FieldType>
VDMt<FieldType>::VDMt(int n, int m, TemplateField<FieldType> *field) {
    this->m_m = m;
    this->m_n = n;
    this->field = field;
    this->m_matrix = new FieldType*[m_n];
    for (int i = 0; i < m_n; i++)
    {
        m_matrix[i] = new FieldType[m_m];
    }
}

template<typename FieldType>
void VDMt<FieldType>::allocate(int n, int m, TemplateField<FieldType> *field) {

    this->m_m = m;
    this->m_n = n;
    this->field = field;
    this->m_matrix = new FieldType*[m_n];
    for (int i = 0; i < m_n; i++)
    {
        m_matrix[i] = new FieldType[m_m];
    }
}

template<typename FieldType>
void VDMt<FieldType>::InitVDMTranspose() {
    vector<FieldType> alpha(m_m);
    for (int i = 0; i < m_m; i++) {
        alpha[i] = field->GetElement(i + 1);
    }

    for (int i = 0; i < m_m; i++) {
        m_matrix[0][i] = *(field->GetOne());
        for (int k = 1; k < m_n; k++) {
            m_matrix[k][i] = m_matrix[k-1][i] * (alpha[k]);
        }
    }
}

/**
 * the function print the matrix
 */
template<typename FieldType>
void VDMt<FieldType>::Print()
{
    for (int i = 0; i < m_m; i++)
    {
        for(int j = 0; j < m_n; j++)
        {
            cout << (m_matrix[i][j]) << " ";

        }
        cout << " " << '\n';
    }

}

template<typename FieldType>
void VDMt<FieldType>::MatrixMult(std::vector<FieldType> &vector, std::vector<FieldType> &answer, int length)
{
    for(int i = 0; i < length; i++)
    {
        // answer[i] = 0
        answer[i] = *(field->GetZero());

        for(int j=0; j < m_m; j++)
        {
            answer[i] += (m_matrix[i][j] * vector[j]);
        }
    }

}


//
template<typename FieldType>
VDMt<FieldType>::~VDMt() {
    for (int i = 0; i < m_n; i++) {
        delete[] m_matrix[i];
    }
    delete[] m_matrix;
}

#endif /* HIM_H_ */
