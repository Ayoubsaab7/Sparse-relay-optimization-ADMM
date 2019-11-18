//
//  functions.h
//  Sparse-relay-optimization-ADMM
//
//  Created by Ayoub Saab  on 2019-11-18.
//  Copyright © 2019 Ayoub Saab. All rights reserved.
//

#ifndef functions_h
#define functions_h

#include "parameters.h"

//preprocessor directives
#include <iostream>
#include <string>
#include <cmath>
#include<vector>
#include<complex>
#include<random>
#include<time.h>

//necessary library for Matrix manipulations
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/MatrixFunctions"
using Eigen::MatrixXd;
using namespace Eigen;

using namespace std;



int vectorSum(vector<unsigned int> V, int startIndex,  int endIndex)
{
    int sum = 0;
    
    if(startIndex<=endIndex)
    {
        for(int i = startIndex; i<=endIndex ; i++)
        {
            sum = sum + V[i];
        }
    }
    
    return sum;
}


double max(const double x,const double y)
{
    double max;
    if(x>=y)
    {
        max = x;
    }
    else{ max = y;}
    
    return max;
}


void Kroneckerproduct(const MatrixXcd  A, const MatrixXcd B, MatrixXcd &C)
{
    long rowA = A.rows();
    long rowB = B.rows();
    long colA = A.cols();
    long colB = B.cols();
    
    long startRow,startCol;
    
    for (int i=0 ; i<rowA ; i++)
    {
        for(int j=0;j<colA ; j++)
        {
            startRow = i*rowB;
            startCol = j*colB;
            for (int k =0 ;k<rowB;k++ )
            {
                for(int l=0 ; l<colB ;l++)
                {
                    C(startRow+k,startCol+l) = A(i,j)*B(k,l);
                }
            }
        }
    }
}


void generateChannels(MatrixXcd& h, MatrixXcd& g, const float sigmaChannel)
{
    std::normal_distribution<double> distribution(0,sigmaChannel);
    //std::default_random_engine generator;
    std::default_random_engine generator(time(NULL));
    
    //determine the size
    int temp=0;
    for (int i=0; i<relays;i++)
    {
        temp = temp+antenna[i];
    }
    
    for (int k=0;k<UEs;k++)
    {
        for( int n=0;n<temp;n++)
        {
            h(n,k) = std::complex<double>(distribution(generator),distribution(generator));    //backward channel
            g(n,k) = std::complex<double>(distribution(generator),distribution(generator)); //forward channel: complex conjugate
        }
    }
    
}


void generatePsi(MatrixXcd & Psi, const MatrixXcd h, const vector<double> sigmaRelay, const vector <unsigned int> N, const vector<unsigned int> Ns)
{
    
    //generate Psi_l, l=0,...,L-1
    for (int l=0;l<relays;l++)
    {
        MatrixXcd Psi_l (antenna[l]*antenna[l],antenna[l]*antenna[l]);
        MatrixXcd temp (antenna[l],UEs);
        MatrixXcd temp2(antenna[l],antenna[l]);
        MatrixXcd eye_l (antenna[l],antenna[l]);
        eye_l.setIdentity();
        
        int offset = vectorSum(N,0,l-1);
        
        //populate the temp matrix:
        for (int k=0;k<UEs;k++)
        {
            for (int j=0;j<antenna[l];j++)
            {
                temp(j,k) = h(j+offset,k);
            }
        }
        
        temp2 = temp*temp.adjoint() + 2*sigmaRelay[l]*sigmaRelay[l]*eye_l;
        
        //compute the Kronecker product of temp2 and eye_l and store it in Psi_l
        Kroneckerproduct(temp2, eye_l , Psi_l);
        
        //populate the Psi matrix with Psi_l, for l=0,1,...,L-1 in a block diagonal fashion
        offset = vectorSum(Ns, 0, l-1);
        for (int i=0;i<antenna[l]*antenna[l]; i++)
        {
            for(int j=0; j<antenna[l]*antenna[l]; j++)
            {
                Psi(i+offset,j+offset) = Psi_l(i,j);
            }
        }
    }
    
}

void generateD(MatrixXcd &D,const MatrixXcd h, const MatrixXcd g, const vector <unsigned int> N, const vector <unsigned int> Ns)
{
    
    for (int l=0;l<relays;l++)
    {
        int offset2 = vectorSum(Ns, 0, l-1);
        MatrixXcd temp (antenna[l]*antenna[l],1);
        for (int k=0 ; k<UEs ; k++)
        {
            
            //get g_kl
            MatrixXcd g_kl (antenna[l],1);
            int offset = vectorSum(N, 0, l-1);
            for (int n=0;n<antenna[l];n++)
            {
                g_kl(n,0) = g(n+offset,k);
            }
            
            for (int j=0 ; j<UEs ; j++)
            {
                //get h_lj
                MatrixXcd h_lj (antenna[l],1);
                for(int m=0 ; m<antenna[l] ; m++)
                {
                    h_lj(m,0) = h(offset+m,j);
                }
                
                //compute the kronecker prodcut of h_lj* and g_kl and store it in temp
                Kroneckerproduct(h_lj.conjugate(), g_kl, temp);
                
                //store all the results in D
                for(int z=0 ; z<antenna[l]*antenna[l] ; z++)
                {
                    D(z+offset2,j+k*UEs) = temp(z,0);
                }
                
            }
        }
    }
    
}


void generateDelta(MatrixXcd& Delta, const MatrixXcd D, const int size2)
{
    for (int k=0;k<UEs;k++)
    {
        MatrixXcd delta_k_j(size2,1);
        
        for (int j=0;j<UEs;j++)
        {
            if(j!=k)
            {
                //populate
                for(int i =0;i<size2;i++)
                {
                    delta_k_j(i,0) = D(i,j+UEs*k);
                }
                //add the outer product of delta_kj
                Delta = Delta + delta_k_j*delta_k_j.adjoint(); // Delta = sum_{k=1}^{UEs} Delta_k
            }
        }
    }
}


void generateG(MatrixXcd& G, const MatrixXcd g,  const vector<double> sigmaRelay,const vector <unsigned int> N, const vector <unsigned int> Ns, const int size2)
{
    for (int k=0;k<UEs;k++)
    {
        MatrixXcd G_k(size2,size2);
        G_k.setZero();
        
        for (int l=0;l<relays; l++)
        {
            MatrixXcd g_kl(antenna[l],1);
            int offset = vectorSum(N,0,l-1);
            //get g_kl
            for (int m=0;m<antenna[l];m++)
            {
                g_kl(m,0) = g(m+offset,k);
            }
            
            //get G_kl
            MatrixXcd G_kl(antenna[l]*antenna[l],antenna[l]*antenna[l]);
            G_kl.setZero();
            
            for (int j=0;j<antenna[l];j++)
            {
                MatrixXcd eye_j (antenna[l],1);
                eye_j.setZero();
                eye_j(j,0) = 1;
                MatrixXcd kron (antenna[l]*antenna[l],1);
                Kroneckerproduct(eye_j, g_kl, kron);
                G_kl = G_kl + kron*kron.adjoint() ;
                
            }
            G_kl = G_kl*2*sigmaRelay[l]*sigmaRelay[l];
            
            //G_k = blkdiag{G_k1, G_k2, ... , G_kL}
            int offset2 = vectorSum(Ns, 0, l-1);
            for (int n=0;n<antenna[l]*antenna[l];n++)
            {
                for (int m=0;m<antenna[l]*antenna[l];m++)
                    G_k(n+offset2,m+offset2) = G_kl(n,m);
            }
        }
        G = G + G_k; //G =
    }
}


#endif /* functions_h */
