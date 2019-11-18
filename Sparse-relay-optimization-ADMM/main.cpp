//
//  main.cpp
//  Sparse-relay-optimization-ADMM
//
//  Created by Ayoub Saab  on 2019-07-08.
//  Copyright © 2019 Ayoub Saab. All rights reserved.
//



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

//system paramters
const unsigned int UEs=2,relays=2;
const unsigned int antenna[relays]={2,2};
const unsigned int power[relays]={2,2};

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

int main()
{
    //initialize the system paramters
    vector<unsigned int> N (relays);
    for(int i=0;i<N.size();i++)
    {
        N.at(i) = antenna[i];
    }
    vector <unsigned int> Ns(relays);
    for (int i=0;i<Ns.size();i++)
    {
        Ns.at(i) = N.at(i)*N.at(i);
    }
    
    //initialize the relay-power budgets
    vector<unsigned int> P(relays);
    for (int i=0;i<P.size();i++)
    {
        P.at(i)=power[i];
    }
    
    //initialize the distortionless constraints
    MatrixXcd c(UEs,1);
    c.setOnes();
    
    //standard deviation at different stations
    const  float sigmaChannel = sqrt(0.5); //channel I/Q std deviation, total = 2*sigmaChannel^2 = 1;
    vector<double> sigmaRelay (relays);
    for(int i=0;i<sigmaRelay.size();i++)
    {
        sigmaRelay.at(i) = sqrt(0.5)*sqrt(0.1);
    }
    //vector<float> sigmaDUE (UEs);
    
    
    //simulation paramters
    const int monteCarlo = pow(10,0);  //number of Monte Carlo simulations
    const double rho = 200.0;         //penalty parameter
    const unsigned int lambda[relays] = {100,100}; //regularization parameter
    double eps_abs = 0;   //absolute tolerance metric
    double eps_rel = pow(10,-4);    //relative tolerance metric
    
    const int size = vectorSum(N,0,relays-1); //sum_{i=1}^L N_l
    const int size2 = vectorSum(Ns,0,relays-1); //sum_{i=1}^L N_l^2
    
    for (int sim = 0 ; sim<monteCarlo ; sim++)
    {
        
        //generate channels, matrix version
        MatrixXcd h(size,UEs);
        MatrixXcd g(size,UEs);
        generateChannels(h,g,sigmaChannel);
        cout<<"Backward channel H: "<<endl;
        cout<<h<<endl<<endl;
        cout<<"Forward channel G: "<<endl;
        cout<<g<<endl<<endl;
        
        //generate Psi: Matrix to express relay-power constraints
        MatrixXcd Psi (size2,size2);
        Psi.setZero();
        generatePsi(Psi,h,sigmaRelay,N,Ns);
        //cout<<"PSI MATRIX: "<<endl;
        //cout<<Psi<<endl;
        
        //generate small delta matrix: D
        MatrixXcd D(size2,UEs*UEs);
        D.setZero();
        generateD(D,h,g,N,Ns);
        //cout<<"D MATRIX:"<<endl;
        //cout<<D<<endl;
        
        //generate PHI MATRIX
        MatrixXcd Phi (size2,UEs);
        for (int k =0 ; k<UEs ; k++)
        {
            for (int i =0;i<size2; i++)
            {
                Phi(i,k) =  D(i,k*UEs + k);
            }
        }
        //cout<<"PHI MATRIX:"<<endl;
        //cout<<Phi<<endl;
        
        //generate DELTA MATRIX
        MatrixXcd Delta (size2,size2);
        Delta.setZero();
        generateDelta(Delta,D,size2);
        //cout<<"DELTA MATRIX:"<<endl;
        //cout<<Delta<<endl;
        
        //generate the G Matrix
        MatrixXcd G(size2,size2);
        G.setZero();
        generateG(G, g, sigmaRelay, N, Ns, size2);
        //cout<<"G MATRIX:"<<endl;
        //cout<<G<<endl;
        
        //generate the Theta Matrix: THETA = G + DELTA
        MatrixXcd Theta(size2,size2);
        Theta = G + Delta;
        //cout<<"THETA MATRIX:"<<endl;
        //cout<<Theta<<endl;
        
        Phi = Psi.pow(-0.5)*Phi;
        //cout<<"New PHI = PSI^(-0.5)*PHI "<<endl;
        //cout<<Phi<<endl;
        
        //begin ADMM algorithm
        MatrixXcd Q(size2,size2);
        MatrixXcd Qinverse(size2,size2);
        MatrixXcd Q_tilde(UEs,UEs);
        MatrixXcd I(size2,size2);
        I.setIdentity();
        MatrixXcd nu (UEs,1);
        MatrixXcd theta(size2,1);
        MatrixXcd w(size2,1);
        MatrixXcd mu(size2,1);
        MatrixXcd rk(size2,1);
        MatrixXcd sk(size2,1);
        MatrixXcd sum_sk(size2,1);
        double eps_pri,eps_dual;
        
        //algorithm initialization
        Q = Psi.pow(-0.5)*Theta*Psi.pow(-0.5) + (rho/2.0)*I; //Q matrix
        Qinverse = Q.inverse();         // Q^-1
        theta.setRandom(); //initialize theta randomly
        mu.setZero();      //initialize the dual variable to all zeros
        
        //initialize the tolerance constraints
        rk = -theta; //primal residual
        sk = (-rho/2)*(theta);  //dual residual
        sum_sk = sk; //sum of dual residuals
        eps_pri = sqrt(size2)*eps_abs + eps_rel*theta.norm();
        eps_dual = sqrt(size2)*eps_abs + eps_rel*mu.norm();
        
        int count = 0;
        
        //algorithm execution
        while( rk.norm() > eps_pri || sk.norm() > eps_dual )
        {
            //ADMM step1, update w
            Q_tilde = Phi.adjoint()*Qinverse*Phi;
            nu = Q_tilde.inverse()*(c+Phi.adjoint()*Qinverse*mu-(rho/2.0)*Phi.adjoint()*Qinverse*theta);
            w = Qinverse*(rho/2.0*theta-mu+Phi*nu);
            
            
            //ADMM step2, update theta by solving L parallel subproblems
            for(int l=0;l<relays;l++)
            {
                //populate a_l
                MatrixXcd a_l (antenna[l]*antenna[l],1);
                int offset = vectorSum(Ns, 0, l-1);
                for(int i=0;i<antenna[l]*antenna[l];i++)
                {
                    a_l(i,0) = (rho/2.0)*w(i+offset,0) + mu(i+offset,0);
                    
                }
                
                if(a_l.norm() - lambda[l] <=0)
                {
                    //set theta_l to zeros
                    for (int j=0;j<antenna[l]*antenna[l];j++)
                    {
                        theta(j+offset,0) = 0;
                    }
                }
                
                else
                {
                    double eta_l;
                    eta_l=max(0.0,(a_l.norm()-lambda[l])/sqrt(P[l])-rho/2.0);
                    for (int i =0 ; i<antenna[l]*antenna[l]; i++)
                    {
                        theta(i+offset,0) = a_l(i,0)*(a_l.norm()-lambda[l])/(a_l.norm()*(rho/2.0+eta_l));
                    }
                }
            }
            
            //ADMM step 3, update dual variable
            mu = mu + (rho/2.0)*(w-theta);
            
            //update loop tolerance metrics
            rk = w-theta;                      //primal residual
            sk = (-rho/2.0)*theta - sum_sk;      //dual residual
            sum_sk = sum_sk + sk;              //sum of dual residuals
            eps_pri = sqrt(size2)*eps_abs + eps_rel*max(w.norm(),theta.norm()) ; //primal tolerance threshold
            eps_dual = sqrt(size2)*eps_abs + eps_rel*mu.norm();                 //dual torelance threshold
            
            
            //compute the objective function at each iteration
            double sumMixedNorm = 0;
            for (int l=0 ; l<relays ; l++)
            {
                MatrixXcd temp (antenna[l]*antenna[l],1);
                int offset = vectorSum(Ns, 0, l-1);
                for (int j=0;j<antenna[l]*antenna[l];j++)
                {
                    temp(j,0) = theta(j+offset,0);
                }
                sumMixedNorm = sumMixedNorm + lambda[l]*temp.norm();
            }
            
            //display the value of the objective function
            //MatrixXcd temp(1,1);
            //temp = theta.adjoint()*Psi.pow(-0.5)*Theta*Psi.pow(-0.5)*theta;
            //cout<<"Sum mixed norm "<<sumMixedNorm<<endl;
            //cout<<"Iteration: "<<count<<", Objective function: "<<temp.real().trace()+sumMixedNorm<<endl<<endl;
            count++;
        }
        //end ADMM algorithm
        
        cout<<"The solution vector is: "<<endl;
        cout<<theta<<endl<<endl;
        //cout<<w<<endl;
        cout<<"The algorithm converged in "<<count<<" iterations."<<endl<<endl;

        
        //compute final-achieved SINR
        //              for l=1:L
        //            W3(idxN(1,l):idxN(2,l),1:N(l))= ...
        //            reshape(Psi(idxNs(1,l):idxNs(2,l),idxNs(1,l):idxNs(2,l))^(-0.5)*theta(idxNs(1,l):idxNs(2,l)),N(l),N(l));
        //        end
        //        sinr_a_c = [sinr_a_c compute_sinr(W3,K,L,N,idxN,G,H,sigma_r,sigma_ue,sigma_due)];
        
    } //end Monte Carlo for-loop
    
    return 0;
    
}


