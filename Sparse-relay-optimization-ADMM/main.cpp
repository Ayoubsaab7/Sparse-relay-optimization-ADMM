//
//  main.cpp
//  Sparse-relay-optimization-ADMM
//
//  Created by Ayoub Saab  on 2019-07-08.
//  Copyright © 2019 Ayoub Saab. All rights reserved.
//


//preprocessor directives
#include "functions.h"
#include <iostream>
using namespace std;


int main()
{
    //initialize the system parameters
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
    vector<double> sigmaRelay (relays);
    for(int i=0;i<sigmaRelay.size();i++)
    {
        sigmaRelay.at(i) = sqrt(0.5)*sigmaR;
    }
    //vector<float> sigmaDUE (UEs);
    
    
    const int size = vectorSum(N,0,relays-1); //sum_{i=1}^L N_l
    const int size2 = vectorSum(Ns,0,relays-1); //sum_{i=1}^L N_l^2
    
    for (int sim = 0 ; sim<monteCarlo ; sim++)
    {
        //generate channels
        MatrixXcd h(size,UEs);
        MatrixXcd g(size,UEs);
        generateChannels(h,g,sigmaChannel);
        cout<<"Backward channel H: "<<endl;
        cout<<h<<endl<<endl;
        cout<<"Forward channel G: "<<endl;
        cout<<g<<endl<<endl;
        
        /* GENERATE THE REQUIRED MATRICES FOR THE OPTIMIZATION FUNCTION:
         
         1)PSI, 
         2)PHI 
         3)DELTA, 
         4)G,
         5)THETA=DELTA+G
         
         */
        
        //generate Psi: Matrix to express relay-power constraints
        MatrixXcd Psi (size2,size2);
        Psi.setZero();
        generatePsi(Psi,h,sigmaRelay,N,Ns);
        
        //generate small delta matrix: D
        MatrixXcd D(size2,UEs*UEs);
        D.setZero();
        generateD(D,h,g,N,Ns);

        //generate PHI MATRIX from D
        MatrixXcd Phi (size2,UEs);
        for (int k =0 ; k<UEs ; k++)
        {
            for (int i =0;i<size2; i++)
            {
                Phi(i,k) =  D(i,k*UEs + k);
            }
        }
        
        //generate DELTA MATRIX
        MatrixXcd Delta (size2,size2);
        Delta.setZero();
        generateDelta(Delta,D,size2);
        
        //generate the G Matrix
        MatrixXcd G(size2,size2);
        G.setZero();
        generateG(G, g, sigmaRelay, N, Ns, size2);
        
        //generate the Theta Matrix: THETA = G + DELTA
        MatrixXcd Theta(size2,size2);
        Theta = G + Delta;
        
        Phi = Psi.pow(-0.5)*Phi;
        //cout<<"New PHI = PSI^(-0.5)*PHI "<<endl;
        //cout<<Phi<<endl;
        
        
        /* BEGIN THE OPTIMIZATION PROCESS THROUGH THE ALTERNATING DIRECTION
         
         METHOD OF MULTIPLIERS (ADMM)
        
         */
        
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


