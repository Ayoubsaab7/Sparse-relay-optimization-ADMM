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
#include<iomanip>
#include<ctime>
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
    c=c*distortionlessConstraint;
    
    //standard deviation at different stations
    vector<double> sigmaRelay (relays);
    for(int i=0;i<sigmaRelay.size();i++)
    {
        sigmaRelay.at(i) = sqrt(0.5)*sigmaR;
    }
    //vector<float> sigmaDUE (UEs);
    
    
    const int size = vectorSum(N,0,relays-1); //sum_{i=1}^L N_l
    const int size2 = vectorSum(Ns,0,relays-1); //sum_{i=1}^L N_l^2
    
    
    /* BEGIN MONTE CARLO SIMULATION */
    
    for (int sim = 0 ; sim<monteCarlo ; sim++)
    {
        /* GENERATE THE WIRELESS CHANNELS H (backward) and G (foreward) */
        
        //generate channels
        MatrixXcd h(size,UEs);
        MatrixXcd g(size,UEs);
        generateChannels(h,g,sigmaChannel);
        cout<<"Backward channel H coefficients: "<<endl;
        cout<<h<<endl<<endl;
        cout<<"Forward channel G coefficients: "<<endl;
        cout<<g<<endl<<endl;
        
        /* GENERATE THE REQUIRED MATRICES FOR THE OPTIMIZATION FUNCTION:
         
         1)PSI, 
         2)PHI 
         3)DELTA, 
         4)G,
         5)THETA=DELTA+G
         
         */
        
        //1) generate Psi: Matrix to express relay-power constraints
        MatrixXcd Psi (size2,size2);
        Psi.setZero();
        generatePsi(Psi,h,sigmaRelay,N,Ns);
        //cout<<"Psi Matrix: "<<endl;
        //cout<<Psi<<endl;
    
        //2a) generate small delta matrix: D
        MatrixXcd D(size2,UEs*UEs);
        D.setZero();
        generateD(D,h,g,N,Ns);

        //2b) generate PHI MATRIX from D
        MatrixXcd Phi (size2,UEs);
        for (int k =0 ; k<UEs ; k++)
        {
            for (int i =0;i<size2; i++)
            {
                Phi(i,k) =  D(i,k*UEs + k);
            }
        }
        //cout<<"Phi Matrix: "<<endl;
        //cout<<Phi<<endl;
        
        //3) generate DELTA MATRIX
        MatrixXcd Delta (size2,size2);
        Delta.setZero();
        generateDelta(Delta,D,size2);
        
        //4) generate the G Matrix
        MatrixXcd G(size2,size2);
        G.setZero();
        generateG(G, g, sigmaRelay, N, Ns, size2);
        
        //5) generate the Theta Matrix: THETA = G + DELTA
        MatrixXcd Theta(size2,size2);
        Theta = G + Delta;
        //cout<<"Theta Matrix: "<<endl;
        //cout<<Theta<<endl;
        
        
        /* BEGIN THE OPTIMIZATION PROCESS THROUGH THE ALTERNATING DIRECTION
         
         METHOD OF MULTIPLIERS (ADMM) ALGORITHM
        
         */
        Phi = Psi.pow(-0.5)*Phi;
    
        std :: clock_t c_start = std :: clock ();
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
        cout<<"Solving..."<<endl;
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
            eps_pri = sqrt(size2)*eps_abs + eps_rel*max(w.norm(),theta.norm()); //primal tolerance threshold
            eps_dual = sqrt(size2)*eps_abs + eps_rel*mu.norm();  //dual torelance threshold
            
            
            //compute the objective function at each iteration
            double objective_function = 0;
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
            
            MatrixXcd temp(1,1);
            temp = w.adjoint()*Psi.pow(-0.5)*Theta*Psi.pow(-0.5)*w;
            objective_function = temp.real().trace() + sumMixedNorm;
            
            //cout<<"Objective function: "<<objective_function<<". Primal Residual: "<<rk.norm()<<". Dual Residual: "<<sk.norm()<<"."<<endl;
            count++;
        }
        //end ADMM algorithm
        std :: clock_t c_end = std :: clock ();
        double elapsed = 1000.0*( c_end - c_start )/CLOCKS_PER_SEC;
        /* END OF OPTIMIZATION ALGORITHM */
        
        /*   DISPLAY RESULTS    */
        cout<<"--------------------------------------------"<<endl;
        cout<<"Status: Solved."<<endl;
        //cout<<"Optimal Value: "<<objective_function<<endl;
        cout<<"Primal Residual: "<<rk.norm()<<endl;
        cout<<"Dual Residual: "<<sk.norm()<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<"The solution vector is: "<<endl;
        cout<<theta<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<"ADMM converged in "<<count<<" iterations."<<endl;
        cout <<"CPU time used : "<<setprecision (2)<<elapsed<< " ms."<<endl;
        cout<<"--------------------------------------------"<<endl;

        //cout<<w<<endl;

        //compute expected SINR
        //create indexing matrices
        MatrixXi idxN(2,relays);
        MatrixXi idxNs(2,relays);
        
        for(int l=0;l<relays;l++)
        {
            idxN(0,l) = vectorSum(N, 0, l-1); //start index
            idxN(1,l) = vectorSum(N, 0, l)-1;   //end index
            
            idxNs(0,l) = vectorSum(Ns,0,l-1);   //start index
            idxNs(1,l) = vectorSum(Ns,0,l)-1;      //end index
        }
        //
        
        double SINR_k[UEs];
        double E_SINR=0;
        double forwardedNoise;
        MatrixXcd temp(1,1);
        MatrixXcd temp2(1,1);
        double desired;
        double MUI;
        
        for(int k=0;k<UEs;k++)
        {
            forwardedNoise=0;
            temp.setZero();
            for (int l=0;l<relays;l++)
            {
                //get the AF-matrix B_l
                MatrixXcd Psi_l(Ns.at(l),Ns.at(l));
                Psi_l.setZero();
                MatrixXcd theta_l(Ns.at(l),1);
                slice(Psi,idxNs(0,l),idxNs(1,l),idxNs(0,l),idxNs(1,l),Psi_l);
                slice(theta,idxNs(0,l),idxNs(1,l),0,0,theta_l);
                MatrixXcd b_l;
                b_l = Psi_l.pow(-0.5)*theta_l;
                Map<MatrixXcd> B_l(b_l.data(), N.at(l),N.at(l));
                
                
                //get the relevant channels
                MatrixXcd g_kl(N.at(l),1);
                MatrixXcd h_lk(N.at(l),1);
                slice(h,idxN(0,l),idxN(1,l),k,k,h_lk);
                slice(g,idxN(0,l),idxN(1,l),k,k,g_kl);
                
                //accumulate desired signal
                temp = temp + g_kl.adjoint()*B_l*h_lk;
                //compute forward relay-noise power
                forwardedNoise = forwardedNoise + 2*sigmaRelay.at(l)*sigmaRelay.at(l)*(g_kl.adjoint()*B_l).squaredNorm();
            }
            
            //compute desired signal's power
            desired=sigmaUE*sigmaUE*temp.squaredNorm();
    
            //compute Multi-User Interference (MUI)
            MUI = 0;
            for(int q=0;q<UEs;q++)
            {
                if(q==k) {continue;}
                else
                {
                    temp2.setZero();
                    //find sum of interference terms from a given interferer 'q'
                    //there are 'L' such copies of this component, one from each relay
                    for (int l=0;l<relays;l++)
                    {
                        //get the AF-matrix B_l
                        MatrixXcd Psi_l(Ns.at(l),Ns.at(l));
                        Psi_l.setZero();
                        MatrixXcd theta_l(Ns.at(l),1);
                        slice(Psi,idxNs(0,l),idxNs(1,l),idxNs(0,l),idxNs(1,l),Psi_l);
                        slice(theta,idxNs(0,l),idxNs(1,l),0,0,theta_l);
                        MatrixXcd b_l;
                        b_l = Psi_l.pow(-0.5)*theta_l;
                        Map<MatrixXcd> B_l(b_l.data(), N.at(l),N.at(l));
                        
                        //get the relevant channels
                        MatrixXcd g_kl(N.at(l),1);
                        MatrixXcd h_lq(N.at(l),1);
                        slice(h,idxN(0,l),idxN(1,l),q,q,h_lq);
                        slice(g,idxN(0,l),idxN(1,l),k,k,g_kl);
                        
                        temp2 = temp2 + g_kl.adjoint()*B_l*h_lq;
                    }
                    MUI  = MUI + sigmaUE*sigmaUE*temp2.squaredNorm();
                }
            }
            
            SINR_k[k] = desired/(MUI + forwardedNoise + sigmaDUE*sigmaDUE);
        }
        cout<<"--------------------------------------------"<<endl;
        for(int i=0;i<UEs;i++)
        {
            E_SINR = E_SINR + SINR_k[i];
            cout<<"SINR @ destination UE "<<i+1<<": "<<10*log10(SINR_k[i])<<" dB."<<endl;
        }
        E_SINR = E_SINR/UEs;
        cout<<"Average SINR: "<<E_SINR<<" dB."<<endl;
        cout<<"--------------------------------------------"<<endl;

        
    } //end Monte Carlo for-loop
    
    
    /* END OF MONTE CARLO SIMULATION */
    
    return 0;
}


