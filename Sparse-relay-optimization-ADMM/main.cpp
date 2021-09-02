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
#include <iomanip>
#include <ctime>

#include "mimo.h"

using namespace std;

int main()
{
    bool display = false;
    bool simulateTxRx = false;

    //declare system parameters
    vector<unsigned int> N(relays);
    vector <unsigned int> Ns(relays);    
    vector<unsigned int> P(relays); //relay-power budgets
    MatrixXcd c(UEs,1); //distortionless constraints vector
    vector<double> sigmaRelay (relays); //stdev @ different stations

    //initialize system parameters
    initializeSystem(N,Ns,P,c,sigmaRelay);
    const int size = vectorSum(N,0,relays-1); //sum_{i=1}^L N_l
    const int size2 = vectorSum(Ns,0,relays-1); //sum_{i=1}^L N_l^2
    
    /* BEGIN MONTE CARLO SIMULATION */
    for (int sim = 0 ; sim<monteCarlo ; sim++){
        /* 
        
        GENERATE THE WIRELESS CHANNELS H (backward) and G (foreward) 
        
        */

        //generate channels
        MatrixXcd h(size,UEs);
        MatrixXcd g(size,UEs);
        generateChannels(h,g,sigmaChannel);
        if (display){
            cout<<"--------------------------------------------"<<endl;
            cout<<"--------------------------------------------"<<endl;
            cout<<"Backward channel H coefficients: "<<endl;
            cout<<"--------------------------------------------"<<endl;
            cout<<h<<endl<<endl;
            cout<<"--------------------------------------------"<<endl;
            cout<<"Forward channel G coefficients: "<<endl;
            cout<<"--------------------------------------------"<<endl;
            cout<<g<<endl<<endl;
            cout<<"--------------------------------------------"<<endl;
            cout<<"--------------------------------------------"<<endl;
        }
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
    
        //2a) generate small delta matrix: D
        MatrixXcd D(size2,UEs*UEs);
        D.setZero();
        generateD(D,h,g,N,Ns);

        //2b) generate PHI MATRIX from D
        MatrixXcd Phi (size2,UEs);
        for (int k =0 ; k<UEs ; k++){
            for (int i =0;i<size2; i++){
                Phi(i,k) =  D(i,k*UEs + k);
            }
        }
        
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
        
        /* 
        
            BEGIN THE OPTIMIZATION PROCESS THROUGH THE ALTERNATING DIRECTION
         
            METHOD OF MULTIPLIERS (ADMM) ALGORITHM
        
         */

        Phi = Psi.pow(-0.5)*Phi;
        
        //begin ADMM algorithm
        std :: clock_t c_start = std :: clock ();

        MatrixXcd Q(size2,size2), Qinverse(size2,size2), Q_tilde(UEs,UEs);
        MatrixXcd I(size2,size2);
        I.setIdentity();
        MatrixXcd nu (UEs,1);
        MatrixXcd theta(size2,1), w(size2,1), mu(size2,1);
        MatrixXcd rk(size2,1),sk(size2,1),sum_sk(size2,1);

        double eps_pri, eps_dual;
        
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
        if(display){
            cout<<"Solving..."<<endl;
        }
        
        //algorithm execution
        Q_tilde = Phi.adjoint()*Qinverse*Phi;
        while( rk.norm() > eps_pri || sk.norm() > eps_dual ){
            
            //ADMM step1, update w
            nu = Q_tilde.inverse()*(c+Phi.adjoint()*Qinverse*mu-(rho/2.0)*Phi.adjoint()*Qinverse*theta);            
            w = Qinverse*(rho/2.0*theta-mu+Phi*nu);
            
            //ADMM step2, update theta by solving L parallel subproblems
            for(int l=0; l<relays; l++){
                //populate a_l
                MatrixXcd a_l (antenna[l]*antenna[l],1);
                int offset = vectorSum(Ns, 0, l-1);
                for(int i=0; i<antenna[l]*antenna[l]; i++){
                    a_l(i,0) = (rho/2.0)*w(i+offset,0) + mu(i+offset,0);
                }
                
                if(a_l.norm() - lambda[l] <=0){
                    //set theta_l to zeros
                    for (int j=0;j<antenna[l]*antenna[l];j++){
                        theta(j+offset,0) = 0;
                    }
                }
                
                else{
                    double eta_l;
                    eta_l=std::max(0.0,(a_l.norm()-lambda[l])/sqrt(P[l])-rho/2.0);
                    for (int i =0 ; i<antenna[l]*antenna[l]; i++){
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
            eps_pri = sqrt(size2)*eps_abs + std::max(w.norm(),theta.norm())*eps_rel; //primal tolerance threshold
            eps_dual = sqrt(size2)*eps_abs + eps_rel*mu.norm();  //dual torelance threshold
            
            //compute the objective function at each iteration
            double objective_function = 0;
            double sumMixedNorm = 0;
            for (int l=0 ; l<relays ; l++){
                MatrixXcd temp (N[l]*N[l],1);
                int offset = vectorSum(Ns, 0, l-1);
                for (int j=0; j<N[l]*N[l]; j++){
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
        
        if (display){
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
            cout<<"--------------------------------------------"<<endl;
        }

        if (simulateTxRx){
                
            //compute expected SINR
            //create indexing matrices
            MatrixXi idxN(2,relays);
            MatrixXi idxNs(2,relays);
            
            for(int l=0;l<relays;l++){
                idxN(0,l) = vectorSum(N, 0, l-1); //start index
                idxN(1,l) = vectorSum(N, 0, l)-1;   //end index
                
                idxNs(0,l) = vectorSum(Ns,0,l-1);   //start index
                idxNs(1,l) = vectorSum(Ns,0,l)-1;      //end index
            }
            //
            
            MatrixXcd temp(1,1);
            MatrixXcd temp2(1,1);
            
            /* SIMULATE TRANSMISSION-RECEPTION    */
            std::default_random_engine generator(time(NULL));
            std::bernoulli_distribution distribution(0.5);
            MatrixXcd sent(UEs,1);
            MatrixXd SINR(UEs,1);
            SINR.setZero();
            //MatrixXcd rec(UEs,1);
            bool x,y;
            
            for(int t=0; t<coherenceTime ; t++){
                //simulate transmition
                cout<<"Simulating Transmission..."<<endl<<endl;
                for(int u = 0; u<UEs ; u++){
                    x = distribution(generator);
                    y = distribution (generator);
                    
                    if (x && y) {sent(u,0)=std::complex<double>(sigmaUE/sqrt(2),sigmaUE/sqrt(2));}
                    if (x && !y) {sent(u,0)=std::complex<double>(sigmaUE/sqrt(2),-sigmaUE/sqrt(2));}
                    if (!x && y) {sent(u,0)=std::complex<double>(-sigmaUE/sqrt(2),sigmaUE/sqrt(2));}
                    if (!x && !y) {sent(u,0)=std::complex<double>(-sigmaUE/sqrt(2),-sigmaUE/sqrt(2));}
                    cout<<"User "<<u+1<<" sent: "<<sent(u,0)<<"."<<endl;
                }
                cout<<endl;
                
                //simulate reception
                cout<<"Simulating Reception..."<<endl<<endl;
                for (int u =0; u<UEs; u++){
                    MatrixXcd forwardedNoise(1,1);
                    forwardedNoise.setZero();
                    temp.setZero();
                    for (int l=0;l<relays;l++){
                        //get the AF-matrix B_l
                        MatrixXcd Psi_l(Ns[l],Ns[l]);
                        Psi_l.setZero();
                        MatrixXcd theta_l(Ns[l],1);
                        slice(Psi,idxNs(0,l),idxNs(1,l),idxNs(0,l),idxNs(1,l),Psi_l);
                        slice(theta,idxNs(0,l),idxNs(1,l),0,0,theta_l);
                        MatrixXcd b_l;
                        b_l = Psi_l.pow(-0.5)*theta_l;
                        Map<MatrixXcd> B_l(b_l.data(), N[l],N[l]);
                        
                        
                        //get the relevant channels
                        MatrixXcd g_kl(N[l],1);
                        MatrixXcd h_lk(N[l],1);
                        slice(h,idxN(0,l),idxN(1,l),u,u,h_lk);
                        slice(g,idxN(0,l),idxN(1,l),u,u,g_kl);
                        
                        //accumulate desired signal
                        temp = temp + g_kl.adjoint()*B_l*h_lk*sent(u,0);
                        
                        //generate noise-at-relay
                        MatrixXcd n_l(N[l],1);
                        std::normal_distribution<double> distribution(0,sigmaRelay[l]);
                        //accumulate forwarded-noise
                        for(int n =0;n<N.at(l);n++){
                            n_l(n,0)=std::complex<double>(distribution(generator),distribution(generator));
                        }
                        forwardedNoise = forwardedNoise + g_kl.adjoint()*B_l*n_l;
                    }
                    
                    //compute Multi-User Interference (MUI)
                    temp2.setZero();
                    for(int q=0;q<UEs;q++){
                        if(q==u) {
                            continue;
                        }
                        else{
                            //find sum of interference terms from a given interferer 'q'
                            //there are 'L' such copies of this component, one from each relay
                            for (int l=0;l<relays;l++){
                                //get the AF-matrix B_l
                                MatrixXcd Psi_l(Ns[l],Ns[l]);
                                Psi_l.setZero();
                                MatrixXcd theta_l(Ns[l],1);
                                slice(Psi,idxNs(0,l),idxNs(1,l),idxNs(0,l),idxNs(1,l),Psi_l);
                                slice(theta,idxNs(0,l),idxNs(1,l),0,0,theta_l);
                                MatrixXcd b_l;
                                b_l = Psi_l.pow(-0.5)*theta_l;
                                Map<MatrixXcd> B_l(b_l.data(), N[l],N[l]);
                                
                                //get the relevant channels
                                MatrixXcd g_kl(N[l],1);
                                MatrixXcd h_lq(N[l],1);
                                slice(h,idxN(0,l),idxN(1,l),q,q,h_lq);
                                slice(g,idxN(0,l),idxN(1,l),u,u,g_kl);
                
                                temp2 = temp2 + g_kl.adjoint()*B_l*h_lq*sent(q,0);
                            }
                        }
                    }
                    //generate local-noise at destination-UE
                    std::normal_distribution<double> distribution(0,sqrt(0.5)*sigmaDUE);
                    MatrixXcd n_d (1,1);
                    n_d(0,0) = std::complex<double>(distribution(generator),distribution(generator));

                    cout<<"User "<<u+1<<" received:"<<endl;
                    cout<<"Desired signal: "<<temp<<", Interference: "<<temp2<<", Forwarded noise: ";
                    cout<<forwardedNoise<<", Local noise: "<<n_d<<endl;
                    cout<<"SINR: "<<10*log10(temp.squaredNorm()/((temp2+forwardedNoise).squaredNorm()+n_d.squaredNorm()))<<" dB."<<endl<<endl;
                    
                    SINR(u,0)=SINR(u,0)+temp.squaredNorm()/((temp2+forwardedNoise).squaredNorm()+n_d.squaredNorm());
                }
                //end reception
                cout<<"--------------------------------------------"<<endl;
            }
            //end coherence time loop
            
            for (int u=0;u<UEs;u++){
                SINR(u,0) = 10*log10(SINR(u,0)/coherenceTime);
                cout<<"Average SINR @User"<<u+1<<" over "<<coherenceTime<<" transmissions: "<<SINR(u,0)<<" dB."<<endl;
            }
        }
    } //end Monte Carlo for-loop

    /* END OF MONTE CARLO SIMULATION */ 
    return 0;
}


