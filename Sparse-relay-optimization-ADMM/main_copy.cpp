//
//  main.cpp
//  Sparse-relay-optimization-ADMM
//
//  Created by Ayoub Saab  on 2019-07-08.
//  Copyright © 2019 Ayoub Saab. All rights reserved.
//

#include "mimo.h"

//preprocessor directives
#include "functions_copy.h"
#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>
using namespace std;


void test(){
        // //compute expected SINR
        // //create indexing matrices
        // MatrixXi idxN(2,relays);
        // MatrixXi idxNs(2,relays);
        
        // for(int l=0;l<relays;l++){
        //     idxN(0,l) = vectorSum(N, 0, l-1); //start index
        //     idxN(1,l) = vectorSum(N, 0, l)-1;   //end index
            
        //     idxNs(0,l) = vectorSum(Ns,0,l-1);   //start index
        //     idxNs(1,l) = vectorSum(Ns,0,l)-1;      //end index
        // }
        // //
        
        //   MatrixXcd temp(1,1);
        //   MatrixXcd temp2(1,1);

/* THIS WAS COMMENTED OUT

//        double SINR_k[UEs];
//        double E_SINR=0;
//        double forwardedNoise;
//        double desired;
//        double MUI;
//        
//        for(int k=0;k<UEs;k++)
//        {
//            forwardedNoise=0;
//            temp.setZero();
//            for (int l=0;l<relays;l++)
//            {
//                //get the AF-matrix B_l
//                MatrixXcd Psi_l(Ns.at(l),Ns.at(l));
//                Psi_l.setZero();
//                MatrixXcd theta_l(Ns.at(l),1);
//                slice(Psi,idxNs(0,l),idxNs(1,l),idxNs(0,l),idxNs(1,l),Psi_l);
//                slice(theta,idxNs(0,l),idxNs(1,l),0,0,theta_l);
//                MatrixXcd b_l;
//                b_l = Psi_l.pow(-0.5)*theta_l;
//                Map<MatrixXcd> B_l(b_l.data(), N.at(l),N.at(l));
//                
//                
//                //get the relevant channels
//                MatrixXcd g_kl(N.at(l),1);
//                MatrixXcd h_lk(N.at(l),1);
//                slice(h,idxN(0,l),idxN(1,l),k,k,h_lk);
//                slice(g,idxN(0,l),idxN(1,l),k,k,g_kl);
//                
//                //accumulate desired signal
//                temp = temp + g_kl.adjoint()*B_l*h_lk;
//                
//                //compute forward relay-noise power
//                forwardedNoise = forwardedNoise + 2*sigmaRelay.at(l)*sigmaRelay.at(l)*(g_kl.adjoint()*B_l).squaredNorm();
//            }
//            
//            //compute desired signal's power
//            desired=sigmaUE*sigmaUE*temp.squaredNorm();
//    
//            //compute Multi-User Interference (MUI)
//            MUI = 0;
//            for(int q=0;q<UEs;q++)
//            {
//                if(q==k) {continue;}
//                else
//                {
//                    temp2.setZero();
//                    //find sum of interference terms from a given interferer 'q'
//                    //there are 'L' such copies of this component, one from each relay
//                    for (int l=0;l<relays;l++)
//                    {
//                        //get the AF-matrix B_l
//                        MatrixXcd Psi_l(Ns.at(l),Ns.at(l));
//                        Psi_l.setZero();
//                        MatrixXcd theta_l(Ns.at(l),1);
//                        slice(Psi,idxNs(0,l),idxNs(1,l),idxNs(0,l),idxNs(1,l),Psi_l);
//                        slice(theta,idxNs(0,l),idxNs(1,l),0,0,theta_l);
//                        MatrixXcd b_l;
//                        b_l = Psi_l.pow(-0.5)*theta_l;
//                        Map<MatrixXcd> B_l(b_l.data(), N.at(l),N.at(l));
//                        
//                        //get the relevant channels
//                        MatrixXcd g_kl(N.at(l),1);
//                        MatrixXcd h_lq(N.at(l),1);
//                        slice(h,idxN(0,l),idxN(1,l),q,q,h_lq);
//                        slice(g,idxN(0,l),idxN(1,l),k,k,g_kl);
//                        
//                        temp2 = temp2 + g_kl.adjoint()*B_l*h_lq;
//                    }
//                    MUI  = MUI + sigmaUE*sigmaUE*temp2.squaredNorm();
//                }
//            }
//            
//            SINR_k[k] = desired/(MUI + forwardedNoise + sigmaDUE*sigmaDUE);
//        }
//        cout<<"--------------------------------------------"<<endl;
//        for(int i=0;i<UEs;i++)
//        {
//            E_SINR = E_SINR + SINR_k[i];
//            cout<<"SINR @ destination UE "<<i+1<<": "<<10*log10(SINR_k[i])<<" dB."<<endl;
//        }
//        E_SINR = E_SINR/UEs;
//        cout<<"Average SINR: "<<10*log10(E_SINR)<<" dB."<<endl;
//        cout<<"--------------------------------------------"<<endl;
//        cout<<"--------------------------------------------"<<endl<<endl;


    THIS WAS COMMENTED OUT */

// /* SIMULATE TRANSMISSION-RECEPTION    */
// //
//         cout<<"Simulating Transmission..."<<endl<<endl;
//         std::default_random_engine generator(time(NULL));
//         std::bernoulli_distribution distribution(0.5);
//         MatrixXcd sent(UEs,1);
//         MatrixXd SINR(UEs,1);
//         SINR.setZero();
//         //MatrixXcd rec(UEs,1);
//         bool x,y;
        
//         for(int t=0; t<coherenceTime ; t++){
//             //simulate transmition
//             for(int u = 0; u<UEs ; u++){
//                 x = distribution(generator);
//                 y = distribution (generator);
                
//                 if (x && y) {sent(u,0)=std::complex<double>(sigmaUE/sqrt(2),sigmaUE/sqrt(2));}
//                 if (x && !y) {sent(u,0)=std::complex<double>(sigmaUE/sqrt(2),-sigmaUE/sqrt(2));}
//                 if (!x && y) {sent(u,0)=std::complex<double>(-sigmaUE/sqrt(2),sigmaUE/sqrt(2));}
//                 if (!x && !y) {sent(u,0)=std::complex<double>(-sigmaUE/sqrt(2),-sigmaUE/sqrt(2));}
//                 cout<<"User "<<u+1<<" sent: "<<sent(u,0)<<"."<<endl;
//             }
//             cout<<endl;
            
//             //simulate reception
//             for (int u =0; u<UEs; u++){
//                 MatrixXcd forwardedNoise(1,1);
//                 forwardedNoise.setZero();
//                 temp.setZero();
//                 for (int l=0;l<relays;l++){
//                     //get the AF-matrix B_l
//                     MatrixXcd Psi_l(Ns.at(l),Ns.at(l));
//                     Psi_l.setZero();
//                     MatrixXcd theta_l(Ns.at(l),1);
//                     slice(Psi,idxNs(0,l),idxNs(1,l),idxNs(0,l),idxNs(1,l),Psi_l);
//                     slice(theta,idxNs(0,l),idxNs(1,l),0,0,theta_l);
//                     MatrixXcd b_l;
//                     b_l = Psi_l.pow(-0.5)*theta_l;
//                     Map<MatrixXcd> B_l(b_l.data(), N.at(l),N.at(l));
                    
                    
//                     //get the relevant channels
//                     MatrixXcd g_kl(N.at(l),1);
//                     MatrixXcd h_lk(N.at(l),1);
//                     slice(h,idxN(0,l),idxN(1,l),u,u,h_lk);
//                     slice(g,idxN(0,l),idxN(1,l),u,u,g_kl);
                    
//                     //accumulate desired signal
//                     temp = temp + g_kl.adjoint()*B_l*h_lk*sent(u,0);
                    
//                     //generate noise-at-relay
//                     MatrixXcd n_l(N.at(l),1);
//                     std::normal_distribution<double> distribution(0,sigmaRelay.at(l));
//                     //accumulate forwarded-noise
//                     for(int n =0;n<N.at(l);n++){
//                         n_l(n,0)=std::complex<double>(distribution(generator),distribution(generator));
//                     }
//                     forwardedNoise = forwardedNoise + g_kl.adjoint()*B_l*n_l;
//                 }
                
//                 //compute Multi-User Interference (MUI)
//                 temp2.setZero();
//                 for(int q=0;q<UEs;q++){
//                     if(q==u) {
//                         continue;
//                     }
//                     else{
//                         //find sum of interference terms from a given interferer 'q'
//                         //there are 'L' such copies of this component, one from each relay
//                         for (int l=0;l<relays;l++){
//                             //get the AF-matrix B_l
//                             MatrixXcd Psi_l(Ns.at(l),Ns.at(l));
//                             Psi_l.setZero();
//                             MatrixXcd theta_l(Ns.at(l),1);
//                             slice(Psi,idxNs(0,l),idxNs(1,l),idxNs(0,l),idxNs(1,l),Psi_l);
//                             slice(theta,idxNs(0,l),idxNs(1,l),0,0,theta_l);
//                             MatrixXcd b_l;
//                             b_l = Psi_l.pow(-0.5)*theta_l;
//                             Map<MatrixXcd> B_l(b_l.data(), N.at(l),N.at(l));
                            
//                             //get the relevant channels
//                             MatrixXcd g_kl(N.at(l),1);
//                             MatrixXcd h_lq(N.at(l),1);
//                             slice(h,idxN(0,l),idxN(1,l),q,q,h_lq);
//                             slice(g,idxN(0,l),idxN(1,l),u,u,g_kl);
            
//                             temp2 = temp2 + g_kl.adjoint()*B_l*h_lq*sent(q,0);
//                         }
//                     }
//                 }
//                 //generate local-noise at destination-UE
//                 std::normal_distribution<double> distribution(0,sqrt(0.5)*sigmaDUE);
//                 MatrixXcd n_d (1,1);
//                 n_d(0,0) = std::complex<double>(distribution(generator),distribution(generator));

//                 cout<<"User "<<u+1<<" received:"<<endl;
//                 cout<<"Desired signal: "<<temp<<", Interference: "<<temp2<<", Forwarded noise: ";
//                 cout<<forwardedNoise<<", Local noise: "<<n_d<<endl;
//                 cout<<"SINR: "<<10*log10(temp.squaredNorm()/((temp2+forwardedNoise).squaredNorm()+n_d.squaredNorm()))<<" dB."<<endl<<endl;
                
//                 SINR(u,0)=SINR(u,0)+temp.squaredNorm()/((temp2+forwardedNoise).squaredNorm()+n_d.squaredNorm());
//             }
//             //end reception
//             cout<<"--------------------------------------------"<<endl;
//         }
//         //end coherence time loop
        
//         for (int u=0;u<UEs;u++){
//             SINR(u,0) = 10*log10(SINR(u,0)/coherenceTime);
//             cout<<"Average SINR @User"<<u+1<<" over "<<coherenceTime<<" transmissions: "<<SINR(u,0)<<" dB."<<endl;
//         }

}

int main()
{
    ofstream outData;
    mimoNetwork_t mimoObject;
    string fileName = "./data/sample_";
    string sep = " | ";

    
    //outData<<"system parameters: UES, RELAYS, ...,"

    string headerLine = "iter" + sep + "Obj." + sep + "Residual (primal)" + sep + "Residual (dual)" + sep;
    /* BEGIN MONTE CARLO SIMULATION */
    for (int sim = 0 ; sim < monteCarlo ; sim++){
        
        outData.open(fileName + std::to_string(sim)+".txt", ios::out);
        outData<<headerLine<<endl; //display header
        outData<<"-----------------------------------------------"<<endl;

        MatrixXcd solution_vector = mimoObject.solve(false,outData);
        //mimoObject.simulateTxRx(solution_vector);
        
        outData.close();
    } //end Monte Carlo for-loop

    /* END OF MONTE CARLO SIMULATION */ 
    return 0;
}


