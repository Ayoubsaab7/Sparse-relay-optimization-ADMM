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
#include <vector>
#include <complex>
#include <random>
#include <time.h>

//necessary library for Matrix manipulations
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/MatrixFunctions"
using Eigen::MatrixXd;
using namespace Eigen;

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


void slice (const MatrixXcd A, int rowA, int rowB, int colA,int colB,MatrixXcd &A_sliced){
    int idx_i=0;
    //test
    //cout<<"rowA: "<<rowA<<" rowB: "<<rowB<<endl;
    //cout<<"colA: "<<colA<<" colB: "<<colB<<endl;
    for (int i=rowA;i<=rowB ; i++){
        int idx_j=0;
        for(int j=colA;j<=colB;j++){
            //cout<<"(i,j) ("<<i<<","<<j<<");"<<endl;
            //cout<<"A(i,j) ="<<A(i,j)<<endl;
            A_sliced(idx_i,idx_j) = A(i,j);
            idx_j++;
        }
        idx_i++;
        //cout<<"incremented the row"<<endl;
    }
}

int vectorSum(const vector<unsigned int> & V, int startIndex,  int endIndex){
    int sum = 0;
    if(startIndex<=endIndex){
        for(int i = startIndex; i<=endIndex ; i++){
            sum = sum + V[i];
        }
    }
    return sum;
}

// void initializeSystem(vector<unsigned int> & N, vector<unsigned int> & Ns, vector<unsigned int> &P, MatrixXcd & c, vector<double>& sigmaRelay){

//     //initialize the system parameters
//     for(int i=0;i<N.size();i++){
//         N.at(i) = antenna[i];
//         Ns.at(i) = N.at(i)*N.at(i);
//     }

//     for (int i=0;i<P.size();i++){
//         P.at(i)=power[i];   //initialize the relay-power budgets
//     }

//     //initialize the distortionless constraints 
//     c.setOnes();
//     c=c*distortionlessConstraint; 

//     //initalize standard deviation @ each relay-station
//     for(int i=0;i<sigmaRelay.size();i++){
//         sigmaRelay.at(i) = sqrt(0.5)*sigmaR;
//     }
// }

void Kroneckerproduct(const MatrixXcd  A, const MatrixXcd B, MatrixXcd &C) {
    long rowA = A.rows();
    long rowB = B.rows();
    long colA = A.cols();
    long colB = B.cols();
    
    long startRow,startCol;
    
    for (int i=0 ; i<rowA ; i++){
        for(int j=0;j<colA ; j++){
            startRow = i*rowB;
            startCol = j*colB;
            for (int k =0 ;k<rowB;k++ ){
                for(int l=0 ; l<colB ;l++){
                    C(startRow+k,startCol+l) = A(i,j)*B(k,l);
                }
            }
        }
    }
}

// void generateChannels(MatrixXcd& h, MatrixXcd& g, const float sigmaChannel){
//     std::normal_distribution<double> distribution(0,sigmaChannel);
//     std::default_random_engine generator(time(NULL));
    
//     //determine the size
//     int temp=0;
//     for (int i=0; i<relays;i++){
//         temp = temp+antenna[i];
//     }
    
//     for (int k=0;k<UEs;k++){
//         for( int n=0;n<temp;n++){
//             h(n,k) = std::complex<double>(distribution(generator),distribution(generator));    //backward channel
//             g(n,k) = std::complex<double>(distribution(generator),distribution(generator)); //forward channel: complex conjugate
//         }
//     }
// }

// void generatePsi(MatrixXcd & Psi, const MatrixXcd h, const vector<double> sigmaRelay, const vector <unsigned int> N, const vector<unsigned int> Ns){
//     //generate Psi_l, l=0,...,L-1
//     for (int l=0;l<relays;l++){
//         MatrixXcd Psi_l (antenna[l]*antenna[l],antenna[l]*antenna[l]);
//         MatrixXcd temp (antenna[l],UEs);
//         MatrixXcd temp2(antenna[l],antenna[l]);
//         MatrixXcd eye_l (antenna[l],antenna[l]);
//         eye_l.setIdentity();
        
//         int offset = vectorSum(N,0,l-1);
        
//         //populate the temp matrix:
//         for (int k=0;k<UEs;k++){
//             for (int j=0;j<antenna[l];j++){
//                 temp(j,k) = h(j+offset,k);
//             }
//         }
        
//         temp2 = sigmaUE*sigmaUE*temp*temp.adjoint() + 2*sigmaRelay[l]*sigmaRelay[l]*eye_l;
        
//         //compute the Kronecker product of temp2 and eye_l and store it in Psi_l
//         Kroneckerproduct(temp2, eye_l , Psi_l);
        
//         //populate the Psi matrix with Psi_l, for l=0,1,...,L-1 in a block diagonal fashion
//         offset = vectorSum(Ns, 0, l-1);
//         for (int i=0;i<antenna[l]*antenna[l]; i++){
//             for(int j=0; j<antenna[l]*antenna[l]; j++){
//                 Psi(i+offset,j+offset) = Psi_l(i,j);
//             }
//         }
//     }
// }

// void generateD(MatrixXcd &D,const MatrixXcd h, const MatrixXcd g, const vector <unsigned int> N, const vector <unsigned int> Ns){
    
//     for (int l=0;l<relays;l++){
//         int offset2 = vectorSum(Ns, 0, l-1);
//         MatrixXcd temp (antenna[l]*antenna[l],1);
//         for (int k=0 ; k<UEs ; k++){
//             //get g_kl
//             MatrixXcd g_kl (antenna[l],1);
//             int offset = vectorSum(N, 0, l-1);
//             for (int n=0;n<antenna[l];n++){
//                 g_kl(n,0) = g(n+offset,k);
//             }
            
//             for (int j=0 ; j<UEs ; j++){
//                 //get h_lj
//                 MatrixXcd h_lj (antenna[l],1);
//                 for(int m=0 ; m<antenna[l] ; m++)
//                 {
//                     h_lj(m,0) = h(offset+m,j);
//                 }
                
//                 //compute the kronecker prodcut of h_lj* and g_kl and store it in temp
//                 Kroneckerproduct(h_lj.conjugate(), g_kl, temp);
                
//                 //store all the results in D
//                 for(int z=0 ; z<antenna[l]*antenna[l] ; z++){
//                     D(z+offset2,j+k*UEs) = temp(z,0);
//                 }
//             }
//         }
//     }
// }

// void generateDelta(MatrixXcd& Delta, const MatrixXcd D, const int size2){
//     MatrixXcd temp(size2,size2);
    
//     for (int k=0;k<UEs;k++){
//         MatrixXcd delta_k_j(size2,1);
        
//         temp.setZero();
//         for (int j=0;j<UEs;j++){
//             if(j!=k){
//                 //populate
//                 for(int i =0;i<size2;i++){
//                     delta_k_j(i,0) = D(i,j+UEs*k);
//                 }
//                 //add the outer product of delta_kj
//                 temp = temp + delta_k_j*delta_k_j.adjoint(); // Delta = sum_{k=1}^{UEs} Delta_k
//             }
//         }
//         temp = temp*sigmaUE*sigmaUE;
//         Delta = Delta + temp;
//     }
// }

// void generateG(MatrixXcd& G, const MatrixXcd g,  const vector<double> sigmaRelay,const vector <unsigned int> N, const vector <unsigned int> Ns, const int size2){
//     for (int k=0;k<UEs;k++) {
//         MatrixXcd G_k(size2,size2);
//         G_k.setZero();
        
//         for (int l=0;l<relays; l++){
//             MatrixXcd g_kl(antenna[l],1);
//             int offset = vectorSum(N,0,l-1);
//             //get g_kl
//             for (int m=0;m<antenna[l];m++){
//                 g_kl(m,0) = g(m+offset,k);
//             }

//             //get G_kl
//             MatrixXcd G_kl(antenna[l]*antenna[l],antenna[l]*antenna[l]);
//             G_kl.setZero();
            
//             for (int j=0;j<antenna[l];j++){
//                 MatrixXcd eye_j (antenna[l],1);
//                 eye_j.setZero();
//                 eye_j(j,0) = 1;
//                 MatrixXcd kron (antenna[l]*antenna[l],1);
//                 Kroneckerproduct(eye_j, g_kl, kron);
//                 G_kl = G_kl + kron*kron.adjoint() ;
//             }
//             G_kl = G_kl*2*sigmaRelay[l]*sigmaRelay[l];
            
//             //G_k = blkdiag{G_k1, G_k2, ... , G_kL}
//             int offset2 = vectorSum(Ns, 0, l-1);
//             for (int n=0;n<antenna[l]*antenna[l];n++){
//                 for (int m=0;m<antenna[l]*antenna[l];m++)
//                     G_k(n+offset2,m+offset2) = G_kl(n,m);
//             }
//         }
//         G = G + G_k; //G =
//     }
// }

#endif /* functions_h */
