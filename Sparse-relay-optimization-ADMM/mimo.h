
//preprocessor directives
#include "parameters.h"
#include "functions.h"

#include <complex>
#include <iostream>
#include <iomanip>

//necessary library for Matrix manipulations
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/MatrixFunctions"
using Eigen::MatrixXd;

using namespace Eigen;
using namespace std;

class mimoNetwork_t{
    private:
        //fixed parameters
        unsigned int K; //UEs
        unsigned int L; //number of relays
        //fixed parameters
        vector<unsigned int> N; //relay antennas
        vector<unsigned int> Ns; //square of relay antennas
        vector<unsigned int> P; //power budgets at each relay
        vector<double> sigmaRelay; //std-dev at each relay station
        MatrixXcd c; //distortionless constraints
        
        //these change for each Monte Carlo simulation
        MatrixXcd h,g; //channels
        //relevant matrices
        MatrixXcd Psi; 
        MatrixXcd Phi;
        MatrixXcd Theta;


    public:
    mimoNetwork_t(){//default constructor
        K = UEs;
        L = relays;
        this->initializeSystem(); //initialize System
    }

    private:
    void initializeSystem(){//initialize the system parameters
        for(int i=0; i<L; i++){
            N.push_back(antenna[i]);
            Ns.push_back(N.at(i)*N.at(i));
        }

        for (int i=0; i<L; i++){
            P.push_back(power[i]); //init. relay-power budgets
        }

        //initialize the distortionless constraints 
        c.setOnes();
        c=c*distortionlessConstraint;

        //initalize standard deviation @ each relay-station
        for(int i=0; i<L; i++){
            sigmaRelay.push_back(sqrt(0.5)*sigmaR);
        }
    }
    
    private:
    bool generateChannels(const float sigmaChannel){
        std::normal_distribution<double> distribution(0,sigmaChannel);
        std::default_random_engine generator(time(NULL));
    
        //determine the size
        int temp=0;
        for (int i=0; i<L;i++){
            temp = temp+N[i];
        }
    
        for (int k=0; k<K; k++){
            for( int n=0;n<temp;n++){
                h(n,k) = std::complex<double>(distribution(generator),distribution(generator));    //backward channel
                g(n,k) = std::complex<double>(distribution(generator),distribution(generator)); //forward channel: complex conjugate
            }
        }
        return true;
    }

    public:
    void showChannels(){
        cout<<"--------------------------------------------"<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<"Backward channel H coefficients: "<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<this->h<<endl<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<"Forward channel G coefficients: "<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<this->g<<endl<<endl;
        cout<<"--------------------------------------------"<<endl;
        cout<<"--------------------------------------------"<<endl;
    }

    private:
    bool createChannels(){
        /* 
        
        GENERATE THE WIRELESS CHANNELS H (backward) and G (foreward) 
        
        */
        const int size = vectorSum(N,0,L-1); //sum_{i=1}^L N_l
        //generate channels
        h = MatrixXcd(size,K);
        g = MatrixXcd(size,K);
        this->generateChannels(sigmaChannel);
        this->showChannels();
        return true;
    }

    private:
    bool generatePsi(){
        //generate Psi_l, l=0,...,L-1
        for (int l=0;l<L;l++){
            MatrixXcd Psi_l (N[l]*N[l],N[l]*N[l]);
            MatrixXcd temp (N[l],K);
            MatrixXcd temp2(N[l],N[l]);
            MatrixXcd eye_l (N[l],N[l]);
            eye_l.setIdentity();
            
            int offset = vectorSum(N,0,l-1);
            
            //populate the temp matrix:
            for (int k=0; k<K; k++){
                for (int j=0;j<N[l];j++){
                    temp(j,k) = h(j+offset,k);
                }
            }
            
            temp2 = sigmaUE*sigmaUE*temp*temp.adjoint() + 2*sigmaRelay[l]*sigmaRelay[l]*eye_l;
            
            //compute the Kronecker product of temp2 and eye_l and store it in Psi_l
            Kroneckerproduct(temp2, eye_l , Psi_l);
            
            //populate the Psi matrix with Psi_l, for l=0,1,...,L-1 in a block diagonal fashion
            offset = vectorSum(Ns, 0, l-1);
            for (int i=0;i<N[l]*N[l]; i++){
                for(int j=0; j<N[l]*N[l]; j++){
                    Psi(i+offset,j+offset) = Psi_l(i,j);
                }
            }
        }
        return true;
    }

    private:
    void generateD(MatrixXcd &D){

        for (int l=0;l<L;l++){
            int offset2 = vectorSum(Ns, 0, l-1);
            MatrixXcd temp (N[l]*N[l],1);
            for (int k=0 ; k<K ; k++){
                //get g_kl
                MatrixXcd g_kl (N[l],1);
                int offset = vectorSum(N, 0, l-1);
                for (int n=0; n<N[l]; n++){
                    g_kl(n,0) = g(n+offset,k);
                }
                
                for (int j=0 ; j<K ; j++){
                    //get h_lj
                    MatrixXcd h_lj (N[l],1);
                    for(int m=0 ; m<N[l] ; m++)
                    {
                        h_lj(m,0) = h(offset+m,j);
                    }
                    
                    //compute the kronecker prodcut of h_lj* and g_kl and store it in temp
                    Kroneckerproduct(h_lj.conjugate(), g_kl, temp);
                    
                    //store all the results in D
                    for(int z=0 ; z<N[l]*N[l] ; z++){
                        D(z+offset2,j+k*K) = temp(z,0);
                    }
                }
            }
        }
    }
    
    private:
    void generateG(MatrixXcd& G, const int size2){
        for (int k=0;k<UEs;k++) {
            MatrixXcd G_k(size2,size2);
            G_k.setZero();
            
            for (int l=0;l<relays; l++){
                MatrixXcd g_kl(antenna[l],1);
                int offset = vectorSum(N,0,l-1);
                //get g_kl
                for (int m=0;m<antenna[l];m++){
                    g_kl(m,0) = g(m+offset,k);
                }

                //get G_kl
                MatrixXcd G_kl(N[l]*N[l],N[l]*N[l]);
                G_kl.setZero();
                
                for (int j=0;j<N[l];j++){
                    MatrixXcd eye_j (N[l],1);
                    eye_j.setZero();
                    eye_j(j,0) = 1;
                    MatrixXcd kron (N[l]*N[l],1);
                    Kroneckerproduct(eye_j, g_kl, kron);
                    G_kl = G_kl + kron*kron.adjoint() ;
                }
                G_kl = G_kl*2*sigmaRelay[l]*sigmaRelay[l];
                
                //G_k = blkdiag{G_k1, G_k2, ... , G_kL}
                int offset2 = vectorSum(Ns, 0, l-1);
                for (int n=0; n<N[l]*N[l]; n++){
                    for (int m=0; m<N[l]*N[l]; m++)
                        G_k(n+offset2,m+offset2) = G_kl(n,m);
                }
            }
            G = G + G_k; //G =
        }
    }
    
    private:
    void generateDelta(MatrixXcd& Delta, const MatrixXcd D, const int size2){
        MatrixXcd temp(size2,size2);
        
        for (int k=0; k<K; k++){
            MatrixXcd delta_k_j(size2,1);
            
            temp.setZero();
            for (int j=0; j<K; j++){
                if(j!=k){
                    //populate
                    for(int i =0;i<size2;i++){
                        delta_k_j(i,0) = D(i,j+K*k);
                    }
                    //add the outer product of delta_kj
                    temp = temp + delta_k_j*delta_k_j.adjoint(); // Delta = sum_{k=1}^{UEs} Delta_k
                }
            }
            temp = temp*sigmaUE*sigmaUE;
            Delta = Delta + temp;
        }
    }
    
    private:
    double ADMM(const int size2){

        Phi = Psi.pow(-0.5)*Phi;
        std :: clock_t c_start = std :: clock ();
        
        //begin ADMM algorithm
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
        //algorithm execution
        cout<<"Solving..."<<endl;
        while( rk.norm() > eps_pri || sk.norm() > eps_dual ){
            //ADMM step1, update w
            Q_tilde = Phi.adjoint()*Qinverse*Phi;
            nu = Q_tilde.inverse()*(c+Phi.adjoint()*Qinverse*mu-(rho/2.0)*Phi.adjoint()*Qinverse*theta);
            cout<<"nu is ok"<<endl;
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
                MatrixXcd temp (antenna[l]*antenna[l],1);
                int offset = vectorSum(Ns, 0, l-1);
                for (int j=0;j<antenna[l]*antenna[l];j++){
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
        cout<<"--------------------------------------------"<<endl;
        
        return elapsed;
    }

    public:
    void solve(){
        
        //STEP 1: generate new channels
        this->createChannels(); //generate new channels
        
        //STEP2: 
        /* GENERATE THE REQUIRED MATRICES FOR THE OPTIMIZATION FUNCTION:
         
         1)PSI, 
         2)PHI 
         3)DELTA, 
         4)G,
         5)THETA=DELTA+G
         
         */
        const int size2 = vectorSum(Ns,0,L-1); //sum_{i=1}^L N_l^2

        //1) generate Psi: Matrix to express relay-power constraints
        Psi = MatrixXcd(size2,size2);
        Psi.setZero();
        this->generatePsi();
        // cout<<"--------------------------------------------"<<endl;
        // cout<<"--------------------------------------------"<<endl;
        // cout<<"Psi Matrix: "<<endl;
        // cout<<Psi<<endl;
    
        //2a) generate small delta matrix: D
        MatrixXcd D(size2,K*K);
        D.setZero();
        this->generateD(D);

        //2b) generate PHI MATRIX from D
        Phi = MatrixXcd(size2,K);
        for (int k =0 ; k<K ; k++){
            for (int i =0;i<size2; i++){
                Phi(i,k) =  D(i,k*K + k);
            }
        }
        // cout<<"--------------------------------------------"<<endl;
        // cout<<"--------------------------------------------"<<endl;
        // cout<<"Phi Matrix: "<<endl;
        // cout<<Phi<<endl;
        
        //3) generate DELTA MATRIX
        MatrixXcd Delta (size2,size2);
        Delta.setZero();
        this->generateDelta(Delta,D,size2);
        
        //4) generate the G Matrix
        MatrixXcd G(size2,size2);
        G.setZero();
        this->generateG(G, size2);
        
        //5) generate the Theta Matrix: THETA = G + DELTA
        Theta = MatrixXcd(size2,size2);
        Theta = G + Delta;
        // cout<<"--------------------------------------------"<<endl;
        // cout<<"--------------------------------------------"<<endl;
        // cout<<"Theta Matrix: "<<endl;
        // cout<<Theta<<endl;


        //STEP3: SOLVE
        /* 
            BEGIN THE OPTIMIZATION PROCESS THROUGH THE ALTERNATING DIRECTION
    
            METHOD OF MULTIPLIERS (ADMM) ALGORITHM
         */
        
        //this->ADMM(size2);
    }
};

