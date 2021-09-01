//
//  parameters.h
//  Sparse-relay-optimization-ADMM
//
//  Created by Ayoub Saab  on 2019-11-18.
//  Copyright © 2019 Ayoub Saab. All rights reserved.
//

#ifndef parameters_h
#define parameters_h

#include<cmath>

/*    system paramters  */
const unsigned int UEs=4;   //number of users
const unsigned int relays=4;//number of relays
const unsigned int antenna[relays]={2,2,2,2}; //number of antennas at each relay
const unsigned int power[relays]={4,4,4,4};//maximum power at each relay


/* channel and noise parameters */
const  float sigmaChannel = sqrt(0.5); //channel I/Q std deviation, total = 2*sigmaChannel^2 = 1;
const float sigmaR = sqrt(0.005); //relay-station noise std deviation
const float sigmaDUE = sqrt(0.0001); //destination-UE noise std deviation
const float sigmaUE = sqrt(2);

/* optimization parameters  */
const unsigned int distortionlessConstraint = 1; //pre-defined constraints
const double rho = 20;    //penalty parameter of ADMM
const unsigned int lambda[relays] = {10,10,10,10}; //regularization parameter for each relay
double eps_abs = 0;   //ADMM absolute tolerance metric
double eps_rel = pow(10,-4);    //ADMM relative tolerance metric


/*simulation paremeters */
const int monteCarlo = pow(10,0);  //number of Monte Carlo simulations
int coherenceTime = 5;    //number of Transmission-reception instances for each MonteCarlo run



#endif /* parameters_h */
