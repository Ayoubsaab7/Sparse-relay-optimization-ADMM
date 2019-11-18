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
const unsigned int UEs=2;   //number of users
const unsigned int relays=2;//number of relays
const unsigned int antenna[relays]={2,2}; //number of antennas at each relay
const unsigned int power[relays]={2,2};//maximum power at each relay

/* channel and noise parameters */
const  float sigmaChannel = sqrt(0.5); //channel I/Q std deviation, total = 2*sigmaChannel^2 = 1;
const float sigmaR = sqrt(0.1); //relay-station noise std deviation
const float sigmaDUE = sqrt(0.01); //destination-UE noise std deviation


/*simulation paremeters */
const int monteCarlo = pow(10,0);  //number of Monte Carlo simulations
const double rho = 200.0;         //penalty parameter of ADMM
const unsigned int lambda[relays] = {100,100}; //regularization parameter for each relay
double eps_abs = 0;   //ADMM absolute tolerance metric
double eps_rel = pow(10,-4);    //ADMM relative tolerance metric


#endif /* parameters_h */
