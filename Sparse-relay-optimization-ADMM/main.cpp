//
//  main.cpp
//  Sparse-relay-optimization-ADMM
//
//  Created by Ayoub Saab on 2019-07-08.
//  Copyright © 2019 Ayoub Saab. All rights reserved.
//


//preprocessor directives
#include <iostream>
#include <iomanip>
#include <ctime>
#include <fstream>
using namespace std;

#include "functions.h"
#include "mimo.h"

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


