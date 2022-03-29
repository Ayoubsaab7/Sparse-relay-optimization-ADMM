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

void slice (const MatrixXcd A, int rowA, int rowB, int colA,int colB,MatrixXcd &A_sliced){
    int idx_i=0;
    for (int i=rowA;i<=rowB ; i++){
        int idx_j=0;
        for(int j=colA;j<=colB;j++){
            A_sliced(idx_i,idx_j) = A(i,j);
            idx_j++;
        }
        idx_i++;
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

#endif /* functions_h */
