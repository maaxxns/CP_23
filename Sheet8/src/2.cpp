#include <iostream>
#include <math.h>
#include <fstream>
#include <Eigen/Dense>
#include "1.hpp"
using namespace std;
using namespace Eigen;
using std::ofstream;

Vector2d box_muller(params x_1, params x_2, int N){
    VectorXd y(N);
    VectorXd Ones(N);
    VectorXd R_sqrt(N/2);
    for(int i = 0; i < N/2; i++){
        Ones[i] = 1.;
    }
    VectorXd v_1 = 2. * iteration_pseudo_rng(N/2, x_1) - Ones;
    VectorXd v_2 = 2. * iteration_pseudo_rng(N/2, x_2) - Ones;
    VectorXd R = v_1 * v_1 + v_2 * v_2;
    for(int i = 0; i < N/2; i++){
        R_sqrt[i] = sqrt(R[i]);
        if(R[i] < 1.){
            y[2 * i] = sqrt(-2 * log(R[i])) * v_1[i] / R_sqrt[i];
            y[2 * i + 1] = sqrt(-2 * log(R[i])) * v_2[i] / R_sqrt[i]; // maybe this will throw segmentation error
        }else cout << "At point " << i << "R > 1";
    }
    return y;
}



int main() {


    return 0;
};