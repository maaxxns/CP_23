#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
using namespace std;


double Hamilton(int& sigma, double& H){
    return sigma * H;
}

class MC{
    private:
    MC(double (*Hamilton)(int, double));
    void equilibrate();
    public:
    
};


int main(){

    double V = 0.5;
    int steps = 10;
    double H = 2.0;

    for(int t = 0; t< steps; t++){
        int sigma = 1;
        double Ham = Hamilton(sigma, H);
    }


    return 0;
}