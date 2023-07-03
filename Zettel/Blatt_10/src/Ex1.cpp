#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <random>
using namespace std;


double Hamilton(int sigma, double H){
    return -sigma * H;
}



int main(){

    ofstream data("build/Ex1_data.csv");
    mt19937 rnd;
    uniform_real_distribution<double> dist(0.0, 1.0);

    double V = 0.5;
    int steps = 10000;
    double H = -5.0;
    double Nr_H_values = 1e4;
    double beta = 1.0;
    double H_step = 10.0/Nr_H_values;

    for(int H_iterator;  H_iterator < Nr_H_values ; H_iterator++){
        int magnetisation_sum = 0;
        int spin = (dist(rnd) < 0.5) ? -1 : 1;
        for(int t = 0; t< steps; t++){


            int old_spin = spin;
            spin *= -1;
            double Delta_E = Hamilton(spin, H) - Hamilton(old_spin, H);
            if(Delta_E < 0 || dist(rnd) < exp(-beta*Delta_E)){
                magnetisation_sum += spin;
            }
            else{
                spin = old_spin;
                magnetisation_sum += spin;
            }
        }
        double magnetisation = static_cast<double>(magnetisation_sum)/steps;

        data << H << ' ' << magnetisation << ' ' << tanh(beta*H) << endl;
        H = H + H_step;

    }
    data.close();

    return 0;
}