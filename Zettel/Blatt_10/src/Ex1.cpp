#include <iostream>
#include <math.h>
#include <fstream>
#include <random>
using namespace std;


double Hamilton(int sigma, double H){
    return -sigma * H;
}



int main(){

    ofstream data("build/Ex1_data.csv");
    mt19937 rnd;
    uniform_real_distribution<double> dist(0.0, 1.0);

    int steps = 10000;
    double H = -5.0; // initial value of the extern B-field
    double Nr_H_values = 1e4; //number of H values
    double beta = 1.0; 
    double H_step = 10.0/Nr_H_values; 

    for(int H_iterator;  H_iterator < Nr_H_values ; H_iterator++){
        int magnetisation_sum = 0;
        int spin = (dist(rnd) < 0.5) ? -1 : 1; //begin with a random spin direction 
        for(int t = 0; t< steps; t++){

            int old_spin = spin; // save old spin
            spin *= -1; // switch spin direction
            double Delta_E = Hamilton(spin, H) - Hamilton(old_spin, H); // calc Delta E
            if(Delta_E < 0 || dist(rnd) < exp(-beta*Delta_E)){
                magnetisation_sum += spin; //accept new spin and add to magnetisation
            }
            else{
                spin = old_spin; // reject new spin
                magnetisation_sum += spin; // add old spin
            }
        }
        double magnetisation = static_cast<double>(magnetisation_sum)/steps; // calc average magnetisation 

        data << H << ' ' << magnetisation << ' ' << tanh(beta*H) << endl;
        H = H + H_step; // set new external field

    }
    data.close();

    return 0;
}