#include <iostream>
#include <math.h>
#include <fstream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream; 

class Ising{
    public:
        Ising(int num_particels, bool random, double T);
    private:
        double ising_move(int i, int j); 
        double calc_H(int i, int j); // calculates the Hamilon for one micro state
        double calc_H_new(int i, int j);
        void calc_Z(); // calc partition function
        double boltzmann(int i, int j);
        int N; // number of rows
        double beta = 1.; // beta = 1/(k_b*T)
        double Z; // Partition function
        double T; // Tempreature
        double V_ij;
        MatrixXi System; // The NxN dimensional system of spins. Spin up is 1, spin down is -1
};

Ising::Ising(int num_particels, bool random, double T):
    N(sqrt(num_particels)),
    System(sqrt(num_particels), sqrt(num_particels)),
    T(T),
    V_ij(1./double(2*num_particels - 1))
{
    if(random){ // random initilization if desired
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                System(i, j) = (double(rand())/RAND_MAX < 0.5) ? 0:1;
            }
        }
    }else{ // else just fill with ones
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                System(i, j) = 1;
            }
        }
    }
}

double Ising::calc_H(int i, int j){
    double H = 0.;
    H = H - System(i,j) * System(i-1, j); 
    H = H - System(i,j) * System(i+1, j); 
    H = H - System(i,j) * System(i, j-1); 
    H = H - System(i,j) * System(i, j+1); 
    return H;
}

double Ising::calc_H_new(int i, int j){
    double H = 0;
    if(System(i, j) == 0){
        H = H - 1 * System(i-1, j); 
        H = H - 1 * System(i+1, j); 
        H = H - 1 * System(i, j-1); 
        H = H - 1 * System(i, j+1); 
        return H;
    }else{
        H = H - (-1) * System(i-1, j); 
        H = H - (-1) * System(i+1, j); 
        H = H - (-1) * System(i, j-1); 
        H = H - (-1) * System(i, j+1); 
    }
}

double Ising::ising_move(int i, int j){
    double Delta_E = calc_H_new(i, j) - calc_H(i,j);
    if(Delta_E < 0){
        System(i,j) = System(i,j) == 1 ? -1 : 1; // change state if energy is better
    }else{
        double p = double(rand())/RAND_MAX;
        if(p < exp(- beta* Delta_E)){
            System(i,j) = System(i,j) == 1 ? -1 : 1;
        }
    }
}

double Ising::boltzmann(int i, int j){
    return 1./Z * exp(-beta * calc_H(i,j));
}

void Ising::calc_Z()
{
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            Z = Z + exp(-beta*calc_H(i, j));
        }
    }
}

int main() {

    return 0;
};