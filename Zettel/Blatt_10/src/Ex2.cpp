#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <random>
using namespace std;

template<unsigned int size>
class Ising_Metropolis{
    double k_BT = 1;
    double beta = 1/k_BT;
    mt19937 rnd;
    public:
    Ising_Metropolis();

    int Spins[size][size];

    void initialise_rnd();
    void initialise_orderly();
    void equlibrate();
    double Energy(int, int);
    

};
template <unsigned int size> Ising_Metropolis<size>::Ising_Metropolis(){

};
template <unsigned int size> void Ising_Metropolis<size>::initialise_rnd(){
    uniform_real_distribution<double> dist(0, 1);
    for(int i=0; i < size; i++){
        for(int j=0; j < size; j++){
            Spins[i][j] = (dist(rnd) < 0.5) ? -1 : 1;
        }

    }
};
template <unsigned int size> void Ising_Metropolis<size>::initialise_orderly(){
    for(int i=0; i < size; i++){
        for(int j=0; j < size; j++){
            Spins[i][j] = 1;
        }

    }
};

template <unsigned int size> void Ising_Metropolis<size>::equlibrate(){
    
    uniform_real_distribution<int> dist_i(0, size-1);
    uniform_real_distribution<int> dist_j(0, size-1);
    int i = dist_i(rnd);
    int j = dist_j(rnd);
    int spin = Spins[i][j];
    int old_spin = spin;
    spin *= -1;
    double Delta_E = Hamilton(spin, H) - Hamilton(old_spin, H);


};
template <unsigned int size> double Ising_Metropolis<size>::Energy(int i, int j){
    int energy_sum;
    int spin_loc = Spins[i][j];
    int spin_top; 
    int spin_down; 
    int spin_left; 
    int spin_right;

    if(i==0 && j==0){
        spin_top = Spins[size-1][j];
        spin_left = Spins[i][size-1];
    if(i==0 && j==size-1){
        spin_top = Spins[size-1][j];
        spin_right = Spins[i][0];
    }
    if(i==size-1 && j==0){
        spin_down = Spins[size-1][j];
        spin_left = Spins[i][0];
    }

        energy_sum = }

    return energy_sum;
}


int main(){


Ising_Metropolis<10> Test1;
Test1.initialise_rnd();


return 0;
}