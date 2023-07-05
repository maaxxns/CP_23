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
    double Energy(int, int, int);
    

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
    double Delta_E = Energy(i,j, spin) - Energy(i,j, old_spin);
    if(Delta_E < 0 || dist(rnd) < exp(-beta*Delta_E)){
        magnetisation_sum += spin;
    }
    else{
        spin = old_spin;
        magnetisation_sum += spin;
    }


};
template <unsigned int size> double Ising_Metropolis<size>::Energy(int i, int j, int spin_loc){
    int energy_sum;
    int spin_top; 
    int spin_down; 
    int spin_left; 
    int spin_right;

    if(i==0){
        spin_top = Spins[size-1][j];
        if(j==0){
            spin_left = Spins[i][size-1];
        }else{
            spin_left = Spins[i][j-1];
        }
        if(j==size-1){
            spin_right = Spins[i][0];
        }else{
            spin_right = Spins[i][j+1];
        }
    }else{
        spin_top = Spins[i-1][j];
        if(j==0){
            spin_left = Spins[i][size-1];
        }else{
            spin_left = Spins[i][j-1];
        }
        if(j==size-1){
            spin_right = Spins[i][0];
        }else{
            spin_right = Spins[i][j+1];
        }
    }
    if(i==size-1){
        spin_down = Spins[0][j];
    }else{
        spin_down = Spins[i+1][j];
    }
    energy_sum = - spin_loc * spin_down - spin_loc * spin_left - spin_loc * spin_top - spin_loc * spin_right;
    return energy_sum;
}


int main(){


Ising_Metropolis<10> Test1;
Test1.initialise_rnd();
for(int i = 0; i<10; i++){
    for(int j = 0; j<10; j++){
            cout << i << ','<< j << ',' << Test1.Energy(i,j) << endl;
    }

}



return 0;
}