#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <random>
using namespace std;

template<int size>
class Ising_Metropolis{
    private:
        random_device rd;
        mt19937 generator;
        uniform_real_distribution<double> p;
        uniform_int_distribution<> dist_i, dist_j;

        void saveToCSV(const std::string& filename);
        double Energy(int, int, int);
    public:
        Ising_Metropolis() 
        : 
            generator(rd()),
            p(0,1),
            dist_i(0,size-1),
            dist_j(0,size-1)
        {};
        double k_BT = 1;
        double beta = 1/k_BT;
        int Spins[size][size];
        void set_k_BT(double);
        void initialise_rnd();
        void initialise_orderly();
        void equilibrate(int);
        void measure(int, int);
        
    

};
template <int size> void Ising_Metropolis<size>::initialise_rnd(){
    for(int i=0; i < size; i++){
        for(int j=0; j < size; j++){
            Spins[i][j] = (p(generator) < 0.5) ? -1 : 1;
        }

    }
};

template <int size> void Ising_Metropolis<size>::set_k_BT(double new_k_BT){
    k_BT = new_k_BT;
    beta = 1/k_BT;

}

template <int size> void Ising_Metropolis<size>::initialise_orderly(){
    for(int i=0; i < size; i++){
        for(int j=0; j < size; j++){
            Spins[i][j] = 1;
        }

    }
};

template <int size> void Ising_Metropolis<size>::equilibrate(int n_sweeps){
    
    for(int sweeps = 0; sweeps < n_sweeps; sweeps++){
        for(int iterator = 0; iterator < size*size; iterator++){
            int i = dist_i(generator);
            int j = dist_j(generator);
            int spin = Spins[i][j];
            int old_spin = spin;
            spin *= -1;
            double Delta_E = Energy(i,j, spin) - Energy(i,j, old_spin);
            if(Delta_E < 0 || p(generator) < exp(-beta*Delta_E)){       
                Spins[i][j] = spin;
            }
            else{
                Spins[i][j] = old_spin;
    
            }
        }
    }
};

template <int size> void Ising_Metropolis<size>::measure(int n_equilibrate, int n_sweeps){

    equilibrate(n_equilibrate);
    for(int sweeps = 0; sweeps < n_sweeps; sweeps++){
        for(int iterator = 0; iterator < size*size; iterator++){
            int i = dist_i(generator);
            int j = dist_j(generator);
            int spin = Spins[i][j];
            int old_spin = spin;
            spin *= -1;
            double Delta_E = Energy(i,j, spin) - Energy(i,j, old_spin);
            if(Delta_E < 0 || p(generator) < exp(-beta*Delta_E)){       
                Spins[i][j] = spin;
            }
            else{
                Spins[i][j] = old_spin;
    
            }
        }
        if(sweeps%(n_sweeps/5)==0){
            saveToCSV("k_BT:"+to_string(int(k_BT))+"_sweep:"+to_string(sweeps/(n_sweeps/5)));
        }

    
    }
};



template <int size> double Ising_Metropolis<size>::Energy(int i, int j, int spin_loc){
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


template<int size> void Ising_Metropolis<size>::saveToCSV(const std::string& filename) {
    std::ofstream outputFile("build/"+filename+".csv");

    if (outputFile.is_open()) {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                outputFile << Spins[i][j];
                if (j != size - 1) {
                    outputFile << ",";
                }
            }
            outputFile << "\n";
        }

        outputFile.close();
        std::cout << "CSV file saved successfully.\n";
    } else {
        std::cout << "Unable to open the file.\n";
    }
};

int main(){

int const size = 100;
Ising_Metropolis<size> Test1;
Test1.initialise_rnd();
Test1.measure(100,1e4);
Test1.set_k_BT(3.0);
Test1.measure(100,1e4);




return 0;
}