#include <iostream>
#include <math.h>
#include <fstream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream; 

class Ising{
    public:
        Ising(int num_particels, bool random, double T, bool warmup);
        void mc_algo(int steps, string path);
        void save_state(ofstream& myfile);
    private:
        void save_energy_per_spin(ofstream& myfile);
        void save_magnetiziation(ofstream& myfile);
        void save_abs_magnetiziation(ofstream& myfile);
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

Ising::Ising(int num_particels, bool random, double T, bool warmup):
    N(sqrt(num_particels)),
    System(int(sqrt(num_particels)), int(sqrt(num_particels))),
    beta(1/T),
    V_ij(1./double(2*num_particels - 1))
{
    if(random){ // random initilization if desired
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                System(i, j) = (double(rand())/RAND_MAX < 0.5)? -1 : 1;
            }
        }
    }else{ // else just fill with ones
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                System(i, j) = 1;
            }
        }
    }
    if(warmup){
        for(int step = 0; step < int(10e4); step++){
            for (int i = 0; i < N;  i++){
                for (int j = 0; j < N; j++){
                    double p = double(rand())/RAND_MAX;
                    if(p < V_ij){
                        ising_move(i, j);
                    }
                }
            }
            cout << "Warmup: " << step/(10e4) * 100 << "%"<< "\r";
        }
        cout << endl;
    }
}

void Ising::mc_algo(int steps, string path){
    string path1 = "bin/" + path;
    ofstream myfile;
    myfile.open(path1.c_str());
    string path2 = "bin/energy_per_spin" + path;
    ofstream myfile2;
    myfile2.open(path2.c_str());
    string path3 = "bin/abs_magnet" + path;
    ofstream myfile3;
    myfile3.open(path3.c_str());
    string path4 = "bin/magnet" + path;
    ofstream myfile4;
    myfile4.open(path4.c_str());
    for(int step = 0; step < steps; step++){
        for (int i = 0; i < N;  i++){
            for (int j = 0; j < N; j++){
                double p = double(rand())/RAND_MAX;
                if(p < V_ij){
                    ising_move(i, j);
                }
            }
        }
        cout << "Processing: " << step/double(steps) * 100 << "%"<< "\r";
        save_energy_per_spin(myfile2);
        save_abs_magnetiziation(myfile3);
        save_magnetiziation(myfile4);
        save_state(myfile);
    }
    cout << endl;
}

double Ising::calc_H(int i, int j){
    double H = 0.;
    H = (i == 0)? H - System(i,j) * System(N - 1, j) : H - System(i,j) * System(i-1, j); 
    H = (i == N - 1)? H - System(i,j) * System(0, j) : H - System(i,j) * System(i+1, j); 
    H = (j == 0)? H - System(i,j) * System(i, N - 1) : H - System(i,j) * System(i, j-1); 
    H = (j == N - 1)? H - System(i,j) * System(i, 0) : H - System(i,j) * System(i, j+1); 
    return H;
}

double Ising::calc_H_new(int i, int j){
    int H = 0;
    if(System(i, j) == -1){
        H = (i == 0)? H - 1 * System(N - 1, j) : H - 1 * System(i-1, j); 
        H = (i == N - 1)? H - 1 * System(0, j) : H - 1 * System(i+1, j); 
        H = (j == 0)? H - 1 * System(i, N - 1) : H - 1 * System(i, j-1); 
        H = (j == N - 1)? H - 1 * System(i, 0) : H - 1 * System(i, j+1); 
    }else{
        H = (i == 0)? H + 1 * System(N - 1, j) : H + 1 * System(i-1, j); 
        H = (i == N - 1)? H + 1 * System(0, j) : H + 1 * System(i+1, j); 
        H = (j == 0)? H + 1 * System(i, N - 1) : H + 1 * System(i, j-1); 
        H = (j == N - 1)? H + 1 * System(i, 0) : H + 1 * System(i, j+1);  
    }
    return H;
}

double Ising::ising_move(int i, int j){
    double Delta_E = calc_H_new(i, j) - calc_H(i,j);
    if(Delta_E < 0){
        System(i,j) = (System(i,j) == 1) ? -1 : 1; // change state if energy is better
    }else{
        double p = double(rand())/RAND_MAX;
        if(p < exp(- beta* Delta_E)){
            System(i,j) = (System(i,j) == 1) ? -1 : 1;
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

void Ising::save_state(ofstream& myfile){
    if(myfile.is_open() == false){
        myfile.open("bin/dummy.csv");
    }
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            myfile << System(i,j) << ", ";
        }
    }
    myfile << endl;
}

void Ising::save_energy_per_spin(ofstream& myfile){
    if(myfile.is_open() == false){
        myfile.open("bin/dummy.csv");
    }
    double H;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            H =+ calc_H(i, j);
        }
    }
    myfile << H/double(N) << ", ";
}


void Ising::save_magnetiziation(ofstream& myfile){
    if(myfile.is_open() == false){
        myfile.open("bin/dummy.csv");
    }
    double m;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            m =+ System(i, j);
        }
    }
    myfile << m/double(N) << ", ";
}

void Ising::save_abs_magnetiziation(ofstream& myfile){
    if(myfile.is_open() == false){
        myfile.open("bin/dummy.csv");
    }
    double m;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            m =+ abs(System(i, j));
        }
    }
    myfile << m/double(N) << ", ";
}
int main() {
    {
        Ising Ising(100*100, false, 300., true);
        Ising.mc_algo(100, "state.csv");
    }
    {
        Ising Ising(100*100, true, 1., true);
        ofstream myfile("bin/state_at_T_1.csv");
        Ising.save_state(myfile);
    }
    {
        Ising Ising(100*100, true, 3., true);
        ofstream myfile("bin/state_at_T_3.csv");
        Ising.save_state(myfile);
    }

    {
        Ising Ising(100*100, true, 1.5, true);
        Ising.mc_algo(10000, "T_1_5.csv");
    }
    {
        Ising Ising(100*100, true, 3, true);
        Ising.mc_algo(10000, "T_3.csv");
    }    
    return 0;

};