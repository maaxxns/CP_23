#include <iostream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <chrono>
using namespace std;
using namespace Eigen;
using std::ofstream;
using std::complex;

// constants hbar = 1.054571817 * 1e-34;


class Schroedinger
{
    public:
       Schroedinger(Vector2d L, double Delta_x, double delta_t, string path1); // constructor
       void time_evolution(double t_max); // applies the time evolution operator to the wavefunction
       void wavefunction(double sigma = 1., double xi_0 = 1); // calcultes the wave function
    
    private:
       Matrix< complex<double>, Dynamic, Dynamic > Hamilton(); // calculates the Hamilton operator
       Matrix< complex<double>, Dynamic, Dynamic > Ones(); // makes a nxm Matrix filled with 1.
       Matrix< complex<double>, Dynamic, Dynamic > A(double delta_t); // makes a the time evolution matrix
       void save(int steps); // saves the calculated entries of the time evolved wave function
       Vector2d L; // length of the array in x 
       double Delta_x; // the discretization of the length in x
       double delta_t;
       vector< VectorXcd > psi;//the wavefunction that solves the schroedinger equation
       VectorXcd psi_n; //entrances of psi
       string path;
       ofstream file;
};

Schroedinger::Schroedinger(Vector2d L, double Delta_x, double delta_t, string path1):
    L(L),
    psi(3),
    delta_t(delta_t),
    Delta_x(Delta_x),
    psi_n(int((L[1] - L[0])/Delta_x))
{
    for (int i = 0; i < psi.size(); i++){
        psi[i] = psi_n; // reserve memory for atelast three entcrances
    }
    path = "bin/" + path1;
    ofstream file;
}

void Schroedinger::time_evolution(double t_max){
    psi[1] = psi[0]; // initial condition
    save(1);
    Matrix< complex<double>, Dynamic, Dynamic > A_1 = A(delta_t);
    for(int steps = 0; steps < int(t_max/delta_t); steps++){
        psi[2] = A_1 * psi[1]; // time evolution matrix
        if(steps % 1 == 0){//save just every xth frame
            save(2); // save it
        }
        psi[1] = psi[2]; // and make it the next step
        cout << "\r" << double(steps)/(int(t_max/delta_t)) *100 << "%";
    }
    file.close();
}

void Schroedinger::save(int steps){
    if(!file.is_open()){ // test if file is open
        file.open(path.c_str()); // if not open the file
    }
    double sum;
    for (int i = 0; i < psi[steps].size(); i++){
        file << pow(abs(psi[steps][i]), 2) << ", "; // now write every entranc of the vector in the file
        sum += pow(abs(psi[steps][i]), 2);
    }
    file << endl; // after that end the line
    if(steps == 1){
        ofstream myfile1;
        myfile1.open("bin/hamilton.csv");
        Matrix< complex<double>, Dynamic, Dynamic > M = Hamilton();
        for (int i = 0; i < M.rows(); i++){
            for(int j = 0; j < M.cols(); j++){
                myfile1 << M(i,j) << ", ";
            }
            myfile1 << endl;
        }
        myfile1.close();
    }
}

Matrix< complex<double>, Dynamic, Dynamic > Schroedinger::A(double delta_t1){ // makes the time evolution matrix
    delta_t = delta_t1;
    complex<double> del_t {delta_t1, 0.};
    Matrix< complex<double>, Dynamic, Dynamic > A;
    Matrix< complex<double>, Dynamic, Dynamic > H = Hamilton();
    Matrix< complex<double>, Dynamic, Dynamic > One = Ones();
    complex<double> i{0., 1.};
    A = (One + i/2. * del_t * H).inverse() * (One - i/2. * del_t * H); // hbar = 1
    return A;
}

Matrix< complex<double>, Dynamic, Dynamic > Schroedinger::Hamilton(){
    const int N = int(abs(L[1]- L[0])/Delta_x);
    Matrix< complex<double>, Dynamic, Dynamic> Hamilton(N,N); // matrix nxm
    //Hamilton.resize(N,N);
    for (int n = 0; n < N; n++){ 
            Hamilton(n, n) = -1./pow(Delta_x, 2) * ( -2. ) + pow(Delta_x, 2) * pow(n, 2); // if the potential is supossed to be symmetrical add + int(L[0]/Delta_x) in pow(n,2)
            if(n > 0){
                Hamilton(n, n - 1) = -1./pow(Delta_x,2)*( 1 );
            }
            if(n < N - 1){
                Hamilton(n, n + 1) = -1./pow(Delta_x,2)*( 1 );
            }
        }       
    return Hamilton;
}

Matrix< complex<double>, Dynamic, Dynamic > Schroedinger::Ones(){ // makes a nxm Matrix filled with 1.
    MatrixXd Ones(int(abs(L[1]- L[0])/Delta_x), int(abs(L[1]- L[0])/Delta_x));
    for (int n = 0; n < int(abs(L[1]- L[0])/Delta_x); n++){ 
        Ones(n, n) = 1.;
    }
    return Ones;
}

void Schroedinger::wavefunction(double sigma, double xi_0){
    complex<double> sum;
    for (int i = 0; i < int((L[1] - L[0])/Delta_x); i++){
        psi[0][i] = pow(1./(2. * M_PI * sigma), 0.25) * exp( -pow((Delta_x * double(i) + L[0]) - xi_0 , 2)/(4 * sigma)); // not sure if i - L[0] is correct
    }
    for (int i = 0; i < int((L[1] - L[0])/Delta_x); i++){
        sum += pow(abs(psi[0][i]),2); // norm with the squared norm
    }
    cout << "Before norm " << sum << endl;
    psi[0] = 1./sqrt(sum) * psi[0]; // apply squared norm
    sum = 0.;
    for (int i = 0; i < int((L[1] - L[0])/Delta_x); i++){
        sum += pow(abs(psi[0][i]),2); // test if its applied correctly
    }
    cout << "After norm " << sum << endl;
}


int main() {
    Vector2d L(-10., 10.);
    double Delta_x = 0.1;
    double t_max = 10.;
    double delta_t = 0.02; // works with x = 0.1 t = 0.002
    string path = "schroedinger.csv";


    Schroedinger Schroedinger(L, Delta_x, delta_t, path);
    Schroedinger.wavefunction(); // make the initial wave function
    Schroedinger.time_evolution(t_max);

    return 0;
}