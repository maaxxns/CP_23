#include <iostream>
#include <math.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream;
using std::complex;

// constants hbar = 1.054571817 * 1e-34;

class Schroedinger
{
    public:
       Schroedinger(); // constructor
       Matrix< complex<double>, Dynamic, Dynamic > Hamilton(); // calculates the Hamilton operator
       void time_evolution(double t_max, double delta_t); // applies the time evolution operator to the wavefunction
       void save(); // saves the calculated entries of the time evolved wave function
       void wavefunction(); // calcultes the wave function
       Matrix< complex<double>, Dynamic, Dynamic > A(); // makes a the time evolution matrix
       Matrix< complex<double>, Dynamic, Dynamic > Ones(); // makes a nxm Matrix filled with 1.
    private:
        Vector2d L; // length of the array in x 
        double t_max; // maximum time 
        double Delta_x; // the discretization of the length in x
        double delta_t; // the discretization of the time
        vector< VectorXd > psi;//the wavefunction that solves the schroedinger equation
};

void Schroedinger::time_evolution(double t_max, double delta_t){
    for (int steps = 0; steps < int(t_max/delta_t); steps++){
        psi[steps + 1] = A() * psi[steps];
    }
}

Matrix< complex<double>, Dynamic, Dynamic > Schroedinger::A(){ // makes a the time evolution matrix
    Matrix< complex<double>, Dynamic, Dynamic > A;
    complex<double> i{0., 1.};
    A = (Ones() - i/2. * delta_t * Hamilton()) * (Ones() + i/2. * delta_t * Hamilton()); // hbar = 1
    return A;
}

Matrix< complex<double>, Dynamic, Dynamic > Schroedinger::Hamilton(){
    MatrixXd Hamilton(int(abs(L[1]- L[0])/Delta_x), int(abs(L[1]- L[0])/Delta_x)); // matrix nxm 
    for (int n = 0; n < int(abs(L[1]- L[0])/Delta_x); n++){ 
        for (int m = 0; m < int(abs(L[1]- L[0])/Delta_x); m++){
            if(n == m){ // delta_nm
                Hamilton(n, m) = -1./pow(Delta_x, 2) * ( -2. ) + pow(Delta_x, 2) * pow(n, 2);
                Hamilton(n, m - 1) = -1./pow(Delta_x,2)*( 1 );
                Hamilton(n, m + 1) = -1./pow(Delta_x,2)*( 1 );
            }
        }       
    }
    return Hamilton;
}

Matrix< complex<double>, Dynamic, Dynamic > Schroedinger::Ones(){ // makes a nxm Matrix filled with 1.
    MatrixXd Ones(int(abs(L[1]- L[0])/Delta_x), int(abs(L[1]- L[0])/Delta_x));
    for (int n = 0; n < int(abs(L[1]- L[0])/Delta_x); n++){ 
        for (int m = 0; m < int(abs(L[1]- L[0])/Delta_x); m++){
            Ones(n, m) = 1.;
        }
    }
    return Ones;
}
int main() {

    return 0;
}