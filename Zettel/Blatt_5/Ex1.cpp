#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <complex>
#include <time.h>

using namespace std;
using namespace Eigen;
typedef MatrixXcd Dynamic2D;

class schrodinger{
    public:
        schrodinger(double Delta_xi, double xi_left_bound, double xi_right_bound, double Delta_tau);
        double Delta_xi;
        double xi_left_bound;
        double xi_right_bound;
        double Delta_tau;
        int dimension = (xi_right_bound-xi_left_bound)/Delta_xi;
        double sigma;
        double xi_0;

        Dynamic2D H;
        Dynamic2D Id_r = Dynamic2D::Zero(dimension, dimension);
        Dynamic2D Id_l = Id_r;
        Dynamic2D Id = Dynamic2D::Identity(dimension, dimension); //identity matrix


        VectorXcd Psi;
        Dynamic2D S_H;



        void time_evolution(double fin_time);
        void Wavefunction();
        double Norm();
        void SH();
        void Hamilton();

    private:


};

schrodinger::schrodinger(double Delta_xi, double xi_left_bound, double xi_right_bound, double Delta_tau)
{
    Delta_xi = Delta_xi;
    xi_left_bound = xi_left_bound;
    xi_right_bound = xi_right_bound;
    Delta_tau = Delta_tau;
}

void schrodinger::Hamilton(){
    
    H = Dynamic2D::Zero(dimension, dimension); // generating matrix with zeros and dimension
    double Delta_xi_2 = Delta_xi*Delta_xi; // short for Delta_xi^2
    Id_r.topRightCorner(dimension-1, dimension-1) = Dynamic2D::Identity(dimension-1, dimension-1); // right id matrix
    Id_l.bottomLeftCorner(dimension-1, dimension-1) = Dynamic2D::Identity(dimension-1, dimension-1); // left id matrix
    H = -1./Delta_xi_2*(Id_r+Id_l);

    for (int n=0; n<dimension; n++){
        H(n,n) = 2./Delta_xi_2 + n*n;
    }

}
void schrodinger::time_evolution(double fin_time){

    int timesteps = (int)fin_time/Delta_tau;
    for(int i = 0; i<= timesteps; i++){
        Psi = S_H * Psi;
    }
    
}

void schrodinger::Wavefunction(){
    double sum = 0;
    for(int j = (int)xi_left_bound/Delta_xi; j <= (int)xi_right_bound/Delta_xi; j++){
        sum += pow(exp(-(j*Delta_xi-xi_0)*(j*Delta_xi-xi_0)/(4.0*sigma)),2) * Delta_xi;
    }
    sum = sqrt(sum);
    for(int j = 0; j <= dimension; j++){
        int shift = xi_right_bound/Delta_xi;
        int k = j - shift;
        Psi(j) = 1/sum * exp(-(k*Delta_xi-xi_0)*(k*Delta_xi-xi_0)/(4*sigma));
    } 
}

double schrodinger::Norm(){
    double Norm=0;
    for(int j=0; j<= dimension; j++){
        Norm += pow(abs(Psi(j)),2) * Delta_xi; 
    }
    return Norm;
}

void schrodinger::SH(){
    S_H = (Id + 1.0i/2. * H * Delta_tau).inverse() * (Id-1.0i/2. * H * Delta_tau); 
}


double Norm_func(VectorXcd Psi, double dimension, double Delta_xi){
    double Norm = 0;
    for(int j=0; j< dimension; j++){
        Norm += pow(abs(Psi(j)),2) * Delta_xi; 
    }
    return Norm;
}

int main(){

    double Delta_xi = 0.1;
    double Delta_tau = 0.02;
    double xi_left_bound = -10;
    double xi_right_bound = 10;
    double xi_0 = 1;
    double sigma = 1;
    double fin_time = 10;

    //schrodinger Objekt(Delta_xi, xi_left_bound, xi_right_bound, Delta_tau);
    //Objekt.sigma = sigma;
    //Objekt.xi_0 = xi_0;

    //Objekt.Hamilton();
    //Objekt.Wavefunction();
    //Objekt.SH();
    //Objekt.time_evolution(fin_time);



    int dimension = (xi_right_bound-xi_left_bound)/Delta_xi;
    Dynamic2D Id_r = Dynamic2D::Zero(dimension, dimension);
    Dynamic2D Id_l = Id_r;
    Dynamic2D Id = Dynamic2D::Identity(dimension, dimension); //identity matrix



    Dynamic2D H = Dynamic2D::Zero(dimension, dimension); // generating matrix with zeros and dimension
    double Delta_xi_2 = Delta_xi*Delta_xi; // short for Delta_xi^2
    Id_r.topRightCorner(dimension-1, dimension-1) = Dynamic2D::Identity(dimension-1, dimension-1); // right id matrix
    Id_l.bottomLeftCorner(dimension-1, dimension-1) = Dynamic2D::Identity(dimension-1, dimension-1); // left id matrix
    H = -1./Delta_xi_2*(Id_r+Id_l);

    for (int n=0; n<dimension; n++){
        H(n,n) = 2./Delta_xi_2 + n*n;
    }



    double sum = 0;
    for(int j = (int)xi_left_bound/Delta_xi; j <= (int)xi_right_bound/Delta_xi; j++){
        sum += pow(exp(-(j*Delta_xi-xi_0)*(j*Delta_xi-xi_0)/(4.0*sigma)),2) * Delta_xi;
    }
    sum = sqrt(sum);

    VectorXcd Psi(dimension);
    for(int j = 0; j < dimension; j++){
        int shift = xi_right_bound/Delta_xi;
        int k = j - shift;
        Psi(j) = 1/sum * exp(-(k*Delta_xi-xi_0)*(k*Delta_xi-xi_0)/(4*sigma));
    } 

    ofstream Psi_csv("Psi.csv");
    Dynamic2D S_H = (Id + 1.0i/2. * H * Delta_tau).inverse() * (Id-1.0i/2. * H * Delta_tau); 
    int timesteps = fin_time/Delta_tau;
    for(int i = 0; i<= timesteps; i++){
        Psi = S_H * Psi;
        VectorXd Psi_2 = Psi.cwiseAbs().cwiseProduct(Psi.cwiseAbs());
        Psi_csv << Psi_2.transpose() << endl;
         //   for(int j = 0; j < dimension; j++){
       //         Psi_csv << Psi_2(j) << ',' << endl;
     //       }
    }
    Psi_csv.close();


    return 0;
}

