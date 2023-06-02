#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <complex>
#include <time.h>

using namespace std;
using namespace Eigen;
typedef MatrixXcd Dynamic2D;



int main(){

    double a = 1.0;
    double b = 1.5;
    double t_start = 0;
    double t_fin = 10;

    double Delta_x = 0.1;
    double Delta_y = Delta_x;
    double Delta_t = 0.02;
    double Delta_x_2 = Delta_x * Delta_x;
    double Delta_y_2 = Delta_y * Delta_y;
    double Delta_t_2 = Delta_t * Delta_t;


    double dim_x = a/Delta_x;
    double dim_y = b/Delta_y;
    MatrixXd u_n = MatrixXd::Zero(dim_y, dim_x);

    //initialise:
    for(int i=0; i < dim_y ; i++){
        for(int j=0; j < dim_y ; j++){
            u_n(i,j) = sin(M_PI/a*j*Delta_x)*sin(2*M_PI/b*Delta_y*i);
        }
    }

    MatrixXd u_nm1 = u_n;
    MatrixXd u_np1 = u_n;
    double timesteps = t_fin/Delta_t;

    fstream T_0("T_0.csv");
    fstream T_1("T_1.csv");
    //iterate
    for(int t=t_start; t <= timesteps; t++){
        if(t == 0)
        {
            T_0 << u_n;
        }
        if(t==timesteps-4)
        {
            T_1 << u_n;
        }
        u_nm1 = u_n;
        u_n = u_np1;
        for(int i=1; i < dim_y-1 ; i++){
            for(int j=1; j < dim_y-1 ; j++){
                u_np1(i,j) = Delta_t_2 * ((u_n(i+1,j) - 2*u_n(i,j) + u_n(i-1,j))/Delta_x_2 + (u_n(i,j+1) - 2*u_n(i,j) + u_n(i,j-1))/Delta_y_2) + 2*u_n(i,j) - u_nm1(i,j);
            }
        }

    }
    T_0.close();
    T_1.close();


    return 0;
}