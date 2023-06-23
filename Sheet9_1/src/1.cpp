#include <iostream>
#include <math.h>
#include <fstream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream;

struct params_intervall{
    double x_0;
    double y_0;
    double z_0;
};

struct params_newton{
    double x_n;
    int i = 0;
};

void save_intervall(params_intervall r, ofstream& myfile){
    if(myfile.is_open()==false){
        myfile.open("bin/intervall.csv");
    }
    myfile << r.x_0 << ", "  << r.y_0 << ", " << r.z_0 << endl;
}

void save_newton(params_newton r, ofstream& myfile){
    if(!myfile.is_open()){
        myfile.open("bin/newton.csv");
    }
    myfile << r.x_n << ", "  << r.i << endl;
}
void half_intervall(double (*f)(double), params_intervall r, double precission){
    double u_n;
    ofstream myfile1;
    while(r.z_0 - r.x_0 > precission){
        save_intervall(r, myfile1);
        if(r.y_0 - r.x_0 > r.z_0 - r.y_0){
            u_n = (r.x_0 + r.y_0) / 2.;
            if(f(u_n) < f(r.y_0)){
                r.z_0 = r.y_0;
                r.y_0 = u_n;
            // r.z_0 = r.z_0;
            }else{
                r.x_0 = u_n;
            }
        }else{
            u_n = (r.y_0 + r.z_0) / 2.;
            if(f(u_n) < f(r.y_0)){
                r.x_0 = r.y_0;
                r.y_0 = u_n;
            // r.z_0 = r.z_0;
            }else{
                r.z_0 = u_n;
            }
        }
    }
    save_intervall(r, myfile1);
    myfile1.close();
}

void newton(double (*f)(double), params_newton r, double precission){
    double h = 1./1000000.;
    ofstream myfile2;
    while(abs(f(r.x_n)) > precission){
        save_newton(r, myfile2);
        double f_ = (f(r.x_n - 2. * h) - 8. * f(r.x_n - h) + 8. * f(r.x_n + h) - f(r.x_n + 2. * h)) / (12. * h);
        r.x_n = r.x_n - f(r.x_n)/(f_);
        r.i = r.i + 1;
        cout << f(r.x_n) << endl;
    }
    save_newton(r, myfile2);
    myfile2.close();
}

double f(double x){
    return x * x - 2;
}

int main() {
    params_intervall r_intervall;
    r_intervall.x_0 = -0.5;
    r_intervall.y_0 = -0.1;
    r_intervall.z_0 = 2.;
    double precission = 1.e-9;
    half_intervall(f, r_intervall, precission);

    params_newton r_newton;
    r_newton.x_n = 1.;
    newton(f, r_newton, precission);

    return 0;
}