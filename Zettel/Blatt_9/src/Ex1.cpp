#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
using namespace std;


double func(double x){
    return (x*x)-2;
};

double func_deriv(double x){
    return 2*x;
};

//Halbierung:
void Intervallhalbierung(double (*func)(double), double x_e, double x, double y, double z){
    ofstream points("./build/halbierung.csv");
    double u;
    int steps = 0;
    points << steps << ' ' << x << ' ' << y << ' ' << z << endl;
    while((z-x)>x_e){
        if(y-x<z-y){
            u = (y+z)/2;
            if(func(u)<func(y)){
                x = y;
                y = u;
                z = z;
            }
            else{
                x = x;
                y = y;
                z = u;
            }
        }
        else{
            u = (x+y)/2;
            if(func(u)<func(y)){
                x = x;
                z = y; 
                y = u;
            }
            else{
                x = u;
                y = y;
                z = z;
            }
        }
         
        steps++;
        points << steps << ' ' << x << ' ' << y << ' ' << z << endl;
    }
    points.close();

};

// numeric calc of derivatives(4 point method)
double firstDerivative(double (*function)(double), double x, double h) {
    double f_x_plus_2h = function(x + 2. * h);
    double f_x_plus_h = function(x + h);
    double f_x_minus_h = function(x - h);
    double f_x_minus_2h = function(x - 2. * h);

    return (-f_x_plus_2h + 8. * f_x_plus_h - 8. * f_x_minus_h + f_x_minus_2h) / (12. * h);
}


double secondDerivative(double (*function)(double), double x, double h) {
    double f_x_plus_h = function(x + h);
    double f_x = function(x);
    double f_x_minus_h = function(x - h);

    return ((f_x_plus_h - 2. * f_x + f_x_minus_h) / (h * h));
}

//Newton:

void Newton(double (*func)(double), double x_e, double x){
    ofstream points("./build/newton.csv");
    double df ,ddf;
    double x_old;
    int steps = 0;
    do{
        x_old = x;
        points << steps << ' ' << x << endl;
        df = firstDerivative(func, x,  0.1);
        ddf = secondDerivative(func, x,  0.1);
        x = x_old - (df/ddf);
        steps++;
    }while(abs(x_old-x) > x_e); //Konvergenzbedingung 
    points.close();
};

int main(){
    double x_0 = -0.5;
    double y_0 = -0.1;
    double z_0 = 2.;
    double x_e = 1e-9;

    double x_0N = 1.;

    Intervallhalbierung(&func,x_e,x_0,y_0,z_0);
    Newton(&func, x_e, x_0N);


    return 0;
}