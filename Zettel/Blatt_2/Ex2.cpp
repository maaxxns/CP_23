#include <iostream>
#include <math.h>
#include <functional>
#include <fstream>
using namespace std;


double f_test(double x, double y){
    return (y*y-x*x)/(y*y+x*x);
}


double Runge_Kutta(function<double(double, double)> f, double y0, double x0, double xn, double N){
    double k1, k2, k3, k4;
    double h = (xn-x0)/N;
    double t = h;
    for(int i=0; i<N;i++){
    k1 = h * f(x0+i*t,y0);
    k2 = h * f(0.5*(x0+i*t+x0+i*t+t),y0 + 0.5 * k1);
    k3 = h * f(0.5*(x0+i*t+x0+i*t+t),y0 + 0.5 * k2);
    k4 = h * f(x0+i*t+t,y0 + k3);
    y0 = y0 + 1./6.*(k1+2*k2+2*k3+k4);
    }
    return y0;
}



int main(){
    double y0 = 1.;
    double x0 = 0.;
    double xn = 0.6;
    double N = 100.;
    double h = (xn-x0)/N;
    double t = h;

    cout << Runge_Kutta(&f_test, y0, x0, xn, N);


    return 0;
}