#include <iostream>
#include <math.h>
using namespace std;

double f_1(double *x) {
    return sin(*x);
}

double two_point(double *x, double h) {
        for (int i=0;i<=sizeof(x); i++){
            double x_plus_h[sizeof(x)];
            double x_minus_h[sizeof(x)];
            x_plus_h[i] = x[i]+h;
            x_minus_h[i] = x[i]-h;
        }
        double dif_f = (f_1(x_plus_h) - f_1(x_minus_h))/2*h;
    return dif_f;
}

int main() {

    return 0;
}