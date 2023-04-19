#include <iostream>
#include <math.h>
#include <fstream>
#include <filesystem>
#include "2.h"
using namespace std;

double trapezoid(double (*f)(double), int N, double a, double b) {
    double h = (b-a)/N; // calculate h
    double x = a;
    double integral_value = integral_value + h/2 * (f(a) + f(b)); 
    for(int i=0; i<=N-2 ; i++) { // we already did the first and last 'iteration' steps so N-2 has to be used
        x = x+h;
        integral_value = integral_value + h*f(x);
    }
    return integral_value;
}

double riemann(double (*f)(double), int N, double a, double b) {
    double h = (b-a)/N; // calculate h
    double x = a;
    double integral_value;
    for(int i=0; i<=N-1; i++) {
        integral_value = integral_value + f(x)*h; // h should be the width of one intervall so there is no need to calculate x_i - x_i-1
        x = x+h;
    }
    return integral_value;
}

double simpson(double (*f)(double), int N, double a, double b) {
    if(N%2==0){ // test if N is a even number
    double h = (b-a)/N; // calculate h
    double x = a;
    double integral_value =  f(a) * 1/3; // the first function value is weigthed with 1/3
    for(int i=1; i<=N-2; i++) {
        x = x+h; // go to the next x value
        if(i%2==0){
                integral_value = integral_value + f(x) * 2/3; // the even function values are weigthed by 2/3
            }
        else{
                integral_value = integral_value + f(x) * 4/3; // the uneven function values are weigthed by 4/3
        }
    }
    integral_value = integral_value + f(b) * 1/3; // last function value is also weigthed by 1/3
    return h*integral_value;
    }else{
        cout << "For Simpson rule N should be an even number" << endl;
        return -1;
    }
}

