#include <iostream>
#include <math.h>
#include <fstream>
//#include <filesystem>
#include "2.h"
#include "2.cc"
using namespace std;

double exponential_function(double x) {
    return (exp(-x) /x);
}

double sin_function(double x) {
    return (x * sin(1/x));
}

double integrate(double (*f)(double), double (*integration_algorithm)(double (*func)(double), int, double, double), double a, double b){
        double temp =  integration_algorithm(f, 1, a, b); // first iteration step is used to compare second itertation step
        for (int N = 1; N<1000; N=N*2) {
            double solution = integration_algorithm(f, N, a, b);
        if (temp - solution <10^(-4)){
            cout << "Necessary N is: " << N << endl;
            return solution;
        }else {
            temp = solution;
        }
    }
}

int main(){
    cout << "Part a) of exercise 3: " << endl;
    double a = 1;
    double b = 100;
    double temp = 1.;
    for(int N = 1; N<1000; N=N*2){ //half the intervall until the change in result becomes smaller than 10^-4 but cap it at 1000 incase of overflow
        double solution = trapezoid(exponential_function, N, a, b);
        if (temp - solution <10^(-4)){
            cout << "Necessary N is: " << N << endl;
            cout << "The soltution is: " << solution << endl;
            break;
        }else {
            temp = solution;
        }
    }

    double solution = integrate(exponential_function, trapezoid, a, b);
    cout << endl << "Solution calculated by the implemented function " << solution << endl;



    //for(int N = 1; N<1000; N=N*2){ //half the intervall until the change in result becomes smaller than 10^-4 but cap it at 1000 incase of overflow
    //    double solution = riemann(exponential_function, N, a, b);
    //    if (temp - solution <10^(-4)){
    //        cout << "Necessary N is: " << N << endl;
    //        cout  << solution;
    //    }else {
    //        temp = solution;
    //    }
    //}
    //double solution = simpson(exponential_function, N, a, b);         

}