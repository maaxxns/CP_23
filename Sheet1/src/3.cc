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
    if (x == 0){
        return 0;
    }else{
    return (x * sin(1/x));
}
}

double integrate(double (*f)(double), double (*integration_algorithm)(double (*func)(double), int, double, double), double a, double b){
        double temp =  1; // first iteration step is used to compare second itertation step
        for (int N = 1; N<10000; N=N*2) {
            double solution = integration_algorithm(f, N, a, b);
        if (abs(temp - solution) < 1e-4){
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

    double solution1 = integrate(exponential_function, trapezoid, a, b);
    cout << "Solution calculated by the implemented trapezoid " << solution1 << endl;

    double solution2 = integrate(exponential_function, riemann, a, b);
    cout << "Solution calculated by the implemented riemann " << solution2 << endl;

    double solution3 = integrate(exponential_function, simpson, a, b);
    cout << "Solution calculated by the implemented simpson " << solution3 << endl;

    cout << "Part b) of exercise 3: " << endl;
    double a2 = 0;
    double b2 = 1;

    double solution4 = integrate(sin_function, trapezoid, a2, b2);
    cout << "Solution calculated by the implemented trapezoid " << solution4 << endl;

    double solution5 = integrate(sin_function, riemann, a2, b2);
    cout << "Solution calculated by the implemented riemann " << solution5 << endl;

    double solution6 = integrate(sin_function, simpson, a2, b2);
    cout << "Solution calculated by the implemented simpson " << solution6 << endl;



    return 0;
}