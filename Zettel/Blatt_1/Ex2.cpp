#include <iostream>
#include <math.h>
#include <functional>
#include <fstream>
using namespace std;


// a)
double trapezoid(function<double(double)> f,double a, double b, int N =100)
{ 

    double h = (b-a)/N;
    cout << "Intervalwidth trap:"<< h << endl;
    double sum = 0;
    for(int i=0; i<N ; i++)
    {
        sum = sum + h/2*(f(a+i*h)+f(a+i*h+h));
    }
    return sum;
}

// b)
double riemann(function<double(double)> f,double a, double b, int N =100)
{
    double h = (b-a)/N;
    cout << "Intervalwidth riemann:"<< h << endl;
    double sum = 0;
    for (int i=0; i<N; i++)
    {
        sum = sum + f(a+i*h+0.5*h)*h;
    }
    return sum;
}

// c)

double simpson(function<double(double)> f,double a, double b, int N =100)
{
    if(N%2 ==1)
    {
        N = N+1;
        cout << "N uneven. Changed N to N+1" << endl;
    }
    
    double h = (b-a)/N;
    cout << "Intervalwidth simpson:"<< h << endl;
    double sum = 0;
    for (int i=0; i<N/2; i++)
    {
        sum = sum + (1.0/3.0)*h*(f(a+(2*i+2)*h)+4*f(a+(2*i+1)*h)+f(a+2*i*h));
    }
    return sum;
}

double testFunc(double x){
    return -x*x*x;
}
double sin_func(double x){
    return sin(x);
}



int main(){



    cout << trapezoid(&testFunc, -3,2,100) << endl;
    cout << riemann(&testFunc, -3,2,100)<< endl;
    cout << simpson(&testFunc, -3,2,100)<< endl;

    cout << trapezoid(&sin_func, -3,2,100) << endl;
    cout << riemann(&sin_func, -3,2,100)<< endl;
    cout << simpson(&sin_func, -3,2,100)<< endl;


    return 0;
}