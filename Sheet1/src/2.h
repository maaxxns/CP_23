#ifndef TWO_H
#define TWO_H


double trapezoid(double (*f)(double), int N, double a, double b);

double riemann(double (*f)(double), int N, double a, double b);

double simpson(double (*f)(double), int N, double a, double b);


#endif