#ifndef ONE_H
#define ONE_H
#include <fstream>
#include <Eigen/Dense>

struct params{
    long long int seed;
    long long int a;
    long long int c;
    long long int m;
};

int pseudo_rng(params parameter);

double pseudo_rng_double(params parameter);

Eigen::VectorXd iteration_pseudo_rng(int iterations, params parameter);

#endif