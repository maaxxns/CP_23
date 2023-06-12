#include <iostream>
#include <math.h>
#include <fstream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream;

struct params{
    long long int seed;
    long long int a;
    long long int c;
    long long int m;
};

long long int pseudo_rng(params parameter){
    return ((parameter.a*parameter.seed + parameter.c) % parameter.m);
};

double pseudo_rng_double(params parameter){
    return (double(pseudo_rng(parameter))/parameter.m);
};

void iteration_pseudo_rng_with_save(int iterations, params parameter, string path){
    long long int r_n_1;
    double r_n_1_float;
    path = "bin/" + path;
    ofstream myfile;
    myfile.open(path.c_str());
    myfile << "# seed " << parameter.seed << ", a " << parameter.a << ", c " << parameter.c << ", m " << parameter.m << endl;
    for (int i = 0; i < iterations; i++){
        r_n_1 = pseudo_rng(parameter);
        parameter.seed = r_n_1;
        r_n_1_float = double(r_n_1)/parameter.m;
        myfile << r_n_1_float << endl;
    }
    myfile.close();
}


int main() {
    unsigned int N = 10e5;
    {
        params a;
        a.seed = 1234;
        a.a = 20;
        a.c = 120;
        a.m = 6075;
        iteration_pseudo_rng_with_save(N, a, "rng_number_a.csv");
    }
    {
        params b;
        b.seed = 1234;
        b.a = 137;
        b.c = 187;
        b.m = 256;
        iteration_pseudo_rng_with_save(N, b, "rng_number_b.csv");
    }
    {
        params c;
        c.seed = 123456789;
        c.a = 65539;
        c.c = 0;
        c.m = 2147483648;
        iteration_pseudo_rng_with_save(N, c, "rng_number_c.csv");
    }
    {
        params d;
        d.seed = 1234;
        d.a = 16807;
        d.c = 0;
        d.m = 2147483648 - 1;
        iteration_pseudo_rng_with_save(N, d, "rng_number_d.csv");
    }
    return 0;
}