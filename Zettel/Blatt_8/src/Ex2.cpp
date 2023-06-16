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

long long int pseudo_rng(params parameter){ // the pseudo rng number as int
    return ((parameter.a*parameter.seed + parameter.c) % parameter.m);
};
 
double pseudo_rng_double(params parameter){ // pseudo rng as double
    return (double(pseudo_rng(parameter))/parameter.m);
};

VectorXd iteration_pseudo_rng(int iterations, params parameter){ // and now make a lot rng numbers
    long long int r_n_1;
    VectorXd y(iterations);
    for (int i = 0; i < iterations; i++){
        r_n_1 = pseudo_rng(parameter);
        parameter.seed = r_n_1;
        //y[i] = double(r_n_1)/parameter.m;
        y[i] = r_n_1;
    }
    return y;
}

Vector2d box_muller(params x_1, params x_2){ // the box mueller scheme
    Vector2d y;
    double v_1 = 2. * pseudo_rng_double(x_1) - 1.;
    double v_2 = 2. * pseudo_rng_double(x_2) - 1.;
    double R = v_1 * v_1 + v_2 * v_2; // should be smaller than 1
    double R_sqrt = sqrt(R);
    if(R < 1.){
        y[0] = sqrt(-2. * log(R)) * v_1 / R_sqrt;
        y[1] = sqrt(-2. * log(R)) * v_2 / R_sqrt;
    }
    return y;
}

double p_1(double x){
    return sin(x)/2.;
}

double p_2(double x){
    return 3 * x * x;
}

double p_2_inverse(double x){
    return 1./6. * pow(x/3., -1./2.);
}

int main() {
    const int N = int(10e5);
    {//a)
        params x_1;
        x_1.seed = 1234;
        x_1.a = 16807;
        x_1.c = 0;
        x_1.m = 2147483648 - 1;

        params x_2;
        x_2 = x_1;
        x_2.seed = 69;

        VectorXd rng_numbers1;
        rng_numbers1 = iteration_pseudo_rng(N/2, x_1);
        VectorXd rng_numbers2;
        rng_numbers2 = iteration_pseudo_rng(N/2, x_2);

        ofstream myfile1;
        myfile1.open("bin/2_box_muller.csv");
        for (int i = 0; i < N/2; i++){
            Vector2d save;
            x_1.seed = rng_numbers1[i];
            x_2.seed = rng_numbers2[i];
            save = box_muller(x_1, x_2);
            myfile1 << save << endl;
        }

        myfile1.close();
    }
    { //c)
        params x_1;
        x_1.seed = 1234;
        x_1.a = 16807;
        x_1.c = 0;
        x_1.m = 2147483648 - 1;

        params x_2;
        x_2 = x_1;
        x_2.seed = 42; // set some starting values

        ofstream myfile2;
        myfile2.open("bin/2_rejection.csv");
        int i = 0;
        double x;
        double y;
        while(i < N){ // make sure that we draw enough numbers
            x_1.seed = pseudo_rng(x_1);
            x = double(x_1.seed) / double(x_1.m);
            x_2.seed = pseudo_rng(x_2);
            y = double(x_2.seed) / double(x_2.m);
            if(y < p_1(x)){
                myfile2 << x << ", " << y << endl;
                i++;
            }
        }
        cout << endl;
        myfile2.close();
    }
    { // d)
        params x_1;
        x_1.seed = 1234;
        x_1.a = 16807;
        x_1.c = 0;
        x_1.m = 2147483648 - 1;

        ofstream myfile3;
        myfile3.open("bin/2_inverse.csv");
        int i = 0;
        double y;
        double x;
        while (i < N){
            x_1.seed = pseudo_rng(x_1);
            x = double(x_1.seed) / double(x_1.m);
            y = p_2(x);
            if(y > 1e-5){
            myfile3 << (y) << ", " << p_2_inverse(y) << endl;
            i++;
            }
        }
        myfile3.close();
    }
    return 0;
};