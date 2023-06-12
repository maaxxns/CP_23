#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
using namespace std;



Eigen::VectorXd Rand_Gen(unsigned long long int r_0, unsigned long long int a, unsigned long long int c, unsigned long long int m, unsigned long long int N){
    Eigen::VectorXd r_vec(N);
    if(N>=m){
        cout << "N must be smaller than m" << endl;
        return r_vec;
    }
    unsigned long long int r = r_0;
    for(unsigned int i = 0; i<N; i++){
        r = (a*r+c)% m ;
        double r_float = (double)r/(double)m;
        r_vec(i) = r_float;
    }
    
    return r_vec;
}


int main(){
    unsigned long long int r_0;
    unsigned long long int a;
    unsigned long long int c;
    unsigned long long int m;
    unsigned long long int N;

    Eigen::VectorXd r_n;

    
    N = 6000;

    ofstream a_data("./build/a_data.csv");
    //a)
    r_0 = 1234;
    a = 20;
    c = 120;
    m = 6075;
    r_n = Rand_Gen(r_0,a,c,m,N);
    a_data << r_n << endl;
    a_data.close();

    N = 250;
    ofstream b_data("./build/b_data.csv");
    //b)
    r_0 = 1234;
    a = 137;
    c = 120;
    m = 256;
    r_n = Rand_Gen(r_0,a,c,m,N);
    b_data << r_n << endl;
    b_data.close();

    
    N = 100000;

    ofstream c_data("./build/c_data.csv");
    //c)
    r_0 = 123456789;
    a = 65539;
    c = 0;
    m = 2147483648;
    r_n = Rand_Gen(r_0,a,c,m,N);
    c_data << r_n << endl;
    c_data.close();

    ofstream d_data("./build/d_data.csv");
    //d)
    r_0 = 1234;
    a = 16807;
    c = 0;
    m = 2147483647;
    r_n = Rand_Gen(r_0,a,c,m,N);
    d_data << r_n << endl;
    d_data.close();



    return 0;
}