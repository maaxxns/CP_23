#include <iostream>
#include <math.h>
#include <fstream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream;



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

Eigen::VectorXd Box_Mueller(unsigned long long int r_0, unsigned long long int a, unsigned long long int c, unsigned long long int m, unsigned long long int N){
    Eigen::VectorXd v_1 = Rand_Gen(r_0, a ,c ,m ,N/2);
    Eigen::VectorXd v_2 = Rand_Gen(r_0+10, a ,c ,m ,N/2);

    Eigen::VectorXd r_vec(N);
    v_1 = 2*v_1 - Eigen::VectorXd::Ones(N/2);
    v_2 = 2*v_2 - Eigen::VectorXd::Ones(N/2);
    for(int i = 0; i<N/2; i++){
        double v1 = v_1(i);
        double v2 = v_2(i);
        double R_2 =  v2*v2+v1*v1;
        if(R_2 < 1.0){
            double R = sqrt(R_2);
            double y_1 = sqrt(-2*log(R_2))*v1/R;
            double y_2 = sqrt(-2*log(R_2))*v2/R;
            r_vec(i) = y_1;
            r_vec(i+N/2) = y_2;

        }
        else{continue;}

    }
    return r_vec;

}



int main() {
    const int N = 1000;
    unsigned long long int r_0 = 1234;
    unsigned long long int a = 16807;
    unsigned long long int c = 0;
    unsigned long long int m = 2147483648 - 1;


    ofstream Box_data("./build/Box_data.csv");
    //a)
    r_0 = 1234;
    a = 20;
    c = 120;
    m = 6075;
    Eigen::VectorXd r_n = Box_Mueller(r_0,a,c,m,N);
    Box_data << r_n << endl;
    Box_data.close();
    return 0;
};