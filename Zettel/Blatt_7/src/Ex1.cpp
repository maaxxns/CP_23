#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
using namespace std;


Eigen::Vector4d Ev_func(Eigen::Vector4d v, Eigen::Matrix4d A , int N){
    Eigen::Vector4d w;
    Eigen::Vector4d Ev;
    int I = v.size();
    for(int i; i<I ;i++ ){
        for(int n = 0; n<N ; n++){
            w = A * v;
            v = w/w.norm();
        }
        Ev[i] = v.dot(A * v); 
        Eigen::Matrix4d V = v * v.transpose();
        A = A - Ev[i]*V;
    }
    return Ev;
}


int main(){
    Eigen::Matrix4d A;
    A << 1 , -2 , -3 , 4,
        -2 , 2 , -1 , 7,
        -3 , -1 , 3 , 6,
        4 , 7 , 6 , 4;
    cout << A;
    
    // a)
    Eigen::Vector4cd Ev_a = A.eigenvalues();


    //b)
    Eigen::Vector4d v;
    Eigen::Vector4d Ev_b;
    int N = 15;
    v << 2 , 2 , 2 , 2; 
    Ev_b = Ev_func(v , A , N);


    cout << Ev_a;
    cout << Ev_b;
    

    return 0;
}