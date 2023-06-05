#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
using namespace std;


Eigen::MatrixXd Build_Matrix(int N){
    Eigen::MatrixXd A(N,N);
    for(int k = 0 ; k<N; k++){
        for(int l = 0 ; l<N; l++){
            A(k,l) = k + l;
            
        }
        A(k,k) = A(k,k) + k;
    }
    return A;

}

Eigen::MatrixXd Householder(Eigen::MatrixXd A){
    double k;
    int N = A.cols(); // N: Dimension
    Eigen::VectorXd v = A.bottomLeftCorner(N-1,1); // def of v-vector
    if(v[0] > 0){
        k = -v.norm(); // choose the same sign of -a_21
    }
    else{
        k = v.norm();
    }
    Eigen::VectorXd e = Eigen::VectorXd::Identity(N-1,1); // id-vector with dim of v
    Eigen::VectorXd u = (v-k*e)/(v-k*e).norm();
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(N-1, N-1); // identity-matrix with dim of v X v
    Eigen::MatrixXd S = Id - 2*u*u.transpose(); //calc S
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(N,N);
    P.bottomRightCorner(N-1, N-1) = S; //bootom right corner is S
    P(0,0) = 1; //top left corner is Id
    
    Eigen::MatrixXd Tri = P.transpose() * A * P; //first step of tridiagonal-matrix

    return Tri;
}

Eigen::MatrixXd Householder_Iteration(Eigen::MatrixXd A){
    int N = A.cols();
    Eigen::MatrixXd Tri = Householder(A); //first itertion

    for(int n = 1 ;n < N-1; n++){
        A = Tri.bottomRightCorner(N-n,N-n); // calc of A'
        Tri.bottomRightCorner(N-n,N-n) = Householder(A); //sets bottom right corner of Tri with Housholder(A')
    }
    return Tri;
}



int main(){
    

    int N = 5;

    Eigen::MatrixXd A;
    Eigen::MatrixXd A_;
    A = Build_Matrix(N);
    A_ = Householder(A);
    cout << A_ << endl;
    Eigen::MatrixXd Tri = Householder_Iteration(A);
    cout << Tri << endl;


    return 0;
}