#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <time.h> 
using namespace std;


Eigen::MatrixXd Build_Matrix(int N){
    Eigen::MatrixXd A(N,N);
    for(int k = 0 ; k<N; k++){
        for(int l = 0 ; l<N; l++){
            A(k,l) = k + l ; 
        }
        A(k,k) = A(k,k) + k ;
    }
    //A(0,0) = 1;
    return A;

}

// a)
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
    P(0,0) = 1.0; //top left corner is Id
    
    
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




//b)
Eigen::MatrixXd Jacobi_Rotation(Eigen::MatrixXd A, int n){
    int N = A.cols();
    Eigen::MatrixXd P = Eigen::MatrixXd::Identity(N,N);
    double t;
    double c;
    double s;
    if(A(n,n) == 0){ //preventing the case of dividing by 0
        c = 0;
        s = 1;
    }
    else{

    t = A(n,n+1)/A(n,n);
    c = 1/(sqrt(t*t+1));
    s = t*c;
    }

    int q,p;
    p = n;
    q = n+1;
    P(p,p) = c;
    P(p,q) = s; 
    P(q,p) = -s;  
    P(q,q) = c; 
    return P;
}
Eigen::MatrixXd QR_Iteration(Eigen::MatrixXd Tri){
    int N = Tri.cols(); // dimension
    for(int i=0; i<N; i++){ //using N iterations
        Eigen::MatrixXd Q_t = Jacobi_Rotation(Tri, 0); // first step
        
        for(int n = 1; n<N-1; n++){
            Q_t = Jacobi_Rotation(Tri, n) * Q_t; // multilpying all Jacobi-Roation matrices (depndend on Tri) to calc Q_t 
        }
        Tri = Q_t * Tri * Q_t.transpose(); // calc new Tri
    }
    
    return Tri;
}   



// c)
Eigen::VectorXd Ev_func(Eigen::VectorXd v, Eigen::MatrixXd A , int N){
    Eigen::VectorXd w;
    int I = v.size();
    Eigen::VectorXd Ev(I);
    
    for(int i; i<I ;i++ ){
        for(int n = 0; n<N ; n++){
            w = A * v;
            v = w/w.norm();
        }
        Ev[i] = v.dot(A * v); 
        Eigen::MatrixXd V = v * v.transpose();
        A = A - Ev[i]*V;
    }
    return Ev;
}


int main(){
    double time1=0.0, tstart;
    
    int N = 10; // dimension N

    //prepare the csv 
    ofstream A_start_10("./build/A_start_10.csv");
    ofstream Tri_10("./build/Tri_10.csv");
    ofstream Diag_10("./build/Diag_10.csv");
    ofstream Ev_House_10("./build/Ev_House_10.csv");
    ofstream Ev_Pot_10("./build/Ev_Pot_10.csv");


    // build matrix A:
    Eigen::MatrixXd A;
    A = Build_Matrix(N);
    A_start_10 << A << endl;
    A_start_10.close();


    tstart = clock();

//a)
    Eigen::MatrixXd Tri = Householder_Iteration(A);
    Tri_10 << Tri << endl;
    Tri_10.close();

//b)
    Eigen::MatrixXd Diag = QR_Iteration(Tri);
    Diag_10 << Diag << endl;
    Diag_10.close();
    
    Eigen::VectorXd Ev_House = Diag.diagonal();
    Ev_House_10 << Ev_House << endl;
    Ev_House_10.close();


    time1 += clock() - tstart;
    time1 = time1/CLOCKS_PER_SEC;
    cout << "  time House = " << time1 << " sec." << endl;


//c)
    tstart = clock();

    Eigen::VectorXd v_0 = Eigen::VectorXd::Ones(N);
    Eigen::VectorXd Ev_pot = Ev_func(v_0, A , 15);
    Ev_Pot_10 << Ev_pot << endl;
    Ev_Pot_10.close();

    time1 += clock() - tstart;
    time1 = time1/CLOCKS_PER_SEC;
    cout << "  time Pot = " << time1 << " sec." << endl;


    return 0;
}