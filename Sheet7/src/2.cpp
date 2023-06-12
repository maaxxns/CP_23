#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream;
class Householder
{
    public:
        Householder(int N);
        void apply();
        void save(string path, int i);
    private:
        VectorXd v(int i);
        VectorXd u(int i);
        MatrixXd S(int i);
        MatrixXd A;
        int N;
};


Householder::Householder(int N):
    A(N, N),
    N(N)
{

    for (int i = 0; i < N; i++){ // A without the delta
        for (int j = 0; j < N; j++){
            A(i ,j) = double(i) + double(j);
        }
    }

    for (int i = 0; i < N; i++){ // the diagonal with the delta
        A(i, i) = 3. * double(i);
    }
}

VectorXd Householder::v(int i){
    VectorXd a(N - i);
    for (int j = (i - 1); j < N - 1; j++){
        a(j - (i - 1)) = A(j + 1, i - 1);
    }
    return a;
}

VectorXd Householder::u(int i){
    VectorXd e_i(N - i);
    VectorXd v_i = v(i);
    double k;
    k = sqrt(v_i.dot(v_i)); // norm
    if(v_i[0] < 0){
        k = -k;
    } 
    for (int j = i; j < N; j++){
        e_i(j - i) = j == i; // zeros but for j == i
    }
    return ((v_i - k * e_i) / sqrt((v_i - k * e_i).dot((v_i - k * e_i)))); // norm
}

MatrixXd Householder::S(int i){
    MatrixXd Ones(N - i, N - i);
    for (int n = 0; n < N - i; n++){
        Ones(n,n) = 1.;
    }
    VectorXd u_i = u(i);
    MatrixXd S_n;
    S_n = Ones - 2.*(u_i * u_i.transpose());
    return S_n;
}

void Householder::apply(){
    for(int i = 0; i < N - 2; i++){
        MatrixXd P = MatrixXd::Identity(N, N);
        MatrixXd S_n = S(i + 1);
        for (int m = (i + 1); m < N; m++){
            for (int n = (i + 1); n < N; n++){
                P(m, n) = S_n(m - (i + 1), n - (i + 1));
            }
        }
        A = (P * A * P); // no need to transpose P as P = P^T
    }
}

void Householder::save(string path, int i){
    path = "bin/" + to_string(i) + path;
    ofstream myfile;
    myfile.open(path.c_str());
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N - 1; j++){
            if(A(i,j) < 0.00001){
                myfile << 0.000 << "\t\t\t";
            }else{
            myfile << A(i, j) << "\t\t\t";
            }
        }
        myfile << A(i, N - 1) << endl;
    }
    myfile.close();
}

int main() {
    {
    int i = 100;
    Householder householder(i);
    householder.save("_N_Matrix_before.csv", i);
    householder.apply();
    householder.save("_N_Matrix.csv", i);
    }
    return 0;
};