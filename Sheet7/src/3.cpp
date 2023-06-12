#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
int main() {
    int N = 1000;
    std::string path = "bin/Matrix_N_" + std::to_string(N) + ".csv";

    Eigen::MatrixXd A(N, N);

    for (int i = 0; i < N; i++){ // A without the delta
        for (int j = 0; j < N; j++){
            A(i ,j) = double(i) + double(j);
        }
    }

    for (int i = 0; i < N; i++){ // the diagonal with the delta
        A(i, i) = 3. * double(i);
    }

    // Tridiagonalize the matrix A
    Eigen::Tridiagonalization< Eigen::MatrixXd > triOfA(A);
    Eigen::MatrixXd T = triOfA.matrixT();
    // Print the tridiagonal matrix T
    std::ofstream myfile;
    myfile.open(path.c_str());
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N - 1; j++){
            myfile << T(i, j) << "\t\t\t";
            }
            myfile << T(i, N - 1) << std::endl;
        }
    myfile.close();
    return 0;
}
