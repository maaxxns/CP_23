#include <iostream>
#include <vector>
#include <math.h> 
#include <fstream>

using namespace std;

//a)

double* Gauss_Seidel(double Phi[], double Rho[], double Delta, int J, int L, double Phi_0y, double Phi_1y, double Phi_x0, double Phi_x1, double epsilon)
{

    // set boundry conditions
    for(int j=0; j<=J; j++){
        Phi[(J+1)*j] = Phi_x0;
        Phi[(J+1)*j+L] = Phi_x1;
    }
    for(int l = 0; l<=L; l++){
        Phi[l] = Phi_0y;
        Phi[(J+1)*J+l] = Phi_1y;
    }

    //algorithm
    double Phi_pre[(J+1)*(L+1)]; //safe pre
    int count_epsilon; //variable to count the number of Phi_pre-Phi < epsilon (If count == L+J, the change of the iteration is small enought)
    int gridpoints = (J-1)*(L-1); // counting the number of grid points without the bounddry points
    int count_it = 0; //variable to count the number of iterations needed

    do{
        count_it ++;
        count_epsilon = 0;
        for(int j = 1; j<= J-1; j++){
            for(int l=1; l<=L-1; l++){
                Phi_pre[(J+1)*j+l] = Phi[(J+1)*j+l];
                Phi[(J+1)*j+l] = 0.25 * (Phi[(J+1)*(j+1)+l]+Phi[(J+1)*(j-1)+l]+Phi[(J+1)*j+(l+1)]+Phi[(J+1)*j+(l-1)])+ 0.25 * Delta * Delta * Rho[(J+1)*j+l];
                if(abs(Phi_pre[(J+1)*j+l]-Phi[(J+1)*j+l]) < epsilon){
                    count_epsilon ++;  // counting the number of gridponts with Phi_(n+1)-Phi_n < epsilon
                }
                
            }
        }
    }while(count_epsilon < gridpoints);
    cout << "Nr. of iterations:" << count_it;

    return Phi;
}

int main(void)
{
    double Delta = 0.05;
    double epsilon = pow(10,-5);
    int J = 1/Delta;
    int L = 1/Delta;
    int N_gridpoints = (J+1)*(L+1);

    //b)

    double Phi[N_gridpoints];
    for (int j = 0; j <= J; j++){
        for (int l = 0; l <= L; l++){
            Phi[J*j+l] = 1.0;
        }
    }
    double Rho[(J+1)*(L+1)];
    for (int j = 0; j <= J; j++){
        for (int l = 0; l <= L; l++){
            Phi[J*j+l] = 0;
        }
    }

    double* Phi_Solution_a = Gauss_Seidel(Phi, Rho, Delta, J, L, 0,0,0,0, epsilon);

    ofstream Phi_file("build/Phi_a");
    for(int k= 0; k<= N_gridpoints; k++){
        Phi_file << Phi_Solution_a[k] << endl;
    }
    Phi_file.close();


    return 0;
}
