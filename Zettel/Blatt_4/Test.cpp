#include <iostream>
#include <vector>
#include <math.h> 
#include <fstream>



using namespace std;

double* Gauss_Seidel(double Delta){
    double array[3] = {1,2,3};
    return array;
}

int main(void){

    int L = 2;
    int J = 2;
    int array1[3][3] = {
        {1,2,3},
        {4,5,6},
        {7,8,9}
    };
    int array2[9] = {9,8,7,6,5,4,3,2,1};
    for (int j =0 ; j<=J ; j++){
        for (int l = 0; l <=L; l++)
        {
            cout << array1[j][l] << array2[l+(J+1)*j] << endl;
            /* code */
        }
        
    }
    double* a = Gauss_Seidel(3);
    //cout << a[1] << endl;

    return 0;

}