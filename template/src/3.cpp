#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

     

int main() {
    int size;
    double var[size];

    std::ofstream myfile1;
    myfile1.open("/mnt/c/Users/Max/Desktop/Computational_Physics/Sheet0/bin/eulera.csv"); // I know the path is not automated but I dont know how that works in c++
    for (int i=0;i<100; i++){
    myfile1 << var[i] << ", "; // write the euler approximation into a .csv file
    }
    myfile1.close(); 


    return 0;
};