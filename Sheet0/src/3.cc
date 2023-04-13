#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int break_count = 0;
double euler_approx [100];
double sym_euler_approx [100];
double recursive_euler(double recursion_value ,double stepwidth){
    if (break_count>100){
        break_count = 0;
        return recursion_value;
    } else {
        euler_approx[break_count] = recursion_value;
        recursion_value = recursion_value* (1-stepwidth); // formula for euler
        break_count = break_count + 1;
        return recursive_euler(recursion_value,stepwidth);
    }
}

double sym_recursive_euler(double recursion_value_n, double recursion_value_n_1 ,double stepwidth){
    if (break_count>100){
        break_count = 0;
        return recursion_value_n;
    } else {
        double temp = recursion_value_n; // make temp var that later give its value to recursion_value_n_1
        sym_euler_approx [break_count] = recursion_value_n;
        recursion_value_n = -2*stepwidth*recursion_value_n + recursion_value_n_1; //formual for the symmetric euler
        recursion_value_n_1 = temp;
        break_count = break_count + 1;
        return sym_recursive_euler(recursion_value_n, recursion_value_n_1,stepwidth);
    }
}       

int main() {
    double time_stepwidth = 0.2;
    double solution = recursive_euler(1, time_stepwidth); // call the recursive euler function

    std::ofstream myfile1;
    myfile1.open("/mnt/c/Users/Max/Desktop/Computational_Physics/Sheet0/bin/eulera.csv"); // I know the path is not automated but I dont know how that works in c++
    for (int i=0;i<100; i++){
    myfile1 << euler_approx[i] << ", "; // write the euler approximation into a .csv file
    }
    myfile1.close(); 

    double sym_solution = sym_recursive_euler(exp(-time_stepwidth), 1, time_stepwidth); // call the recursive symmetrical euler function

    std::ofstream myfile2;
    myfile2.open("/mnt/c/Users/Max/Desktop/Computational_Physics/Sheet0/bin/symeulera.csv");
    for (int i=0;i<100; i++){
    myfile2 << sym_euler_approx[i] << ", "; //write the approximation into a csv file
    }
    myfile2.close(); 

    double solution_b = recursive_euler(1-time_stepwidth, time_stepwidth); // same as above just with new start values

    std::ofstream myfile3;
    myfile3.open("/mnt/c/Users/Max/Desktop/Computational_Physics/Sheet0/bin/eulerb.csv");
    for (int i=0;i<100; i++){
    myfile3 << euler_approx[i] << ", ";
    }
    myfile3.close(); 

    double sym_solution_b = sym_recursive_euler(1-time_stepwidth, 1, time_stepwidth);

    std::ofstream myfile4;
    myfile4.open("/mnt/c/Users/Max/Desktop/Computational_Physics/Sheet0/bin/symeulerb.csv");
    for (int i=0;i<100; i++){
    myfile4 << sym_euler_approx[i] << ", ";
    }
    myfile4.close(); 
    
    return 0;
};