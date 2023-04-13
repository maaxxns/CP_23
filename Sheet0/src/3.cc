#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int break_count = 0;

double recursive_euler(double recursion_value ,double stepwidth){
    break_count = break_count + 1;
    if (break_count>100){
        break_count = 0;
        return recursion_value;
    } else {
        cout << recursion_value <<",";
        std::ofstream myfile;
        myfile.open ("C:/Users/Max/Desktop/Computational_Physics/Sheet0/bin/euler.csv");
        myfile << recursion_value << ", ";
        myfile.close(); // I know this is very unefficient but its the easiest way I could think of 
        recursion_value = recursion_value* (1-stepwidth); // formula for euler
        return recursive_euler(recursion_value,stepwidth);
    }
}

double sym_recursive_euler(double recursion_value_n, double recursion_value_n_1 ,double stepwidth){
    break_count = break_count + 1;
    if (break_count>100){
        break_count = 0;
        return recursion_value_n;
    } else {
        double temp = recursion_value_n; // make temp var that later give its value to recursion_value_n_1
        cout << recursion_value_n <<",";
        recursion_value_n = -2*stepwidth*recursion_value_n + recursion_value_n_1; //formual for the symmetric euler
        recursion_value_n_1 = temp;
        return sym_recursive_euler(recursion_value_n, recursion_value_n_1,stepwidth);
    }
}       

int main() {
    double time_stepwidth = 0.2;
    cout <<"Here starts the Euler" << endl;
    double solution = recursive_euler(1, time_stepwidth); // I later noticed that the function could have been a void, but now I dont want to change it
    cout << endl << endl <<"Here starts the symmetrical Euler" << endl;
    double sym_solution = sym_recursive_euler(exp(-time_stepwidth), 1, time_stepwidth);
    cout << endl << endl <<"Here starts the b) part of 3" << endl;
    cout <<"Here starts the Euler" << endl;
    double solution_b = recursive_euler(1-time_stepwidth, time_stepwidth);
    cout << endl << endl <<"Here starts the symmetrical Euler" << endl;
    double sym_solution_b = sym_recursive_euler(1-time_stepwidth, 1, time_stepwidth);
    return 0;
};