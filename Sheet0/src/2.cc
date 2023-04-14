#include <iostream>
#include <math.h>
using namespace std;

double functiona(double x) {
    if (x==1.0) {
        cout << "x should not be 1";
        return 1;
    } else
    return (1/sqrt(x) - 1/sqrt(x-1));
}

double functiona_new(double x) {
        if (x==1.0) {
        cout << "x should not be 1";
        return 1;
    } else
    return ((sqrt(x+1) -sqrt(x))/(sqrt(x+1)*sqrt(x)));
}

double functionb(double x) {
    return ((1-cos (x))/sin (x));
}

double functionb_new(double x) {
    return tan(x/2);
}


double functionc(double x, double delta) {
    return sin(x + delta) - sin(x);
}

double functionc_new(double x, double delta) {
    return (-sin(x) + cos(x)*sin(delta) + cos(delta)*sin(x));
}

double relativ_error(double x_1,double x_2) {
    return abs(x_1 - x_2)/x_2;
}

int main() {
    double output1;
    double output2;
    double error;
    double start_1;
    double start_2;
    double start_3;
    double delta_1;
    
    cout << "choose a value for x_1\n";
    cin >> start_1;
    output1 = functiona(start_1);
    cout << "function a gives the value: ";
    cout << output1 << " \n";
    output2 = functiona_new(start_1);
    cout << "The redefined function a gives the value: " << output2 << endl;
    error = relativ_error(output1, output2);
    cout << "The relativ error between the definitions is: " << error << endl;

    cout << "choose a value for x_2\n";
    cin >> start_2; // is this bad programming in c?
    output1 = functionb(start_2);
    cout << "function b gives the value: ";
    cout << output1 << " \n";
    output2 = functionb_new(start_2);
    cout << "The redefined function a gives the value: " << output2 << endl;
    error = relativ_error(output1, output2);
    cout << "The relativ error between the definitions is: " << error << endl;


    cout << "choose a value for x_3\n";
    cin >> start_3; // is this bad programming in c?
    cout << "choose a value for delta\n";
    cin >> delta_1; // is this bad programming in c?
    output1 = functionc(start_3, delta_1);
    cout << "function b gives the value: ";
    cout << output1 << " \n";
    output2 = functionc_new(start_3, delta_1);
    cout << "The redefined function a gives the value: " << output2 << endl;
    error = relativ_error(output1, output2);
    cout << "The relativ error between the definitions is: " << error << endl;
    return 0;
}