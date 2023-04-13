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

double functionb(double x) {
    return ((1-cos (x))/sin (x));
}

double functionc(double x, double delta) {
    return sin(x + delta) - sin(x);
}

int main() {
    double end;
    double start;
    double delta;
    
    cout << "choose a value for x_1\n";
    cin >> start;
    end = functiona(start);
    cout << "function a gives the value: ";
    cout << end << " \n";

    cout << "choose a value for x_2\n";
    cin >> start; // is this bad programming in c?
    end = functionb(start);
    cout << "function b gives the value: ";
    cout << end << " \n";

    cout << "choose a value for x_3\n";
    cin >> start; // is this bad programming in c?
    cout << "choose a value for delta\n";
    cin >> delta; // is this bad programming in c?
    end = functionc(start, delta);
    cout << "function b gives the value: ";
    cout << end << " \n";

    return 0;
}