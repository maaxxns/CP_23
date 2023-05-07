#include <iostream>
#include <math.h>
#include <functional>
#include <fstream>
#include <vector>
using namespace std;



double *F_func(double r1[], double r2[], double m1, double m2){
    static double F_array[2] = {0,0};
    for(int i=0; i<=1; i++){
        if(r1[i]==r2[i]){
            F_array[i] = 0;
        }
        else{
        F_array[i] = -m1*m2*(r1[i]-r2[i])/pow(pow(r1[0]-r2[0],2)+pow(r1[1]-r2[1],2),3./2.);
        }
    }

    return F_array;
}

void euler(double r1[], double r2[], double m1, double m2, double v1[], double v2[], double h, int T_grenze){
    double *a; 
    ofstream Ex_3_euler_1("bin/Ex_3_euler_1.csv");
    ofstream Ex_3_euler_2("bin/Ex_3_euler_2.csv");
    Ex_3_euler_1 <<  r1[0] << ',' << r1[1] << ',' << v1[0] << ',' << v1[1] << endl;
    Ex_3_euler_2 <<  r2[0] << ',' << r2[1] << ',' << v2[0] << ',' << v2[1] << endl;

    for(int t;t<=T_grenze/h;t++){
        a = F_func(r1,r2,m1,m2);
        r1[0] = r1[0] + v1[0]*h;
        r1[1] = r1[1] + v1[1]*h;
        r2[0] = r2[0] + v2[0]*h;
        r2[1] = r2[1] + v2[1]*h;
        v1[0] = v1[0] + a[0]/m1*h;
        v1[1] = v1[1] + a[1]/m1*h;
        v2[0] = v2[0] + a[0]/m2*h;
        v2[1] = v2[1] + a[1]/m2*h;
        Ex_3_euler_1 <<  r1[0] << ',' << r1[1] << ',' << v1[0] << ',' << v1[1] << endl;
        Ex_3_euler_2 <<  r2[0] << ',' << r2[1] << ',' << v2[0] << ',' << v2[1] << endl;
    }
    Ex_3_euler_1.close();
    Ex_3_euler_2.close();
}


int main(){

    double r1[] = {0,1};
    double r2[] = {0,-0.5};
    double v1[] = {0.8,0};
    double v2[] = {-0.4,0};
    double m1 = 1;
    double m2 = 2;
    int T = 100; 
    double h = 0.1;

    euler(r1,r2,m1,m2,v1,v2,h,T);


    return 0;
}