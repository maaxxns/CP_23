#include <iostream>
#include <math.h>
#include <functional>
#include <fstream>
using namespace std;

float two_point(float x, float h, function<float(float)> f){
    return (f(x+h)-f(x-h))/(2*h);
}
float sin_func(float x){
    return sin(x);
}
float cos_func(float x){
    return cos(x);
}

float relativ_error(float num,float ana) {
    return abs(num - ana)/ana;
}
//b)
float two_point2(float x, float h, function<float(float)> f){
    return (two_point(x+h,h,f)-two_point(x-h,h,f))/(2*h);
}
//c)
float four_point(float x, float h, function<float(float)> f){
    return (-f(x+2*h)+8*f(x+h)-8*f(x-h)+f(x-2*h))/(12*h);

}


int main(){
    float x = M_PI/4; //random x position (but not pi/2, there are numerical problems)

    ofstream Ex_1_a("bin/Ex_1_a.csv");
    for(float h=1e-5;h <= 10;h = h*2){
        Ex_1_a << h << "," << two_point(x,h, &sin_func) << endl; //safe the h's and numeric values in Ex_1_a
    }
    Ex_1_a.close();

    float h = 0.1;

    float N = 100; // number of steps in the interval [-pi,pi] 

    ofstream Ex_1_a_2("bin/Ex_1_a_2.csv");
    for(float x=-M_PI;x <= M_PI;x+=2*M_PI/N){
        Ex_1_a_2 << x << "," << relativ_error(cos(x), two_point(x,h,&sin_func)) << endl;
    }
    Ex_1_a_2.close();


//b)

    ofstream Ex_1_b("bin/Ex_1_b.csv");
    for(float h=1e-5;h <= 10;h = h*2){
        Ex_1_b << h << "," << two_point(x,h, &sin_func) << endl; //safe the h's and numeric values in Ex_1_a
    }
    Ex_1_b.close();



        //ca. same h as before


    ofstream Ex_1_b_2("bin/Ex_1_b_2.csv");
    for(float x=-M_PI;x <= M_PI;x+=2*M_PI/N){
        Ex_1_b_2 << x << "," << relativ_error(-sin(x), two_point2(x,h,&sin_func)) << endl;
    }
    Ex_1_b_2.close();

//c)



    ofstream Ex_1_c("bin/Ex_1_c.csv");
    for(float x=-M_PI;x <= M_PI;x+=2*M_PI/N){
        Ex_1_c << x << "," << relativ_error(cos(x), four_point(x,h,&sin_func)) << endl;
    }
    Ex_1_c.close();


    return 0;
}