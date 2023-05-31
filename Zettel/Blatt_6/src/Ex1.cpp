#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;


class Bifurkationsdiagramm{
    public:
    Bifurkationsdiagramm(double);
    double tolerance = 0.0001; // tolerance of beeing "equal"
    double x_0;
    double r_max_log = 4;
    double r_max_kub = 3;
    double Delta_r = 0.001;
    unsigned int count = 0;
    double N;




    void logistisch_warm(double r);
    int count_fix(double r);
    void save_fixpoints(double r_max);

    private: 
        double x;

        vector<double> fix_point;

        void kubisch_warm(double r);
        
        



};

Bifurkationsdiagramm::Bifurkationsdiagramm(double x_0): x_0(x_0)
{

}

void Bifurkationsdiagramm::logistisch_warm(double r){
    x = x_0;
    for(unsigned int i = 0; i < N; i++){
        x = r*x*(1-x);
    }
}

void Bifurkationsdiagramm::kubisch_warm(double r){
    x = x_0;
    for(unsigned int i = 0; i < N; i++){
        x = r*x - x*x*x;
    }
}

int Bifurkationsdiagramm::count_fix(double r){
    fix_point = {};
    
    do{
    fix_point.push_back(x);
    x = r*x*(1-x);
    count++;
    }while(!(x-tolerance < fix_point[0] && fix_point[0] < x+tolerance));

    return count;
}

void Bifurkationsdiagramm::save_fixpoints(double r_max){
    
    ofstream r_data_csv ("build/r_data_csv.csv");

    count = 0;
    int n = 0;

    while(count <= 64){
        fix_point = {};
        double r = Delta_r * n;
        count = 0;
        logistisch_warm(r);
        do{
        fix_point.push_back(x);
        x = r*x*(1-x);
        r_data_csv << r << ',' << x << endl;
        count++;
        n++;
        }while(!(x-tolerance < fix_point[0] && fix_point[0] < x+tolerance));
        r_data_csv.close();

    }




}

    



//void Bifurkationsdiagramm::count_fix(){
//    unsigned int i = 0; 
//
//    do{
//        if(!(x-tolerance < fix_point[i][0] < x+tolerance)){
//            fix_point.push_back({x,fix_point[i][1]+1});
//            i++;
//        }
//        x = r*x*(1-x);
//    }while(fix_point[i-1][1] <= 10);
//}

int main(){

    double r = 3.5698;

    double x_0 = 0.5;


    Bifurkationsdiagramm logistisch(x_0);
    logistisch.N = 30;
    logistisch.save_fixpoints(4);
    



    return 0;
}