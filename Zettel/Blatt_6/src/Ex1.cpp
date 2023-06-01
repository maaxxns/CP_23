#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;


class Bifurkationsdiagramm{
    public:
    Bifurkationsdiagramm(double);
    double tolerance = 0.01; // tolerance of beeing "equal"
    double x_0;
    double r_max_log = 4;
    double r_max_kub = 3;
    double Delta_r = 0.001;
    unsigned int count = 0;
    double N;
    vector<double> r_inf;




    
    int count_fix(double r);
    void save_fixpoints_log(double r_max);
    void save_fixpoints_kub(double r_max);

    private: 
        double x;

        vector<double> fix_point;
        void logistisch_warm(double r, double x_0);

        void kubisch_warm(double r, double x_0);
        
        



};

Bifurkationsdiagramm::Bifurkationsdiagramm(double x_0): x_0(x_0)
{

}

void Bifurkationsdiagramm::logistisch_warm(double r,double x_0){
    x = x_0;
    for(unsigned int i = 0; i < N; i++){
        x = r*x*(1-x);
    }
}

void Bifurkationsdiagramm::kubisch_warm(double r, double x_0){
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

void Bifurkationsdiagramm::save_fixpoints_log(double r_max){
    vector<double> x_0 = {0.1,0.3,0.5,0.7,0.9};
    ofstream r_inf_data("./build/r_inf_log.csv");
    ofstream r_double_data("./build/r_double_log.csv");

    for(int i=0; i<size(x_0);i++){
        ofstream r_data_csv("./build/r_data_log"+to_string(i)+".csv");

        count = 0;
        double count_pre;
        int n = 0;

        while(count <= 64){
            fix_point = {};
            count_pre = count;
            double r = Delta_r * n;
            count = 0;
            logistisch_warm(r,x_0[i]);
            do{
                fix_point.push_back(x);
                x = r*x*(1-x);
                r_data_csv << r << "," << x << endl;
                count++;
                n++;
                if(count > 64){
                    r_inf_data << x_0[i] << "," << r << endl;
                    break;
                }
                
            }while(!(x-tolerance < fix_point[0] && fix_point[0] < x+tolerance));
            if(count == 2*count_pre){
                cout << x_0[i] << ':' << r << ',' << count << endl;
                
            }
        }
        
        r_data_csv.close();
    }
    r_inf_data.close();
    r_double_data.close();

}



void Bifurkationsdiagramm::save_fixpoints_kub(double r_max){
    vector<double> x_0 = {0.1,0.3,0.5,0.7,0.9};
    ofstream r_inf_data("./build/r_inf_kub.csv");
    ofstream r_double_data("./build/r_double_kub.csv");

    for(int i=0; i<size(x_0);i++){
        ofstream r_data_csv("./build/r_data_kub"+to_string(i)+".csv");

        count = 0;
        int n = 0;

        while(count <= 64){
            fix_point = {};
            double r = Delta_r * n;
            count = 0;
            kubisch_warm(r,x_0[i]);
            do{
                fix_point.push_back(x);
                x = r*x - x*x*x;
                r_data_csv << r << "," << x << endl;
                count++;
                n++;
                if(count > 64){
                    r_inf_data << x_0[i] << "," << r << endl;
                    break;
                }
            }while(!(x-tolerance < fix_point[0] && fix_point[0] < x+tolerance));
        }
        r_data_csv.close();
    }
    r_inf_data.close();
    r_double_data.close();

}

int main(){

    double r = 3.5698;

    double x_0 = 0.5;


    Bifurkationsdiagramm logistisch(x_0);
    logistisch.N = 30;
    logistisch.save_fixpoints_log(4);

    Bifurkationsdiagramm kubisch(x_0);
    kubisch.N = 30;
    kubisch.save_fixpoints_kub(3);
    



    return 0;
}