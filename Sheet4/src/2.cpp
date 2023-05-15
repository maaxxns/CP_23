#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream;

class Diffusion
{
    public:
        Diffusion (Vector2d L, double t_max, double delta_t, double Delta_x);
    // functions
        void FTCS(double t_max, double Delta_x, vector<double> Starting_condition);
        void save(const string& filename);

    private:
        Vector2d L;
        double t_max;
        double Delta_x;
        double delta_t;
        Matrix<double, Dynamic, Dynamic> u; //u(x, t) a Matrix in which collumns the position x for every time t is saved
        // so row = x collumns = t 
        // Matrix indices are accesed by the (,) operator not with [][] as in python

};

Diffusion::Diffusion(Vector2d L, double t_max, double delta_t, double Delta_x):
    L(L),
    delta_t(delta_t),
    Delta_x(Delta_x),
    t_max(t_max),
    u(int((L[1] - L[0]) / Delta_x), int(t_max/delta_t)) // alocate space 
{
    cout << "Collumns " << u.cols() << ", rows " << u.rows() << endl;
    cout << "u(" << u.cols() << ", " << u.rows() << ")" << endl;
    //u.resize(int(t_max/delta_t), int((L[1] - L[0]) / Delta_x));
}

void Diffusion::save(const string& filename){
    ofstream myfile1;
    string path = "bin/" + filename;
    myfile1.open(path);
    for (int j = 0; j < u.rows(); j++){
        if(j == 0){
            myfile1 << "# rows are time, collumns are space" << endl;
        }
        for (int i = 0; i < u.cols(); i++ ){
            myfile1 << u(j, i) << ", ";
        }
        myfile1 << endl;
    }
    myfile1.close();
}

void Diffusion::FTCS( double t_max, double Delta_x, vector<double> Starting_condition){
    for (int i = 0; i < int((L[1] - L[0])/Delta_x); i++){
        u(i, 0) = Starting_condition[i]; // starting condition for the system so u(x,0)
    }
    for (int steps = 0; steps < int(t_max/delta_t) - 1; steps++){ // time iteration
        // I think this can also be solved by matrix multiplikation which would make it way faster
        for (int i = 0; i < int((L[1] - L[0])/Delta_x); i++){ // i have to watch out for the edges
        
        if(i > 0 && i < (L[1] - L[0])/Delta_x - 1){
            u(i, steps + 1) = u(i, steps) + delta_t/pow(Delta_x,2) * (u(i + 1, steps) - 2* u(i, steps) + u(i - 1, steps)); // FCTS time step
        }else{
            if(i == 0){
                u(i, steps + 1) = u(i, steps); // isolating boundary
            }
            if(i == (L[1] - L[0])/Delta_x - 1){
                u(i, steps + 1) = u(i, steps); // isoltation boundary
            }
        }
        }
    }
}


int main(void) {
    const Vector2d L = {0.,1.};
    const double t_max = 0.1;
    const double t_max2 = 0.005;
    const double Delta_x = 0.01;
    const double delta_t_good = 1e-06; // smaller than stability criteon
    const double delta_t_bad = 6e-05; // bigger than stability criteon
    const double delta_t_almostbad = 4e-05;

    { // Part a with u(x,0) = 1
        Diffusion Diffusion( L, t_max, delta_t_good, Delta_x);
        vector<double> Starting_condition(int((L[1] - L[0]) / Delta_x));
        for(int i = 0; i < int((L[1] - L[0]) / Delta_x); i++){
            Starting_condition[i] = 1.;
        }

        Diffusion.FTCS(t_max, Delta_x, Starting_condition);
        Diffusion.save("Diffusion_a.csv");
    }


    vector<double> Starting_delta(int((L[1] - L[0]) / Delta_x));
    for(int i = 0; i < int((L[1] - L[0]) / Delta_x); i++){
        Starting_delta[i] = 0;
    }


    { // part b with u(x,0) = delta(x - 0.5) above condition
        Starting_delta[(int((L[1] - L[0]) / (2.*Delta_x)))] = 1.; // in the middle should be a 1
        Diffusion Diffusionb(L, t_max2, delta_t_almostbad, Delta_x);
        Diffusionb.FTCS(t_max2, Delta_x,Starting_delta);
        Diffusionb.save("Diffusion_b_good.csv");
    }


    { // part b with u(x,0) = delta(x - 0.5) below condition
        Starting_delta[(int((L[1] - L[0]) / (2.*Delta_x)))] = 1.; // in the middle should be a 1
        Diffusion Diffusionb(L, t_max2, delta_t_bad, Delta_x);
        Diffusionb.FTCS(t_max2, Delta_x,Starting_delta);
        Diffusionb.save("Diffusion_b_bad.csv");
    }


    { // part c with u(x,0)=Haevyside(x - 0.5) above condition
        vector<double> starting_Haevyside(int((L[1] - L[0]) / Delta_x));
        for(int i = 0; i < int((L[1] - L[0]) / (2.*Delta_x)); i++){
            starting_Haevyside[i] = 0;
        }
        for(int i = int((L[1] - L[0]) / (2.*Delta_x)); i < int((L[1] - L[0]) / (Delta_x)); i++){
            starting_Haevyside[i] = 1;
        }
        Diffusion Diffusionc(L, t_max, delta_t_good, Delta_x);
        Diffusionc.FTCS(t_max, Delta_x, starting_Haevyside);
        Diffusionc.save("Diffusion_c_haevy.csv");
    }


    {
        vector<double> weird_function(int((L[1] - L[0]) / Delta_x));
        for(int i = 0; i < int((L[1] - L[0]) / Delta_x); i++){
            for(int n = 1; n <= 9; n++){
                if(i == int(0.1 * n / Delta_x)){
                    weird_function[i] =+ 1./9.;
                }
            }
        }

        Diffusion Diffusionc(L, t_max2, delta_t_good, Delta_x);
        Diffusionc.FTCS(t_max2, Delta_x, weird_function);
        Diffusionc.save("Diffusion_c_weird.csv");
    }
    return 0;
}