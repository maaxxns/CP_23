#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream;

class Lorentz
{
    public:
        Lorentz(double r, double sigma, double b, double delta_t, double t_max, string path, Vector3d start); // Constructor
        void Runge_kutta();

    private:
        void save();
        Vector3d k_1();
        Vector3d k_2(Vector3d k_1);
        Vector3d k_3(Vector3d k_2);
        Vector3d k_4(Vector3d k_3);

        double r_1;
        double sigma;
        double b;
        double delta_t;
        double t_max;
        vector<Vector3d> R; // timewise vector of three elements that represent X, Y, Z
        ofstream file;
        string path;
};

Lorentz::Lorentz(double r, double sigma, double b, double delta_t, double t_max, string path, Vector3d start):
    r_1(r),
    sigma(sigma),
    b(b),
    delta_t(delta_t),
    t_max(t_max),
    R(2),
    path(path)
{
    R[0] = start;
    ofstream file;
};

void Lorentz::Runge_kutta()
{
    double t_n_1_2;
    for (int steps = 0; steps < int(t_max/delta_t); steps++){
        Vector3d k_1_i = k_1();
        Vector3d k_2_i = k_2(k_1_i);
        Vector3d k_3_i = k_3(k_2_i);
        Vector3d k_4_i = k_4(k_3_i);

        R[1] = R[0] + 1./6. * (k_1_i + 2* k_2_i + 2*k_3_i + k_4_i);
        save();
        R[0] = R[1];
    }
    save(); // save last position
}

void Lorentz::save(){
    if(!file.is_open()){ // test if file is open
        file.open(path.c_str()); // if not open the file
    }
    file << R[0][0] << ", " << R[0][1] << ", " << R[0][2]; // now write every entrance of the vector in the file
    // we save as X, Y, Z endl;
    file << endl; // after that end the line
}

Vector3d Lorentz::k_1(){ // do I really need to do this component wise? I dont think so but at the moment it seems easier
    Vector3d k_1;
    k_1[0] = -sigma* R[0][0] + sigma * R[0][1]; // -sigma * X +  sigma * Y
    k_1[1] = - R[0][0]*R[0][2] + r_1 * R[0][0] - R[0][1]; // - X * Z + r*X - Y
    k_1[2] = R[0][0] * R[0][1] - b * R[0][2]; //  X * Y - b * Z 
    return delta_t * k_1;
}

Vector3d Lorentz::k_2(Vector3d k_1){
    Vector3d k_2;
    k_2[0] = -sigma* (R[0][1] + delta_t/2. - R[0][0] + delta_t/2. * k_1[0]); // -sigma * X +  sigma * Y
    k_2[1] = (-( R[0][1] + delta_t/2. * k_1[1] ) + (R[0][0] + delta_t / 2.) * (r_1 - (R[0][2] + delta_t / 2.))); // - X * Z + r*X - Y
    k_2[2] = ((R[0][0]+ delta_t / 2.) * (R[0][1] + delta_t / 2.) - b * (R[0][2] + delta_t / 2.* k_1[2])); //  X * Y - b * Z 
    return delta_t * k_2;
}

Vector3d Lorentz::k_3(Vector3d k_2){
    Vector3d k_3;
    k_3[0] = sigma*(R[0][1] + delta_t / 2. - R[0][0] + delta_t / 2. * k_2[0]);
    k_3[1] = (-(R[0][1] + delta_t / 2. * k_2[1]) + (R[0][0]+ delta_t / 2.)*(r_1 - (R[0][2] + delta_t / 2.)));
    k_3[2] = ((R[0][0]+ delta_t / 2.) * (R[0][1] + delta_t / 2.) - b * (R[0][2] + delta_t / 2. * k_2[2]));
    return delta_t * k_3;
}

Vector3d Lorentz::k_4(Vector3d k_3){
    Vector3d k_4;
    k_4[0] = sigma*(R[0][1] + delta_t - R[0][0] + delta_t * k_3[0]);
    k_4[1] = (-(R[0][1]+ delta_t * k_3[1]) + (R[0][0]+ delta_t) * (r_1 - (R[0][2] + delta_t)));
    k_4[2] = ((R[0][0] + delta_t ) * (R[0][1] + delta_t) - b * (R[0][2] + delta_t * k_3[2]));
    return delta_t * k_4;
}

int main() {

    {
        double sigma = 10.;
        double b = 8./3.;
        double r = 20.;
        double t_max = 100;
        double delta_t = 0.01;
        Vector3d start;
        start << 0, 0, 0;
        string path = "bin/Lorentz_r_is_20.csv";
        Lorentz Lorentz(r, sigma, b, delta_t, t_max, path, start);
        Lorentz.Runge_kutta();
    }
    {
        double sigma = 10.;
        double b = 8./3.;
        double r = 28.;
        double t_max = 100;
        double delta_t = 0.01;
        Vector3d start;
        start << 0, 0, 0;
        string path = "bin/Lorentz_r_is_28.csv";
        Lorentz Lorentz(r, sigma, b, delta_t, t_max, path, start);
        Lorentz.Runge_kutta();
    }

    return 0;
}