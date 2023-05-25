#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream;

struct Start_values{
    MatrixXf Starting_condition;
    double v_0;
};

class Waveequation
{
    public:
        Waveequation (Vector2d L, double t_max, double delta_t, double Delta_x, double Delta_y, const string& filename); // constructor to initilize the 1D Strip and the alocate memory for the Matrix
    // functions
        void Integration(double t_max, double Delta_x, double Delta_y, Start_values Starting_condition); // Actual measurment function

    private:
        Vector2d L; // length of the array in x and y Strip
        double t_max; // maximum time 
        double Delta_x; // the discretization of the length in x
        double Delta_y; // the discretization of the length in y
        double delta_t; // the discretization of the time
        vector<MatrixXd> u;//u = {u_0(x, y), u_1(x, y), u_2(x,y),...} a vctor of matrices that saves the functions value at space (x,y) at time t
        MatrixXd u_t;
        void save(int steps); // save to csv function
        // Matrix indices are accesed by the (,) operator not with [][] as in python
        ofstream myfile1;
        ofstream myfile2;
        string filename;
};

Waveequation::Waveequation(Vector2d L, double t_max, double delta_t, double Delta_x, double Delta_y, const string& filename):
    L(L),
    delta_t(delta_t),
    Delta_x(Delta_x),
    Delta_y(Delta_y),
    t_max(t_max),
    filename(filename),
    u(3), // alocate space 
    u_t(int(L[0]/Delta_x), int(L[1]/Delta_y))
{
    for(int i = 0; i < u.size(); i++){
        u[i] = u_t;
    }
    cout << u_t.rows() <<" x "<< u_t.cols() << endl;
    ofstream myfile1;
    ofstream myfile2;
}

void Waveequation::save(int steps){
    string path = "bin/" + filename;
    if(!myfile1.is_open()){ // test if file is open
        myfile1.open(path.c_str()); // if not open the file
    }
    for (int i = 0; i < u[steps].rows(); i++ ){
        for (int k = 0; k < u[steps].cols(); k++){ // i am not 100% certain about the u.cols()
            myfile1 << u[steps](i, k) << ", ";
        }
    }
    myfile1 << endl;
}

void Waveequation::Integration( double t_max, double Delta_x, double Delta_y, Start_values Starting_values){
    for (int i = 0; i < Starting_values.Starting_condition.rows(); i++){
        for(int k = 0; k < Starting_values.Starting_condition.cols(); k++){
            u[1](i, k) = Starting_values.Starting_condition(i, k); // starting condition for the system so u(x,y,0)
        }
    }

    for (int i = 0; i < Starting_values.Starting_condition.rows(); i++){
        for(int k = 0; k < Starting_values.Starting_condition.cols(); k++){
            u[0](i, k) = -2 * delta_t * Starting_values.v_0 + u[1](i, k); // t = 0 is what t=-1 is in the script
        }
    }

    for (int steps = 1; steps < int(t_max/delta_t) - 1; steps++){ // time iteration
        for(int i = 0; i < int(L[0]/Delta_x); i++){
            u[1](i, 0) = 0; // fixed frame
            u[1](i, int(L[1]/Delta_y) - 1) = 0;
        }
        for(int j = 0; j < int(L[1]/Delta_y); j++){
            u[1](0, j) = 0; // fixed frame
            u[1](int(L[0]/Delta_x) - 1, j) = 0;            
        }
        for(int i = 1; i < int(L[0]/Delta_x) - 1; i++){ // the x coordinate for the membrane
            for(int j = 1; j < int(L[1]/Delta_y) - 1; j++){ // y coordinate of the membrane
                u[1 + 1](i, j) = 2*u[1](i, j) - u[1 - 1](i, j) 
                + pow(delta_t,2)
                * ((u[1](i + 1, j) - 2*u[1](i, j) + u[1](i - 1, j))/pow(Delta_x, 2) 
                + (u[1](i, j + 1) - 2*u[1](i, j) + u[1](i, j - 1))/pow(Delta_y, 2)); // the scheme from the sheet with c=1
            }
        }
        if(steps % 100 == 0){
            save(2);
        }
        u[0] = u[1];
        u[1] = u[2];
        cout << "\r" << double(steps)/((t_max/delta_t) - 1) *100 << "%";
    }
    cout << endl;
    myfile1.close();
    myfile2.close();
}


int main(void) {
    const Vector2d L {1., 1.5};
    const double t_max = 4;
    const double Delta_x = 0.01;
    const double Delta_y = 0.01; // Delta_x = Delta_y
    const double delta_t = 1e-05; // smaller than stability criteon

    {
        MatrixXf Starting_condition(int(L[0]/Delta_x), int(L[1]/Delta_y));
        for (int i = 0; i < int(L[0]/Delta_x); i++){
            for(int j = 0; j < int(L[1]/Delta_y); j++){
                Starting_condition(i,j) = sin(M_PI/L[0] * Delta_x * i) * sin(2 * M_PI/L[1] * Delta_y * j);
            }
        }
        double v_0 = 0;
        Start_values starting;
        starting.Starting_condition = Starting_condition; // u(x,y,0)
        starting.v_0 = v_0; // du(x,y,t=0)/dt 

        Waveequation waveequation(L, t_max, delta_t, Delta_x, Delta_y, "wiggle_that_membrane.csv");
        waveequation.Integration(t_max, Delta_x, Delta_y, starting);
    }
    return 0;
}
