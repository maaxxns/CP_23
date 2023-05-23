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
        Waveequation (Vector2d L, double t_max, double delta_t, double Delta_x, double Delta_y); // constructor to initilize the 1D Strip and the alocate memory for the Matrix
    // functions
        void Integration(double t_max, double Delta_x, double Delta_y, Start_values Starting_condition); // Actual measurment function
        void save(const string& filename); // save to csv function

    private:
        Vector2d L; // length of the array in x and y Strip
        double t_max; // maximum time 
        double Delta_x; // the discretization of the length in x
        double Delta_y; // the discretization of the length in y
        double delta_t; // the discretization of the time
        vector<MatrixXd> u;//u = {u_0(x, y), u_1(x, y), u_2(x,y),...} a vctor of matrices that saves the functions value at space (x,y) at time t
        MatrixXd u_t;
        // Matrix indices are accesed by the (,) operator not with [][] as in python
};

Waveequation::Waveequation(Vector2d L, double t_max, double delta_t, double Delta_x, double Delta_y):
    L(L),
    delta_t(delta_t),
    Delta_x(Delta_x),
    Delta_y(Delta_y),
    t_max(t_max),
    u(int(t_max/delta_t)), // alocate space 
    u_t(int(L[0]/Delta_x), int(L[1]/Delta_y))
{
    for(int i = 0; i < u.size(); i++){
        u[i] = u_t;
    }
    cout << u_t.rows() <<" x "<< u_t.cols() << endl;
}

void Waveequation::save(const string& filename){
    ofstream myfile1;
    ofstream myfile2;
    string path = "bin/" + filename;
    myfile1.open(path.c_str());
    for (int t = 0; t < u.size(); t++){ // time
        if(t == 0){
            myfile1 << "# rows ar time, the rest of the data is dummed in one big line" << endl;
        }
        for (int i = 0; i < u[t].rows(); i++ ){
            for (int k = 0; k < u[t].cols(); k++){ // i am not 100% certain about the u.cols()
                myfile1 << u[t](i, k) << ", ";
            }
        }
        myfile1 << endl;
    }
    myfile2.close();
    string path2 = "bin/compressed" + filename;
    myfile2.open(path2.c_str());
    int iterator = int(u.size()/100);
    for (int t = 0; t < 100; t++){ // time
        if(t == 0){
            myfile2 << "# rows ar time, the rest of the data is dummed in one big line" << endl;
        }
        for (int i = 0; i < u[t].rows(); i++ ){
            for (int k = 0; k < u[t].cols(); k++){ // i am not 100% certain about the u.cols()
                myfile2 << u[t*iterator](i, k) << ", ";
            }
        }
        myfile2 << endl;
    }
    myfile2.close();

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
            u[steps](i, 0) = 0; // fixed frame
            u[steps](i, int(L[1]/Delta_y) - 1) = 0;
        }
        for(int j = 0; j < int(L[1]/Delta_y); j++){
            u[steps](0, j) = 0; // fixed frame
            u[steps](int(L[0]/Delta_x) - 1, j) = 0;            
        }
        for(int i = 1; i < int(L[0]/Delta_x) - 1; i++){ // the x coordinate for the membrane
            for(int j = 1; j < int(L[1]/Delta_y) - 1; j++){ // y coordinate of the membrane
                u[steps + 1](i, j) = 2*u[steps](i, j) - u[steps - 1](i, j) 
                + pow(delta_t,2)
                * ((u[steps](i + 1, j) - 2*u[steps](i, j) + u[steps](i - 1, j))/pow(Delta_x, 2) 
                + (u[steps](i, j + 1) - 2*u[steps](i, j) + u[steps](i, j - 1))/pow(Delta_y, 2)); // the scheme from the sheet with c=1
            }
        }
    }
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

        Waveequation waveequation(L, t_max, delta_t, Delta_x, Delta_y);

        waveequation.Integration(t_max, Delta_x, Delta_y, starting);
        waveequation.save("wiggle_that_membrane.csv");
    }
    return 0;
}