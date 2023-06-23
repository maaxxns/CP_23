#include <iostream>
#include <math.h>
#include <fstream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream; 

struct params{
    double omega;
    double c_1;
    double c_2;
};

double f(double x, double y){
    return ((x - 1.9) * (x - 1.9) + (y - 2.1) * (y - 2.1) + 2. * cos(4. * x + 2.8) + 3. * sin(2. * y + 0.6));
}

class Swarm{

    public:
        Swarm(int N, Vector2d r_max, Vector2d v_max); // constructor for N particles in the box r_i^2
        void find_minimum(int Iterations, params constants, string path);
    private:
        int N;
        MatrixXd r; // particle positions for N particles
        MatrixXd v; // particle velocity for N particles
        Vector2d r_max;
        Vector2d v_max;
        MatrixXd r_persbest;
        Vector2d r_globalbest;
        void update_vel(int i, params constants);
        void update_pos();
        void calc_globalbest();
        void save(ofstream& myfile);
};

Swarm::Swarm(int N, Vector2d r_max, Vector2d v_max):
    N(N),
    r(N, N),
    r_max(r_max),
    v_max(v_max),
    r_persbest(N, N)
{
    for (int i = 0; i < 2; i++){
        r(i, 0) = (double(rand())/double(RAND_MAX)* (r_max[1] - r_max[0])) + r_max[0]; // draw a number between 0 and 1, then move it into the intervall
        r(i, 1) = (double(rand())/double(RAND_MAX)* (r_max[1] - r_max[0])) + r_max[0]; // this is for the y value
        v(i, 0) = (double(rand())/double(RAND_MAX)* (v_max[1] - v_max[0])) + v_max[0]; // draw a random number between 0 and one and move it into the velocity intervall
        v(i, 1) = (double(rand())/double(RAND_MAX)* (v_max[1] - v_max[0])) + v_max[0];
    }
    r_persbest = r;
    r_globalbest[0] = r(0, 0); // initialize r_globalbest with some random value
    r_globalbest[1] = r(0, 1);
    calc_globalbest(); // now calculate the actual value
}

void Swarm::find_minimum(int Iterations, params constants, string path){
    ofstream myfile1;
    myfile1.open(path);
    for (int i = 0; i < Iterations; i++){
        update_pos();
        for (int j = 0; j < N; j++){
            update_vel(j, constants);
        }
        calc_globalbest();
        save(myfile1);
    }
    myfile1.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Private functions

///////////////////////////////////////////////////////////////////////////////////////////////////////////

void Swarm::update_pos(){ // update the particle position
    r = r + v;
}

void Swarm::update_vel(int i, params constants){ // calculate the new velocities 
    v = constants.omega*v + constants.c_1 * r(i,0) * (r_persbest - r) + constants.c_2 * r(i,1) * (r_globalbest - r);
}

void Swarm::calc_globalbest(){ // calculates the global best value 
    double globalfind;
    for (int i = 0; i < N; i++){
        globalfind = f(r(i, 0), r(i, 1));
        if(globalfind < f(r_globalbest[0], r_globalbest[1])){
            r_globalbest[0] = r(i, 0);
            r_globalbest[1] = r(i, 1);
        }
    }
}

void Swarm::save(ofstream& myfile){
    for (int i = 0; i < N; i++){
        myfile << r
    } 
}

int main() {

    return 0;
};