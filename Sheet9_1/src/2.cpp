#include <iostream>
#include <math.h>
#include <fstream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
using std::ofstream; 

class Swarm{

    public:
        Swarm(int N, Vector2d r_max, Vector2d v_max); // constructor for N particles in the box r_i^2
    private:
        MatrixXd r; // particle positions for N particles
        MatrixXd v; // particle velocity for N particles
        Vector2d r_max;
        Vector2d v_max;
        MatrixXd r_persbest;
        Vector2d r_globalbest;
        void update_vel();
        void update_pos();
};

Swarm::Swarm(int N, Vector2d r_max, Vector2d v_max):
    r(N, N),
    r_max(r_max),
    v_max(v_max)
{
    for (int i = 0; i < 2; i++){
        r(i, 0) = (double(rand())/double(RAND_MAX)* (r_max[1] - r_max[0])) + r_max[0]; // draw a number between 0 and 1, then move it into the intervall
        r(i, 1) = (double(rand())/double(RAND_MAX)* (r_max[1] - r_max[0])) + r_max[0]; // this is for the y value
        v(i, 0) = (double(rand())/double(RAND_MAX)* (v_max[1] - v_max[0])) + v_max[0]; // draw a random number between 0 and one and move it into the velocity intervall
        v(i, 1) = (double(rand())/double(RAND_MAX)* (v_max[1] - v_max[0])) + v_max[0];
    }
}

void Swarm::update_pos(){
    r = r + v;
}

void Swarm::update_vel(double omega, double c_1, double c_2){
    v = omega*v + c_1 * r_(0) * (r_persbest - r) + c_2 * r_2 * (r_globalbest - r)
}

int main() {

    return 0;
};