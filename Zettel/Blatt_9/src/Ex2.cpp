#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <random>
using namespace std;

class particle_swarm{
    public:
        particle_swarm(vector<double> r_vec,vector<double>, int);
        vector<vector<double>> r ={};
        vector<vector<double>> v ={};
        vector<double> r_global;
        void iterate(int N_iterations, double (*f)(vector<double>&));
    private:
        int N_particles; //nr of particles
        vector<double> r_intervall; //intervall of start velocities
        vector<double> v_intervall; //intervall of start positions

        vector<vector<double>> r_pers;
        

        double w = 0.8;
        double c1 = 0.1;
        double c2 = 0.1;

        void initialise();
        


};

particle_swarm::particle_swarm(vector<double> r_inter, vector<double> v_inter, int N_){
    r_intervall = r_inter;
    v_intervall = v_inter;
    N_particles = N_;
    initialise();

}; 

void particle_swarm::initialise(){
    mt19937 rnd;
    uniform_real_distribution<double> r_dist(r_intervall[0], r_intervall[1]);
    uniform_real_distribution<double> v_dist(v_intervall[0], v_intervall[1]);
    for(int i=0; i<N_particles; i++){
            r.push_back({r_dist(rnd), r_dist(rnd)});
            v.push_back({v_dist(rnd), v_dist(rnd)});
    }
};  


void particle_swarm::iterate(int N_iterations, double (*f)(vector<double>&)){
    mt19937 rnd;
    uniform_real_distribution<double> r12_dist(0, 1);
    double r1, r2;
    ofstream data_global("./build/data_global.csv");
    ofstream data_perso("./build/data_perso.csv");
    ofstream data_r("./build/data_r.csv");
    ofstream data_v("./build/data_v.csv");
    ofstream data_global_r("./build/data_global_r.csv");
    ofstream data_perso_r("./build/data_perso_r.csv");

    r_pers = r; //using first values as pers min
    r_global = r_pers[0]; //first setting global minimum on any value 
    //find right r_global:
    for(int i=0; i<N_particles; i++){
        if(f(r_pers[i]) < f(r_global)){ //if f(r_pers) is smaller than f(r_global): r_global -> r_pers
            r_global = r_pers[i];
        }



        //data saving by just using same "for":
        data_perso << f(r_pers[i]) << ' ';
        data_perso_r << r_pers[i][0] << ',' << r_pers[i][1] << ' ';
        data_r << r[i][0] << ' ' << r[i][1] << ' ';
        data_v << v[i][0] << ' ' << v[i][1] << ' ';

    }
    data_global << f(r_global) << endl;
    data_global_r << r_global[0] << ' ' << r_global[1] << endl;
    data_r << endl;
    data_v << endl;
    data_perso << endl;
    data_perso_r << endl;

    for(int n = 0; n<N_iterations; n++){
        for(int i=0; i<N_particles; i++){
            vector<vector<double>> r_old  = r; //save a copy of r
            r1 = r12_dist(rnd);
            r2 = r12_dist(rnd);
            //r iteration:
            r[i][0] = r_old[i][0] + v[i][0];
            r[i][1] = r_old[i][1] + v[i][1];


            //v iteration:
            v[i][0] = w*v[i][0] + c1*r1*(r_pers[i][0]-r_old[i][0]) + c2*r2*(r_global[0]-r_old[i][0]);
            v[i][1] = w*v[i][1] + c1*r1*(r_pers[i][1]-r_old[i][1]) + c2*r2*(r_global[1]-r_old[i][1]);

            //update the r_pers:
            if(f(r[i])<f(r_pers[i])){ //if f(r) is smaller than f(r_pers): r_pers -> r
                r_pers[i] = r[i];
            }
            //update r_global:
            if(f(r_pers[i]) < f(r_global)){
                r_global = r_pers[i];
            }

            //saving new generated data:
            data_perso << f(r_pers[i]) << ' ';
            data_perso_r << r_pers[i][0] << ',' << r_pers[i][1] << ' ';
            data_r << r[i][0] << ' ' << r[i][1] << ' ';
            data_v << v[i][0] << ' ' << v[i][1] << ' ';
        }
        data_r << endl;
        data_v << endl;
        data_perso << endl;
        data_perso_r << endl;
        data_global << f(r_global) << endl;
        data_global_r << r_global[0] << ' ' << r_global[1] << endl;
    }
    data_global.close();
    data_global_r.close();
    data_perso.close();
    data_perso_r.close();
    data_r.close();
    data_v.close();
};

double func(vector<double>& r){
    return (r[0]-1.9)*(r[0]-1.9)+(r[1]-2.1)*(r[1]-2.1)+2*cos(4*r[0]+2.8)+3*sin(2*r[1]+0.6);
}


int main(){
    int N = 10;
    vector<double> r_intervall = {-5.,5.}; //intervall of start velocities
    vector<double> v_intervall = {-1.,1.};
    particle_swarm swarm1(r_intervall, v_intervall, N);
    swarm1.iterate(100, func);
    cout << "r_global: " << swarm1.r_global[0] << ',' << swarm1.r_global[1] << endl;
    cout << "f(r_global)=" <<  func(swarm1.r_global) << endl;


    return 0;
}