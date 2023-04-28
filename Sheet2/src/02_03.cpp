#include <iostream>
#include <math.h>
#include <fstream>
#include <ctime> // time_t

using namespace std;

using std::string;
using std::ofstream;
using std::fstream;
int i = 0;

struct parameter {
    double r_1[2];
    double r_2[2];
    double r_minus_1_1[2];
    double r_minus_1_2[2];
    double v_1[2];
    double v_2[2];
    double mass1;
    double mass2;
    int T;
};


double Newton_Gravitation(double r, double R, double mass){
    if(r==0) {
        return 0;
    }
    return (-(r/pow(R,3) * mass));
}

void csv_print(double trajectory[4], string filename){
    ofstream myfile;
    const string path = ("bin/" + filename + ".csv"); // I would rather have one big array which I write into the file but I ran into porblems with that. As the array became to big and I got an overflow
    myfile.open(path.c_str(), fstream::app);
    for (int j=0; j<4; j++){
    myfile << trajectory[j]; // write the euler approximation into a .csv file
    if (j<3){
        myfile << ", ";
    };
    }
    myfile << endl;
    myfile.close(); 

}

double Distance(double r_1_x, double r_1_y, double r_2_x, double r_2_y){ //Distance between r_1 and r_2
    return pow(pow((r_1_x-r_2_x),2) + pow((r_1_y-r_2_y),2),1./2.); 
}

parameter euler(double (*F)(double, double, double), parameter euler_parameter, double h, double trajectory[4]){ //F is the Force e.g. Newton Gravitation, r_n the Value of the position, v_n the velocity, h the stepwidth wich should be in the same Order of magnitude as r_n, T is the maximum time to which euler integrates, trajectory saves the r_n states to a trajector array of desired size
    double a_1[2]; //acce of mass 1
    double a_2[2]; //acce of mass 2
    double v_n_1[2]; // next velo of mass 1
    double v_n_2[2];    // next velo of mass2
    double R = Distance(euler_parameter.r_1[0], euler_parameter.r_1[1], euler_parameter.r_2[0], euler_parameter.r_2[1]);
    trajectory[0] = euler_parameter.r_1[0];
    trajectory[1] = euler_parameter.r_1[1];
    trajectory[2] = euler_parameter.r_2[0];
    trajectory[3] = euler_parameter.r_2[1];
    csv_print(trajectory, "euler");
    for (int j = 0; j<2; j=j+1){
    a_1[j] = F(euler_parameter.r_1[j], R, euler_parameter.mass2); // calculate the acceleration for r_1
    a_2[j] = F(euler_parameter.r_2[j], R, euler_parameter.mass1); // calculate the acceleration for r_2

    v_n_1[j] = euler_parameter.v_1[j] + a_1[j]*h; // next velocity
    v_n_2[j] = euler_parameter.v_2[j] + a_2[j]*h; // next velocity

    euler_parameter.r_1[j] = euler_parameter.r_1[j] + euler_parameter.v_1[j]*h; //next position
    euler_parameter.r_2[j] = euler_parameter.r_2[j] + euler_parameter.v_2[j]*h; //next position

    euler_parameter.v_1[j] = v_n_1[j]; // save next velocity into struct
    euler_parameter.v_2[j] = v_n_2[j];
    }
    
    i = i + 1; 
    if (i>int(euler_parameter.T/abs(h))){
        return euler_parameter;
    }
    return euler(F, euler_parameter, h, trajectory);
}

parameter verlet(double (*F)(double, double, double), parameter verlet_parameter, double h, double trajectory[4]){
    double R = Distance(verlet_parameter.r_1[0], verlet_parameter.r_1[1], verlet_parameter.r_2[0], verlet_parameter.r_2[1]);
    double a_1[2];
    double a_2[2];
    double r_plus_1_1[2]; // r_{1, n+1}
    double r_plus_1_2[2]; // r_{2, n+1}
    
    trajectory[0] = verlet_parameter.r_1[0];
    trajectory[1] = verlet_parameter.r_1[1];
    trajectory[2] = verlet_parameter.r_2[0];
    trajectory[3] = verlet_parameter.r_2[1];
    csv_print(trajectory, "verlet");

    for (int j = 0; j<2; j=j+1){
    a_1[j] = F(verlet_parameter.r_1[j], R, verlet_parameter.mass2); // calculate the acceleration for r_1
    a_2[j] = F(verlet_parameter.r_2[j], R, verlet_parameter.mass1); // calculate the acceleration for r_2
    
    if(i == 0){ // start parameters r_{n-1} 
        verlet_parameter.r_minus_1_1[j] = verlet_parameter.r_1[j] - verlet_parameter.v_1[j]*h + 1./2. *a_1[j] * pow(h,2);
        verlet_parameter.r_minus_1_2[j] = verlet_parameter.r_2[j] - verlet_parameter.v_2[j]*h + 1./2. *a_2[j] * pow(h,2);
    }

    r_plus_1_1[j] = 2.*verlet_parameter.r_1[j] - verlet_parameter.r_minus_1_1[j] + a_1[j] * pow(h,2); // next position
    r_plus_1_2[j] = 2.*verlet_parameter.r_2[j] - verlet_parameter.r_minus_1_2[j] + a_2[j] * pow(h,2); 
    
    //verlet_parameter.v_1[j] = (r_plus_1_1[j] - verlet_parameter.r_minus_1_1[j]); // I dont even need these but who knows, maybe they come in handy later
    //verlet_parameter.v_2[j] = (r_plus_1_2[j] - verlet_parameter.r_minus_1_2[j]);

    verlet_parameter.r_minus_1_1[j] = verlet_parameter.r_1[j]; // make the n the n-1 step
    verlet_parameter.r_1[j] = r_plus_1_1[j];    // make the n+1 the n step so basically we shift everything by one
    verlet_parameter.r_minus_1_2[j] = verlet_parameter.r_2[j]; // same here just for r_2
    verlet_parameter.r_2[j] = r_plus_1_2[j];
    }
    i = i + 1;

    if (i>int(verlet_parameter.T/abs(h)) - 1){
        return verlet_parameter;
    }
    return verlet(F, verlet_parameter, h, trajectory);
}

int main(){
        //initialize the parameter for the euler problem
    parameter Parameter;
    parameter parameter_end; // output off the functions
    Parameter.T = 1000;
    Parameter.mass1 = 1.;
    Parameter.mass2 = 2.;
    Parameter.r_1[0] = 0.; 
    Parameter.r_1[1] = 1.;
    Parameter.r_2[0] = 0.;
    Parameter.r_2[1] = -0.5;
    Parameter.v_1[0] = 0.8;
    Parameter.v_1[1] = 0.;
    Parameter.v_2[0] = -0.4;
    Parameter.v_2[1] = 0.;
    
    double h = 0.01;
    double trajectory[4]; // Trajectory should fit as many 
    remove("bin/euler.csv");  // I use fstream::add to write into the csv file so I have to delet the old csv when I start the script again
    time_t begin,end;
    time (&begin);



    
    parameter_end = euler(Newton_Gravitation, Parameter, h, trajectory); // The callstack of recursive function grows too large to process and causes a segmentation error
    // I have to find a way to stop the recursion at some point and call it again.






    time (&end);
    double runtime_euler = difftime(end,begin);
    cout << "The euler scheme took " << runtime_euler << " seconds to process" << endl;
    i = 0;
    parameter_end = euler(Newton_Gravitation, parameter_end, -h, trajectory); // Back integration

    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    
    double trajectory1[4];
    remove("bin/verlet.csv");   // I use fstream::add to write into the csv file so I have to delet the old csv when I start the script again
    i = 0;
    time (&begin);
    parameter_end = verlet(Newton_Gravitation, Parameter, h, trajectory1);
    time (&end);
    double runtime_verlet = difftime(end,begin);
    cout << "The verlet scheme took " << runtime_verlet<< " seconds to process" << endl;
    i = 0;
    parameter_end = verlet(Newton_Gravitation, parameter_end, -h, trajectory1); // Back integration
    ///////////////////////////////////////////////////////////////////////////////////////////////////////


    return 0;
}