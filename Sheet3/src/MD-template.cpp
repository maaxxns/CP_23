#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <fstream>

using namespace std;
using namespace Eigen;
using std::ofstream;

// =================================================================================================
//                      PROGRAM STRUCTURE
//
//                          ===========
//                          | Dataset |
//                          ===========
//                              |
//                           ========
//                           | Data |
//                           ========
//                              |
// ==============             ======                =============
// | Thermostat | ----------- | MD | -------------- | Potential |
// ==============             ======                =============
//                              |
//                           ========
//                           | main |
//                           ========
//
// - The MD class contains the primary logic of the algorithm
// - Thermostat and Potential class are separated from it to allow different
// thermostats and potentials can be implemented by inheritance and used flexibly.
// can be used
// - The Data class stores the data stored in the MD simulation and
// takes care of the storage
// - Dataset is a data set consisting of time, temperature, ...
// - Data holds several datasets and some more data that are not time-resolved
// stored
// - Instead of using getter and setter the members of Data
// and Dataset are public, since they are simple data containers.
// - main() calls MD with the parameters needed for the task parts
//
// Notes on the vectors used:
// - For performance reasons, we use Vector2d instead of VectorXd.
// However, this makes the possible transition to 3d states more cumbersome than with VectorXd
// - For lists of data std::vector is used.
// =================================================================================================


// ================================ Potential-class ================================================

// Virtual class from which concrete potentials can be inherited
// (Here only Lennard-Jones necessary, but so you can quickly implement other potentials).

struct parameter {
    Vector2d r_1;
    Vector2d r_minus_1_1;
    Vector2d v_1;
    double mass = 1;
};

double Distance(double r_1_x, double r_1_y, double r_2_x, double r_2_y){ //Distance between r_1 and r_2
    return sqrt(pow((r_1_x-r_2_x),2) + pow((r_1_y-r_2_y),2)); 
}


parameter verlet(Vector2d F, parameter verlet_parameter,int i, double h, double L){
    Vector2d a_1;
    Vector2d r_plus_1_1; // r_{1, n+1}
    
    for (int j = 0; j<2; j=j+1){
    a_1[j] = F[j]; // calculate the acceleration for r_1
    
    if(i == 0){ // start parameters r_{n-1} 
        verlet_parameter.r_minus_1_1[j] = verlet_parameter.r_1[j] - verlet_parameter.v_1[j]*h + 1./2. *a_1[j] * pow(h,2);
    }

    r_plus_1_1[j] = 2.*verlet_parameter.r_1[j] - verlet_parameter.r_minus_1_1[j] + a_1[j] * pow(h,2); // next position
    
    if(r_plus_1_1[j] > L){ // periodicall boundaries
        r_plus_1_1[j] = r_plus_1_1[j] - L*floor(r_plus_1_1[j]/L);
    }
    verlet_parameter.v_1[j] = (r_plus_1_1[j] - verlet_parameter.r_minus_1_1[j]); 

    verlet_parameter.r_minus_1_1[j] = verlet_parameter.r_1[j]; // make the n the n-1 step
    verlet_parameter.r_1[j] = r_plus_1_1[j];    // make the n+1 the n step so basically we shift everything by one
    }
    return verlet_parameter;
}

class Potential
{
    public:
        virtual double      V               ( double r2 ) const = 0;  // Virtual function
        virtual Vector2d    F               ( Vector2d r ) const = 0; // Virtual function
};

class PotentialLJ: public Potential
{
    public:
        double              V               ( double r2 ) const;  // Overwrites virtual function
        Vector2d            F               ( Vector2d r ) const; // Overwrites virtual function
};

// For the potential, the square of the vector length is sufficient, which saves a root calculation.
double PotentialLJ::V ( double r2 ) const
{
    return (4.*(pow(1./r2, 12) - pow(1./r2, 6)));
}

Vector2d PotentialLJ::F ( Vector2d r ) const
{
    Vector2d potential;
    for (int i = 0; i<2; i++){
        potential = {(4.*(pow(1./ r[0], 12) - pow(1./r[0], 6))), (4.*(pow(1./r[1], 12) - pow(1./r[1], 6)))};
    }
    return potential;
}

// ------------------------------ End of Potential-class -------------------------------------------
// ================================ Thermostat class ===============================================

// Virtual class from which concrete thermostats can be inherited
class Thermostat
{
    public:
        virtual void        rescale         ( vector<Vector2d>& v, double T ) const = 0;
};

// No thermostat
class NoThermostat: public Thermostat
{
    public:
        void                rescale         ( vector<Vector2d>& v, double T ) const {} // does nothing
};

// Isokinetic thermostat for task d)
class IsokinThermostat: public Thermostat
{
    public:
        void                rescale         ( vector<Vector2d>& v, double T ) const;
};

void IsokinThermostat::rescale( vector<Vector2d>& v, double T ) const
{
    double T_t;
    int N_f = 2*v.size() - 2; // degrees of freedom
    for (int i = 0; i < v.size();){
        T_t += 2./N_f * (pow(v[i][0], 2) + pow(v[i][1],2));
        v[i] = v[i] * pow(T/T_t, 2);
    }
}

// ------------------------------ End of Thermostat class ------------------------------------------
// ================================ Data-Structs ===============================================

// Data set for time resolved data
// (structs basically the same as class, but all members are public by default)
struct Dataset
{
    double                  t, T, Ekin, Epot;
    Vector2d                vS;
};

// Return data of the MD simulation
// Data data(n) constructor; reserves memory and fills pair correlation function with 0s
struct Data
{
    vector<Dataset> datasets; // Time-resolved datasets.
    vector<double> rBin, g;   // Averaged pair correlation function
    vector<Vector2d> r;       // snapshot of the final position
                              // For task e) it may be useful to use r instead
                              // in the time-resolved datasets instead

                            Data            ( uint n, uint numBins, double binSize );
    void                    save            ( const string& filenameSets,
                                              const string& filenameG,
                                              const string& filenameR ) const;
};

Data::Data( uint n, uint numBins, double binSize ):
    datasets( n ),  // Initializer list, because it calls constructors of the members
    rBin( numBins ),
    g( numBins, 0. ),
    r( 0 )
{
}

void Data::save ( const string& filenameSets, const string& filenameG, const string& filenameR ) const
{   
    const string filenames[3] = {filenameSets, filenameG, filenameR}; // I dont really get this but anyway
    for (int j = 0; j < 3 ; j++){
        ofstream Sets(filenames[j]);
        for(int i = 0; i<=datasets.size(); i++){
            if(i == 0){
                Sets << "#t, T, Ekin, Epot, vSx, vSy" << endl; 
            }
        Sets << datasets[i].t << ',' << datasets[i].T << ',' << datasets[i].Ekin << ',' << datasets[i].Epot << ',' << datasets[i].vS[0] << "," << datasets[i].vS[1] << endl;
    }
    }
}

// ------------------------------ End of Data-Structs ------------------------------------------
// ================================ MD-Class ===============================================

class MD
{
    public:
                            MD              ( double L, uint N, uint particlesPerRow, double T,
                                            Potential& potential, Thermostat& thermostat,
                                            uint numBins = 1000 );

        void                equilibrate     ( const double dt, const unsigned int n );
        Data                measure         ( const double dt, const unsigned int n );

    private:
        vector<Vector2d>    r, v;
        double              L;
        uint                N;
        Potential&          potential;
        Thermostat&         thermostat;
        double              t = 0.;

        uint                numBins;
        double              binSize;

        // Particles are moved in box [0,L]x[0,L].
        void                centerParticles ();

        // Calculations of important measured variables
        double              calcT           () const;
        double              calcEkin        () const;
        double              calcEpot        () const;
        Vector2d            calcvS          () const;
        Dataset             calcDataset     () const;

        // Calculation of the acceleration
        // To avoid redundant calculations, it may be useful to update the histogram
        // when calculating the accelerations, so it is passed here as a reference.
        vector<Vector2d>    calcAcc         ( vector<double>& hist ) const;

        // Calculation of the distance vector between particle r[i] and closest mirror particle of r[j].
        Vector2d            calcDistanceVec ( uint i, uint j ) const;
};

// Initialization of the system via constructor
MD::MD( double L, uint N, uint particlesPerRow, double T,
        Potential& potential, Thermostat& thermostat,
        uint numBins ):
        L(L), N(N), potential(potential), thermostat(thermostat),
        numBins(numBins), binSize(L/N),
        r(N), v(N)
{
        // Initialize random number generator
        mt19937 rnd;
        uniform_real_distribution<> dist(-1, 1);

        // Initialize positions

        // calculate the distance between each particle
        const double distance = sqrt(2.0 * L * L / (N* sin(2.0 * M_PI / N)));
        
        // fill the box
        for (int i = 0; i < particlesPerRow; i++) {
            for (int j = 0; j < particlesPerRow; j++) {
                double x = (i + 0.5) * distance; // we dnt want particles on the edge of the box 
                double y = (j + 0.5) * distance; // because of the way we defined the equal distance
                if(x < L && y < L){ // test the boundary conditions
                r[i] = {x, y};
                }
            }
        }

        // Initialize velocities with random numbers
        for (int i = 0; i < N; i++) {
            for (int j=0; j<2; j++) {
                v[i] = {dist(rnd), dist(rnd)}; // this just gives us random velocities from -1 - 1
            }
        }
        // calculate the velocity of the center of mass which is a simple sum over all velocities
        double vxcm, vycm; // define velocity of the center of mass

        Vector2d vcm = calcvS();

        // Rescale velocities to have a non moving center of mass
        for (int i = 0; i < N; i++) {
            v[i][0] -= vcm[0];
            v[i][1] -= vcm[1];
        }

        // Now rescale velocities for given T 


}

// Integration without data acquisition for pure equilibration
void MD::equilibrate ( const double dt, const unsigned int n )
{
    /*TODO*/
    vector<int> vec(10);
    for ( int i: vec )
    {
        cout << i << "\t";
    }
}

Data MD::measure ( const double dt, const unsigned int t_end )
{
    Data data(int(t_end/dt), 2, 2.); // number bins??
    for (int steps = 0; steps < int(t_end/dt); steps++){ // the actual time steps
        Dataset dataset;
        vector<Vector2d> force_i; // contains the force on every particle 
        parameter pos;
        vector<Vector2d> r_minus1(N); //
        if(steps == 0){ // initial parameters
            dataset.Ekin = calcEpot(); 
            dataset.Ekin = calcEkin(); 
            dataset.T = calcT();
            dataset.t = steps*dt; 
            dataset.vS = calcvS();
            data.datasets[steps] = dataset;
            }
        for (int i = 0; i < N; i++){ // sum over all particles
            Vector2d force_ij;
            force_ij[0] = 0.;
            force_ij[1] = 0.;
            for(int k = 0; k < N; k++){ // sum over all particles but one 
                if(i != k){
                    for (int l = -1; l < 1; l++){
                        Vector2d L_vec = {l*L, l*L};
                        Vector2d r_dist = calcDistanceVec(i, k);
                        double r1_2 = Distance(r_dist[0], r_dist[1], l * L, l * L);
                        if(r1_2 < (L/2.) and r1_2 != 0){ // cutoff r1_2 < (L/2.)
                            force_ij += -(calcDistanceVec(i, k) + L_vec)/(r1_2) * (potential.V(r1_2) - potential.V(L/2)); // force of particle j on particle i with boundary conditions
                        }else {force_ij += Vector2d {0,0};} // coutoff gives us force = 0
                        force_ij += force_ij;  // add force to total force of partcle i// segmentation error
                        }
                    }
                    force_ij += force_ij;
                }
                force_i.push_back(force_ij);
            }
        for (int i = 0; i < N; i++){ // calculate new position of all particles
            if(i != 0){
                pos.r_minus_1_1 = r_minus1[i]; // print in the last position of particle i 
                } 
            pos.v_1 = v[i];
            pos.r_1 = r[i];
            pos = verlet(force_i[i], pos, steps, dt, L); // pos also contains last position
            r[i] = pos.r_1;
            r_minus1[i] << pos.r_minus_1_1; // save last position of particle i
            v[i] = pos.v_1; // save the new positon in the r vector
        }
        dataset.Ekin = calcEpot(); 
        dataset.Ekin = calcEkin(); 
        dataset.T = calcT(); // 
        dataset.t = steps*dt; 
        dataset.vS = calcvS();
        data.datasets[steps] = dataset;
    }   
    return data;
}

void MD::centerParticles()
{
    for (int i = 0; i < N; i++){
        for (int j = 0; j < 2 ; j++){
            if (r[i][j] < 0 or r[i][j] > L){
                r[i][j] = r[i][j] - L*floor(r[i][j]/L); // recenter particles in box
            }
        }
    }
}

double MD::calcT() const
{
    int N_f = 2 * N - 2;// degrees of freedom
    double T = 2./N_f * calcEkin(); // this is T*k_B in units of epsilon
    return T;
}

double MD::calcEkin() const
{
    double Ekin = 0;
    for (int i = 0; i < N; i++){
        double power = (pow(abs(v[i][0]), 2) + pow(abs(v[i][1]), 2)); // power = v^2 as we would need to overload the function pow otherwise
        Ekin += 1./2. * power;
    } 
    return Ekin;
}

double MD::calcEpot() const
{   
    double Epot;
    for (int i = 0; i < N; i++){
    double r2 = (pow(r[i][0], 2) + pow(r[i][1],2));
    Epot += potential.V(r2);
    }
    return Epot;
}

Vector2d MD::calcvS() const // if this really is center of mass this should be called calcvCM()
{
    double vxcm, vycm;
    for (int i = 0; i < N; i++) {
    vxcm = vxcm + v[i][0];
    vycm = vycm + v[i][1];
    }
    return {vxcm, vycm};
}

//Dataset MD::calcDataset() const
//{
//    /*TODO*/
//}

Vector2d MD::calcDistanceVec( uint i, uint j ) const
{
    return r[i] - r[j];
}

//vector<Vector2d> MD::calcAcc( vector<double>& hist ) const
//{
//    /*TODO*/
//}

// ------------------------------ End of MD-class ------------------------------------------


int main(void)
{
    PotentialLJ      LJ;
    NoThermostat     noThermo;
    IsokinThermostat isoThermo;

    const int n                 = 2;
    const uint N                = n*n;
    const double L              = 2*n;
    const double distance = sqrt(2.0 * L * L / (N* sin(2.0 * M_PI / N)));
    const uint partPerRow       = ceil(L/distance);
    const int numBins           = 3; // ???????????????????????????????????????????????????

    // b) Equilibration test
    {
        const double T          = 1;
        const double dt         = 0.01;
        const uint t_end        = int(T/dt);// number of steps?? or what??

        MD md( L, N, partPerRow, T, LJ, noThermo, numBins );
        md.measure( dt, t_end ).save( "b)set.dat", "b)g.dat", "b)r.dat" );
    }

    // c) Pair correlation function
    //string TstringVec[3] = { "0.01", "1", "100" };
    //for ( auto& Tstring: TstringVec )
    //{
    //    const double T          = stod(Tstring);
    //    const double dt         = /*TODO*/
    //    const uint equiSteps    = /*TODO*/
    //    const uint steps        = /*TODO*/
//
    //    MD md( L, N, partPerRow, T, LJ, noThermo, numBins );
    //    md.equilibrate( dt, equiSteps );
    //    md.measure( dt, steps ).save( "c)set" + Tstring + ".dat", "c)g" + Tstring + ".dat", "c)r" + Tstring + ".dat" );
    //}
//
    //// d) Thermostat
    //{
    //    /*TODO*/
    //}


    return 0;
}

