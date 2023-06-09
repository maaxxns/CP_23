#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <fstream>

using namespace std;
using namespace Eigen;
using std::ofstream;

// sources of errors
// 1. distance is calculated wrong
// 2. the potential is wrong or has the wrong sign 
// 3. I dont look correctly in the neighbouring boxes
// 4. the velocity is calculated wrong in verlet
// 5. the boundary condiations a not implemented correctly 
// 6. 


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

double Distance(double r_1_x, double r_1_y, double r_2_x, double r_2_y){ //Distance between r_1 and r_2
    return sqrt(pow((r_1_x-r_2_x),2) + pow((r_1_y-r_2_y),2)); 
}

double length(Vector2d vec){ // returns the length of vector
        double distance = sqrt(pow(vec[0], 2) + pow(vec[1], 2));
    return distance;
}

class Potential
{
    public:
        virtual double      V               ( double r2 ) const = 0;  // Virtual function
        virtual double      V2              ( double r2 ) const = 0;
        virtual Vector2d    F               ( Vector2d r ) const = 0; // Virtual function
};

class PotentialLJ: public Potential
{
    public:
        double              V               ( double r2 ) const;  // Overwrites virtual function
        double              V2              ( double r2 ) const; 
        Vector2d            F               ( Vector2d r ) const; // Overwrites virtual function
};

// For the potential, the square of the vector length is sufficient, which saves a root calculation.
double PotentialLJ::V ( double r2 ) const
{
    return 48./r2 * (pow(1./r2, 12) - 1./2. * pow(1./r2, 6));
    //return (4.*(pow(1./r2, 12) - pow(1./r2, 6)));
}

double PotentialLJ::V2 ( double r2 ) const
{
    return 4.* (pow(1./r2, 12) - pow(1./r2, 6));
    //return (4.*(pow(1./r2, 12) - pow(1./r2, 6)));
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
    for (int i = 0; i < v.size(); i++){
        T_t += 1./N_f * (pow(v[i][0], 2) + pow(v[i][1],2));
    }
    for (int i = 0; i < v.size(); i++){
        v[i] = v[i] * sqrt(T/T_t);
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
    g( numBins),
    r( 0 )
{
}

void Data::save ( const string& filenameSets, const string& filenameG, const string& filenameR ) const
{   
    const string filenames[3] =  {"bin/" + filenameSets,"bin/" + filenameG,"bin/" + filenameR}; // I dont really get this but anyway
    ofstream myfile1;
    ofstream myfile2;
    for (int j = 0; j < 3 ; j++){
        if(j == 0){
        myfile1.open(filenames[j]);
            for(int i = 0; i<=datasets.size() -1; i++){
                if(i == 0){
                    myfile1 << "#t, T, Ekin, Epot, vSx, vSy" << endl; 
                }
                myfile1 << datasets[i].t << ',' << datasets[i].T << ',' << datasets[i].Ekin << ',' << datasets[i].Epot << ',' << datasets[i].vS[0] << "," << datasets[i].vS[1] << endl;
            }
            myfile1.close();
        }   
        if(j == 1){
            myfile2.open(filenames[j]);
            for (int i = 0; i < g.size(); i++){
                myfile2 << g[i] << endl;
            }
            myfile2.close();
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

        void                equilibrate     ( const double dt, const unsigned int n, double T );
        Data                measure         ( const double dt, const unsigned int n, double T );
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
        vector<double>      paircorrfun     () const;

        // Calculation of the acceleration
        vector<Vector2d>    calcAcc         () const;

        // Calculation of the distance vector between particle r[i] and closest mirror particle of r[j].
        Vector2d            calcDistanceVec ( uint i, uint j ) const;
};

// Initialization of the system via constructor
MD::MD( double L, uint N, uint particlesPerRow, double T,
        Potential& potential, Thermostat& thermostat,
        uint numBins ):
        L(L), N(N), potential(potential), thermostat(thermostat),
        numBins(numBins), binSize(1),
        r(N), v(N)
{
        // Initialize random number generator
        mt19937 rnd;
        uniform_real_distribution<> dist(-1, 1);

        // bin size
        binSize = L/(2*numBins);
        // Initialize positions

        // calculate the distance between each particle
        const double distance = sqrt(2.0 * L * L / (N* sin(2.0 * M_PI / N)));
        
        // fill the box
        int m = 0;
        for (int i = 0; i < particlesPerRow; i++) {
            for (int j = 0; j < particlesPerRow; j++) {
                double x = (i + 0.5) * distance; // we dnt want particles on the edge of the box 
                double y = (j + 0.5) * distance; // because of the way we defined the equal distance
                x = x - L * floor(x / L);
                y = y - L * floor(y / L);
                if(x < L && y < L && m < N){ // test the boundary conditions
                    r[m] = {x, y};
                    m++;
                }
            }
        }

        // Initialize velocities with random numbers
        for (int i = 0; i < N; i++) {
                v[i] = {dist(rnd), dist(rnd)}; // this just gives us random velocities from -1 - 1
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
        thermostat.rescale(v, T);

        // some safety
        if(calcvS().squaredNorm() > 0.000001){
            cout << "The Center of mass velocity is non zero!" << endl;
            cout << calcvS().squaredNorm() << endl;
        }
        if(calcT() < 1.){
            cout << "Temperature below 1" << endl;
            cout << calcT() << endl;
        }

}

// Integration without data acquisition for pure equilibration
void MD::equilibrate ( const double dt, const unsigned int t_end, double T )
{
    for (int steps = 0; steps < t_end/dt; steps++){ // the actual time steps
    vector<Vector2d> force_i; // contains the force on every particle 
    vector<Vector2d> force_i_plus_1; // contains the force on every particle at postion r_{n+1}
        { // Velocity verlet
            force_i = calcAcc(); // calc the force on every particle
            for (int i = 0; i < N; i++){ // all particles
                for (int j = 0; j < 2; j++){ // all dimensions
                    r[i][j] = r[i][j] + v[i][j] * dt + 1./2. * force_i[i][j] * pow(dt,2); // get all the new positions r
                    centerParticles();
                }
            }
            force_i_plus_1 = calcAcc(); // get the forces at the new positions
            for (int i = 0; i < N; i++){ // all particles
                for(int j = 0; j < 2; j++){ // all dimensions
                    v[i][j] = v[i][j] + 1./2. * ( force_i_plus_1[i][j] + force_i[i][j]) * dt; // get v_{n+1}
                }
            }
        } // end of velocity verlet
        thermostat.rescale(v, t);
    }
    cout << "system is kind of in equilibrium" << endl;
}

Data MD::measure ( const double dt, const unsigned int t_end, double T )
{
    equilibrate(dt, int(t_end/5), T);
    Data data(t_end/dt, numBins, binSize); // number bins??
    vector<Vector2d> r_minus1(N); //
    vector<double> p_l; // all the pair found by the pair correlation over t
    for (int steps = 0; steps < t_end/dt; steps++){ // the actual time steps
        Dataset dataset;
        vector<Vector2d> force_i; // contains the force on every particle 
        vector<Vector2d> force_i_plus_1; // contains the force on every particle at postion r_{n+1}
        if(steps == 0){ // initial parameters
            dataset.Epot = calcEpot(); 
            dataset.Ekin = calcEkin(); 
            dataset.T = calcT();
            dataset.t = steps*dt; 
            dataset.vS = calcvS();
            data.datasets[steps] = dataset;
            }
            { // Velocity verlet
                force_i = calcAcc(); // calc the force on every particle
                for (int i = 0; i < N; i++){ // all particles
                    for (int j = 0; j < 2; j++){ // all dimensions
                        r[i][j] = r[i][j] + v[i][j] * dt + 1./2. * force_i[i][j] * pow(dt,2); // get all the new positions r
                        centerParticles();
                    }
                }
                force_i_plus_1 = calcAcc(); // get the forces at the new positions
                for (int i = 0; i < N; i++){ // all particles
                    for(int j = 0; j < 2; j++){ // all dimensions
                        v[i][j] = v[i][j] + 1./2. * ( force_i_plus_1[i][j] + force_i[i][j]) * dt; // get v_{n+1}
                    }
                }
            } // end of velocity verlet
        for (int i = 0; i < numBins; i++){
            data.g[i] = data.g[i] + paircorrfun()[i]/double(N * (L/numBins) * t_end * (N / (L * L))); // add to the overall count and average
        }
        dataset.Epot = calcEpot(); 
        dataset.Ekin = calcEkin(); 
        dataset.T = calcT(); // 
        dataset.t = steps*dt; 
        dataset.vS = calcvS();
        data.datasets[steps] = dataset;
        thermostat.rescale(v, T);
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

vector<double> MD::paircorrfun() const // calculates the number of pair at time t
{
    
    double r_ik;
    vector<double> pairs_in_bin(numBins);
    for (int i = 0; i < N; i++){ // all particles i
        for (int k = 0; k < N; k++){ // all particles k
            for (int j = 0; j < 2; j++){ // all dimensions
                if(r[i][j] > 0 && r[i][j] < L/2.){ // boundary conditions
                    if(i != k){ // dont overcount particle i = k
                        r_ik = length(calcDistanceVec(i, k));
                        for(int l = 1; l <= numBins; l++){
                                if((l-1.)*binSize < r_ik && r_ik < l*binSize){
                                    pairs_in_bin[l - 1] =+ 1; // if we found a pair we add one to that positon
                                }
                            }
                        }
                    }
                }
            }
        }
    return pairs_in_bin;
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
        for (int k = 0; k < N; k++){
            if(i != k){
                for (int n_x = -1; n_x < 2; n_x++){
                    for (int n_y = -1; n_y < 2; n_y++){
                        Vector2d L_vec = {n_x*L, n_y*L}; // this vector brings us in the neighbouring boxes
                        Vector2d r_dist = calcDistanceVec(i, k);
                        double r1_2 = length(r_dist + L_vec); // |r_ij + nL|
                            if(r1_2 < (L/2.) and r1_2 != 0){
                                Epot += potential.V2(r1_2);
                            }
                    }
                }
            }
        }
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
    return {1./N *vxcm, 1./N *vycm};
}

//Dataset MD::calcDataset() const
//{
//    /*TODO*/
//}

Vector2d MD::calcDistanceVec( uint i, uint j ) const
{
    return (r[i] - r[j]);
}

vector<Vector2d> MD::calcAcc() const
{
    vector<Vector2d> force_i;
    for (int i = 0; i < N; i++){ // sum over all particles
    Vector2d force_ij;
    force_ij[0] = 0.;
    force_ij[1] = 0.;
    for(int k = 0; k < N; k++){ // sum over all particles but one 
        if(i != k){
            for (int n_x = -1; n_x < 1; n_x++){ // lets look in all the virtual boxes around the normal box
                for (int n_y = -1; n_y < 1; n_y++){
                    Vector2d L_vec = {n_x*L, n_y*L}; // this vector brings us in the neighbouring boxes
                    Vector2d r_dist = calcDistanceVec(i, k);
                    double r1_2 = length(r_dist + L_vec); // |r_ij + nL|
                        if(r1_2 < (L/2.) and r1_2 != 0){ // cutoff r1_2 < (L/2.)
                            force_ij += (calcDistanceVec(i, k) + L_vec)/(r1_2) * (potential.V(r1_2) - potential.V(L/2.)); // force of particle j on particle i with boundary conditions
                        }else {force_ij += Vector2d {0,0};} // cutoff gives us force = 0
                    }
                }
            }
        }
        force_i.push_back(force_ij); // add forces of particle i two the overall force vector
    }
    return force_i;
}

// ------------------------------ End of MD-class ------------------------------------------


int main(void)
{
    PotentialLJ      LJ;
    NoThermostat     noThermo;
    IsokinThermostat isoThermo;

    const int n                 = 16;
    const uint N                = n*n;
    const double L              = 2*n;
    const double distance = sqrt(2.0 * L * L / (N* sin(2.0 * M_PI / N)));
    const uint partPerRow       = n;
    const int numBins           = 100; // ???????????????????????????????????????????????????


    { // test with n = 4
    string TstringVec[3] = { "001", "1", "100" };
    for ( auto& Tstring: TstringVec )
    {
        const double T          = 1;
        const double dt         = 0.001;
        const uint t_end        = 10;// number of steps?? or what??

        MD md4( L, N, partPerRow, T, LJ, isoThermo, numBins );
        md4.measure( dt, t_end , T).save( "bN" + to_string(N) + to_string(T) + ".dat","bg_N" + to_string(N) + to_string(T) +".dat", "dummy" );
    }
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

