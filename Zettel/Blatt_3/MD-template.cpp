#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <numeric>
#include <fstream>
#include <random>

using namespace std;
using namespace Eigen;

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
// own functions:



// ================================ Potential-class ================================================

// Virtual class from which concrete potentials can be inherited
// (Here only Lennard-Jones necessary, but so you can quickly implement other potentials).
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
    return 4*(pow((1/r2),6)-pow((1/r2),3)); //r2 is the square of r = |\vec(r)|
    /*TODO*/
}



Vector2d PotentialLJ::F ( Vector2d r) const
{   
    double r2 = r.dot(r); //r2 is squared l2 norm of r
    return r*48/(r2)*(pow((1/r2),6)-0.5*pow((1/r2),3)); //out of Kierfeld-Skript (5.4)
    /*TODO*/
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
    /*TODO*/
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

                            Data            ( uint n, uint numBins );
    void                    save            ( const string& filenameSets,
                                              const string& filenameG,
                                              const string& filenameR ) const;
};

Data::Data( uint n, uint numBins):
    datasets( n ),  // Initializer list, because it calls constructors of the members
    rBin( numBins ),
    g( numBins, 0. ),
    r( 0 )
{
}

void Data::save ( const string& filenameSets, const string& filenameG, const string& filenameR ) const
{
    ofstream Sets(filenameSets);
    ofstream G(filenameG);
    ofstream R(filenameR);

    for(int i; i<=datasets.size(); i++) {
    Sets << datasets[i].t << ',' << datasets[i].T << ',' << datasets[i].Ekin << ',' << datasets[i].Epot << ',' << datasets[i].vS << endl;
    }
    for(int i; i<=g.size(); i++){
    G << g[i] << endl;
    }
    for(int i; i<=r.size(); i++){
    R << r[i] << endl;
    }
    
    Sets.close();
    G.close();
    R.close();
    /*TODO*/
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
        vector<Vector2d>    calc_force();
};

// Initialization of the system via constructor
MD::MD( double L, uint N, uint particlesPerRow, double T,
        Potential& potential, Thermostat& thermostat,
        uint numBins ):
    L(L),
    N(N),
    potential( potential ),
    thermostat( thermostat ),
    numBins( numBins ),
    binSize( L/numBins /*TODO*/ )
{   
    //make a grid with a particle in every square

    for(int i = 0; i<particlesPerRow; i++){
        for(int j = 0; j<particlesPerRow; j++){
            r[i+j][0] = i * L/particlesPerRow;
            r[i+j][1] = j * L/particlesPerRow;
        }
    }

    //number generator
    mt19937 rnd;
    uniform_real_distribution<> dist(-1, 1);
    Vector2d v_CoM = {0,0};

    //set the random velocities for t=0
    for(int k=0; k<N; k++){
        double random_vx = dist(rnd);
        double random_vy = dist(rnd);
        v[k]= {random_vx, random_vy};
        v_CoM = v_CoM + v[k];
    }
    
    //subtract the Center of Mass velocity
    for(int k=0; k<N; k++){
        v[k]= v[k] - 1/N * v_CoM;

    }


    //rescale
    double alpha = pow(T*2*N, 0.5);
    for(int k=0; k<N; k++){
        v[k]= v[k]*alpha;

    }


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

Data MD::measure ( const double dt, const unsigned int n )
{
    Data data(n,numBins);

    

    //verlet:
    vector<Vector2d> a_n;
    vector<Vector2d> a_n1;
    for(int steps=0; steps<n; steps++){
        data.datasets[steps].t = steps * dt; //time
        data.datasets[steps].T = calcT();
        data.datasets[steps].Ekin = calcEkin();
        data.datasets[steps].Epot = calcEpot();
        a_n = calc_force(); //calculate forces of the actual positions
        for(int i=0; i<N; i++){
            r[i] = r[i] + v[i]*dt + 0.5 * a_n[i]*dt*dt;
            if(r[i][0]>L){
                r[i][0] = r[i][0] - L; //periodic conditions
            }
            if(r[i][1]>L){
                r[i][1] = r[i][1] - L; //periodic conditions
            }
        }
        a_n1 = calc_force(); //calculate forces in new positions
        for(int i=0; i<N; i++){
            v[i] = v[i] + 0.5 *(a_n1[i]+a_n[i])*dt;
        }
    }

    data.r= r; //save r vector
    return data;
}
    
vector<Vector2d> MD::calc_force()
{
    vector<Vector2d> F;
    Vector2d r_distance; //distancevector between partciles i and j
    Vector2d nL; // nL vector
    Vector2d r_nL; // short for: r_distance + nL
    double r_nL_2; // square of r_nL
    for(int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
            r_distance = r[j] - r[i];
                for(int n_x=-1; n_x<=1 ; n_x++){ //n_x and n_y element of {-1,0,1}
                    for(int n_y=-1; n_y<=1; n_y++){
                        if(i==j && n_x==0 && n_y==0){
                            F[i] = F[i]; // if i = j and n_x = n_y = 0: dont add something
                        }
                        else{
                            nL = {n_x*L, n_y*L};
                            r_nL = r_distance + nL;
                            r_nL_2 = r_nL.dot(r_nL);
                            if(r_nL_2 <= L*L/4.0){ //if the square of the distance is smaller than the square of L/2: add something
                                F[i] = F[i] + (r_nL)/(r_nL_2)*48*(pow(1/r_nL_2,6)-0.5*pow(1/r_nL_2,3)); //calculating all forces on particle i
                            }
                            else{
                                F[i] = F[i]; // if the distance is too far: dont add something
                            }
                        }
                    }
                }
        }
        F[i] = -F[i]; // change the direction of the force

    }
    return F;
}

void MD::centerParticles()
{
    
    /*TODO*/
}

double MD::calcT() const
{   
    return 2*calcEkin()/(2*N-2);
    /*TODO*/
}

double MD::calcEkin() const
{   
    double sum;
    for (int i = 0; i < N; i++)
    {
        sum += v[i].dot(v[i]);
    }
    return 0.5*sum;
}

double MD::calcEpot() const
{
    double Epot = 0.0;
    vector<double> V;
    Vector2d r_distance; //distancevector between partciles i and j
    Vector2d nL; // nL vector
    Vector2d r_nL; // short for: r_distance + nL
    double r_nL_2; // square of r_nL
    for(int i = 0; i<N; i++){
        for (int j = 0; j<N; j++){
            r_distance = r[j] - r[i];
                for(int n_x=-1; n_x<=1 ; n_x++){ //n_x and n_y element of {-1,0,1}
                    for(int n_y=-1; n_y<=1; n_y++){
                        if(i==j && n_x==0 && n_y==0){
                            V[i] = V[i]; // if i = j and n_x = n_y = 0: dont add something
                        }
                        else{
                            nL = {n_x*L, n_y*L};
                            r_nL = r_distance + nL;
                            r_nL_2 = r_nL.dot(r_nL);
                            if(r_nL_2 <= L*L/4.0){ //if the square of the distance is smaller than the square of L/2: add something
                                V[i] = V[i] + potential.V(r_nL_2); //calculating all potentials on particle i
                            }
                            else{
                                V[i] = V[i]; // if the distance is too far: dont add something
                            }
                        }
                    }
                }
        }
        Epot = Epot + V[i];
    }

    return 0.5*Epot;

}

Vector2d MD::calcvS() const
{  
    Vector2d v_CoM; 
    for (int i = 0; i < N; i++)
    {
        v_CoM = v_CoM + v[i];
    
    }
    return v_CoM/N;
}

//Dataset MD::calcDataset() const
//{
//    /*TODO*/
//}

//Vector2d MD::calcDistanceVec( uint i, uint j ) const
//{
//    /*TODO*/
//}

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
    const uint n = 4;

    const uint partPerRow       = n;
    const uint N                = n*n;/*TODO*/
    const double L              = 2*n;/*TODO*/
    const int numBins           = 100;/*TODO*/

    // b) Equilibration test
    {
        const double T          = 1/*TODO*/;
        const double dt         = 1/*TODO*/;
        const uint steps        = 100/*TODO*/;

        MD md( L, N, partPerRow, T, LJ, noThermo, numBins );
        md.measure( dt, steps ).save( "b)set.dat", "b)g.dat", "b)r.dat" );
    }

    // c) Pair correlation function
    string TstringVec[3] = { "0.01", "1", "100" };
    for ( auto& Tstring: TstringVec )
    {
        const double T          = stod(Tstring);
        //const double dt         = /*TODO*/;
        //const uint equiSteps    = /*TODO*/;
        //const uint steps        = /*TODO*/;

        //MD md( L, N, partPerRow, T, LJ, noThermo, numBins );
        //md.equilibrate( dt, equiSteps );
        //md.measure( dt, steps ).save( "c)set" + Tstring + ".dat", "c)g" + Tstring + ".dat", "c)r" + Tstring + ".dat" );
    }

    // d) Thermostat
    {
        /*TODO*/
    }


    return 0;
}

