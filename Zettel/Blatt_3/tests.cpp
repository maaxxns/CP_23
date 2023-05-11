#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <numeric>

using namespace std;
using namespace Eigen;


double l2_norm(Vector2d const& u) {
    double s_product = u.dot(u);
    return sqrt(s_product);
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

double PotentialLJ::V ( double r2 ) const
{   
    return 4*(pow((1/r2),12)-pow((1/r2),6));
    /*TODO*/
}



Vector2d PotentialLJ::F ( Vector2d r ) const
{   
    double l2_r = l2_norm(r); //l2_r is the l2 norm of vector r 
    return r*48/(l2_r*l2_r)*(pow((1/l2_r),12)-0.5*pow((1/l2_r),6)); //out of Kierfeld-Skript (5.4)
    /*TODO*/
}

int main(){
    PotentialLJ testObj;
    Vector2d u = {3,6};
    vector<Vector2d> test(10);


    //cout << l2_norm(u) << endl;
    //cout << testObj.F(u) << endl;
    for(int i =0; i<10; i++){
        test[i] = {i,2*i};

    }
    u[1] = 2;
    cout << u;

    return 0;
}