#include <iostream>
#include <random>
#include <cmath>
#include <vector>

using namespace std;

// Define the constants
const double epsilon = 1.0;
const double sigma = 1.0;
const double rc = 0.5;
const double h = 0.01;
const double k_B = 1.0;

// Define the Lennard-Jones potential function
double LJ_potential(double r) {
    double r6 = pow(r, 6);
    double r12 = pow(r6, 2);
    return 4 * epsilon * (r12 / pow(sigma, 12) - r6 / pow(sigma, 6));
}

// Define the function to compute the forces and energies
void compute_forces(vector<double>& x, vector<double>& y, vector<double>& fx, vector<double>& fy, double& potential_energy, double& kinetic_energy) {
    int N = x.size();
    double L = 2.0 * sigma;
    potential_energy = 0.0;
    fx.assign(N, 0.0);
    fy.assign(N, 0.0);

    for (int i = 0; i < N; i++) {
        for (int j = i+1; j < N; j++) {
            double dx = x[j] - x[i];
            double dy = y[j] - y[i];
            // Apply periodic boundary conditions
            if (dx > L/2.0) {
                dx -= L;
            }
            else if (dx < -L/2.0) {
                dx += L;
            }
            if (dy > L/2.0) {
                dy -= L;
            }
            else if (dy < -L/2.0) {
                dy += L;
            }
            double r = sqrt(dx*dx + dy*dy);
            if (r < rc) {
                double f = LJ_potential(r) / r;
                fx[i] += f * dx;
                fy[i] += f * dy;
                fx[j] -= f * dx;
                fy[j] -= f * dy;
                potential_energy += LJ_potential(r);
            }
        }
    }

    // Compute the kinetic energy
    kinetic_energy = 0.0;
    for (int i = 0; i < N; i++) {
        kinetic_energy += 0.5 * (fx[i]*fx[i] + fy[i]*fy[i]);
    }
}

// Define the function to rescale velocities to set the temperature
void rescale_velocities(vector<double>& vx, vector<double>& vy, double target_temperature) {
    int N = vx.size();
    double current_temperature = 0.0;
    for (int i = 0; i < N; i++) {
        current_temperature += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i]);
    }
    current_temperature /= k_B * (N - 1);
    double scale_factor = sqrt(target_temperature / current_temperature);
    for (int i = 0; i < N; i++) {
        vx[i] *= scale_factor;
        vy[i] *= scale_factor;
    }
}

// Define the main function
int main() {
    // Set the simulation parameters
    int n = 20;
    int N = n*n;
    double L = 2.0 * sigma * n;
    double T = 1.0;
    double dt = h;

    // Initialize the particle positions and velocities
    vector<double> x(N), y(N), vx(N);