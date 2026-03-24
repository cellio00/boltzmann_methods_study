// -----------------------------------------------------------------------------
// Compile-time / run parameters for toy 1D Boltzmann solver
// -----------------------------------------------------------------------------
#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>


// Grid
constexpr int nx = 100;          // total number of stored zones
constexpr int nghost = 1;        // ghost zones on each side
constexpr int nx_valid = 100;    // physical zones
constexpr double dx = 1.0 / nx;  // requested choice

// Time integration
constexpr double courant_factor = 0.4;
constexpr double tend = 0.4;
enum Integrator_Type {
    forward_euler = 0,
    rk2 = 1,
    rk4 = 2
};

// Transport parameters
constexpr double mu0 = 1.0;      // discrete propagation cosine(theta)
constexpr double c_light = 1.0;  // massless relativistic neutrino limit, set c=1
constexpr double kappa = 0.0;    // absorption / emissivity
constexpr double sigma_s = 0.5;  // isotropizing scattering opacity
constexpr double f_eq = 0.0;     // equilibrium distribution target

// Initial conditions
constexpr double fplusL  = 1.2;
constexpr double fplusR  = 0.0;
constexpr double fminusL = 0.0;
constexpr double fminusR = 1.34;

// Boundary conditions
enum BoundaryCondition {
    periodic = 0,
    reflecting = 1,
    outflow = 2
};

// Manually set boundary condition type and integration type at compile time
constexpr BoundaryCondition boundary_condition = periodic;
constexpr Integrator_Type integrator = forward_euler;