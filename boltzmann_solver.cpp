// -----------------------------------------------------------------------------
// Boltzmann Toy Solver Computational Project, Carrie Elliott, Physics 643
// -----------------------------------------------------------------------------
// It would be longer to say all of the assumptions made by this code,
// and it in no way aims to make meaningful physical assertions about neutrino transport.
// However, the goal is to hopefully demonstrate *why* solving the neutrino transport problem,
// even without oscillation, let alone with, is going to be a tough undertaking! 
#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <cassert>

#include "boltzmann_rhs.h"
#include "grid.h"
#include "parameters.h"
#include "state.h"
#include "closure.h"

using namespace std;

int main() {
  // evolved variables
  // state[0][i] = f_plus
  // state[1][i] = f_minus
  array<array<double, nx>, 2> state;

  // moments
  // moments[0][i] = J , number density
  // moments[1][i] = H , number flux 
  // moments[2][i] = K , number pressure? pressure tensor 
  array<array<double, nx>, 3> moments;

  // RHS
  array<array<double, nx>, 2> dstate_dt;

  // characteristic speed for CFL
  double max_signal_speed;

  // set initial conditions for the +/- distribution functions
  state = set_initial_conditions<nx>(fplusL, fplusR, fminusL, fminusR);

  // apply BCs before first use
  apply_boundary_conditions<nx, nghost>(state);

  // compute initial moments
  moments = compute_moments<nx>(state);

  // set up output
  int it = 0;
  double t = 0.0;
  ofstream output;
  output.open("output.dat");
  output << "# it t ix x f_plus f_minus J H K" << endl;
  print<nx>(output, it, t, dx, state, moments);

  // start time integration
  bool end = false;
  while (!end) {

    // define timestep size
    double max_signal_speed = c_light * mu0;
    double dt = courant_factor * dx / max_signal_speed;

    // Define time bounds
    if (t + dt > tend) {
        dt = tend - t;
        end = true;
    }

    // Define RHS of the 1D Boltzmann equation 
    auto rhs = [&](const State& s) {
        auto m = compute_moments<nx>(s);
        double speed;
        return compute_rhs<nx>(s, m, dx, &speed);
    };

    // evolve in time via function in boltzmann_rhs
    apply_integrator(state, dt, rhs);

    // apply BCs via function in boltzmann_rhs
    apply_boundary_conditions<nx, nghost>(state);

    // output
    auto moments = compute_moments<nx>(state);

    it++;
    t += dt;

    print<nx>(output, it, t, dx, state, moments);
}

  output.close();
  return 0;
}