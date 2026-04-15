// -----------------------------------------------------------------------------
// Boltzmann Toy Solver Computational Project, Carrie Elliott, Physics 643
// -----------------------------------------------------------------------------
// Two-Flavor Neutrino Vacuum Oscillations in 1D with Two-Angle Discrete Ordinates
//
// This code solves the quantum kinetic equations for neutrino flavor evolution:
//   ∂ρ/∂t + μc ∂ρ/∂x = -i[H, ρ] + C[ρ]
//
// where ρ is the 2×2 flavor density matrix, H is the Hamiltonian, and C is collisions.
//
#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <cassert>

#include "boltzmann_rhs.h"
#include "grid.h"
#include "parameters.h"
#include "closure.h"
#include "integrators.h"

using namespace std;

int main() {
    cout << "======================================================\n";
    cout << " Two-Flavor Neutrino Oscillation Transport Solver\n";
    cout << "======================================================\n";
    cout << "Grid: " << nx_valid << " zones, dx = " << dx << "\n";
    cout << "Time: tend = " << tend << ", CFL = " << courant_factor << "\n";
    cout << "Mixing angle: θ₁₂ = " << theta12_deg << "°\n";
    cout << "Δm²₂₁ = " << delta_m2_21_eV2 << " eV²\n";
    cout << "Neutrino energy: E = " << E_neutrino_MeV << " MeV\n";
    cout << "Oscillation length: L_osc = " << compute_L_osc(E_neutrino_MeV) << " cm\n";
    cout << "Vacuum frequency: ω = " << compute_omega_vac(E_neutrino_MeV) << " cm⁻¹\n";
    cout << "======================================================\n\n";

    // State variables
    // state[direction][zone] where direction: 0 = +mu, 1 = -mu
    array<array<DensityMatrix, nx>, 2> state;

    // Moments: moments[flavor][moment][zone]
    // flavor: 0 = electron, 1 = heavy (x)
    // moment: 0 = J, 1 = H, 2 = K
    array<array<array<double, nx>, 3>, 2> moments;

    // Set initial conditions
    state = set_initial_conditions<nx>(fplusL, fplusR, fminusL, fminusR);
    
    cout << "Initial conditions set:\n";
    cout << "  Left side: f_ee = " << fplusL << " (pure electron flavor)\n";
    cout << "  Right side: f_xx = " << fplusR << " (pure heavy flavor)\n\n";

    // Apply boundary conditions before first use
    apply_boundary_conditions<nx, nghost>(state);

    // Compute initial moments
    moments = compute_moments<nx>(state);

    // Set up output
    int it = 0;
    double t = 0.0;
    ofstream output;
    output.open("output.dat");
    print_header(output);
    print<nx>(output, it, t, dx, state, moments);
    
    cout << "Starting time integration...\n";
    
    // Start time integration
    bool end = false;
    while (!end) {

        // Define timestep size based on CFL condition
        double max_signal_speed = c_light * mu0;
        double dt = courant_factor * dx / max_signal_speed;

        // Adjust final timestep to hit tend exactly
        if (t + dt > tend) {
            dt = tend - t;
            end = true;
        }

        // Define RHS function for the integrator
        auto rhs = [&](const array<array<DensityMatrix, nx>, 2>& s) {
            auto m = compute_moments<nx>(s);
            double speed;
            return compute_rhs_with_oscillations<nx>(s, m, dx, dt, &speed);
        };

        // Evolve state using chosen integrator
        if constexpr (integrator == rk2) {
            Integrators::rk2(state, dt, rhs);
        }
        else if constexpr (integrator == rk4) { 
            Integrators::rk4(state, dt, rhs);
        }
        else if constexpr (integrator == forward_euler) {
            Integrators::forward_euler(state, dt, rhs);
        }

        // Apply boundary conditions after evolution
        apply_boundary_conditions<nx, nghost>(state);

        // Compute moments for output
        moments = compute_moments<nx>(state);

        // Update time
        it++;
        t += dt;

        // Output every timestep (or could do every N steps)
        print<nx>(output, it, t, dx, state, moments);
        
        // Progress indicator
        if (it % 10 == 0) {
            cout << "  it = " << it << ", t = " << t << " / " << tend << "\r" << flush;
        }
    }

    cout << "\nSimulation complete!\n";
    cout << "  Total iterations: " << it << "\n";
    
    // Check conservation (trace of density matrix should be conserved)
    double initial_trace_sum = 0.0;
    double final_trace_sum = 0.0;
    
    auto initial_state = set_initial_conditions<nx>(fplusL, fplusR, fminusL, fminusR);
    for (int i = nghost; i < nx - nghost; i++) {
        initial_trace_sum += initial_state[0][i].trace() + initial_state[1][i].trace();
        final_trace_sum += state[0][i].trace() + state[1][i].trace();
    }

// Implementing some physical safeguards; i.e. check that the resule remains traceless
cout << "\n======================================================\n";
cout << "CONSERVATION CHECK\n";
cout << "======================================================\n";


// Compute final trace
double final_trace = 0.0;
for (int i = nghost; i < nx - nghost; i++) {
    final_trace += state[0][i].trace() + state[1][i].trace();
}

// Compute initial trace
auto initial_state_ref = set_initial_conditions<nx>(fplusL, fplusR, fminusL, fminusR);
double initial_trace = 0.0;
for (int i = nghost; i < nx - nghost; i++) {
    initial_trace += initial_state_ref[0][i].trace() + initial_state_ref[1][i].trace();
}

cout << "Initial trace: " << initial_trace << "\n";
cout << "Final trace:   " << final_trace << "\n";
cout << "Change:        " << final_trace - initial_trace << "\n";

if constexpr (boundary_condition == reflecting) {
    // Reflecting: trace should be exactly conserved
    double error = abs(final_trace - initial_trace) / initial_trace;
    cout << "Relative error: " << error << "\n";
    
    if (error < 1e-8) {
        cout << " EXCELLENT: Conservation to < 1e-8\n";
    } else if (error < 1e-6) {
        cout << " GOOD: Conservation to < 1e-6\n";
    } else {
        cout << " POOR: Conservation error too large!\n";
        cout << "   → Check for bugs in evolution equations\n";
    }
    
} else if constexpr (boundary_condition == outflow) {
    // Outflow: trace should decrease, never increase
    double ratio = final_trace / initial_trace;
    cout << "Ratio (final/initial): " << ratio << "\n";
    
    if (ratio > 1.01) {
        cout << " ERROR: Trace INCREASED by " << (ratio-1)*100 << "%\n";
        cout << "   → Particles being created (BUG!)\n";
        cout << "   → Check commutator signs\n";
    } else if (ratio > 0.99) {
        cout << "  WARNING: Trace barely changed\n";
        cout << "   → Beams haven't exited yet (run longer?)\n";
    } else {
        cout << "✓ Trace decreased to " << ratio*100 << "% of initial\n";
        cout << "  (Expected: particles exited domain)\n";
    }
    
} else if constexpr (boundary_condition == periodic) {
    cout << " WARNING: Periodic BC not appropriate for collision!\n";
    cout << "   → Use reflecting (for testing) or outflow (for physics)\n";
}

cout << "======================================================\n";
    output.close();
    return 0;
}
