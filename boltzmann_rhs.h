// Do the physics here! 
// Inspired by S. Richers' and collaborators' implementation in EMU: https://github.com/AMReX-Astro/Emu

#pragma once
#include <array>
#include <algorithm>
#include <cmath>
#include "integrators.h"
#include "state.h"
#include "grid.h"
#include "parameters.h"
using namespace std;

// -----------------------------------------------------------------------------
// Setup Boundary Conditions - Need Ghost Cells for Periodic!
// -----------------------------------------------------------------------------

template<int nx, int nghost>
void apply_boundary_conditions(array<array<DensityMatrix, nx>, 2>& state) {
    if constexpr (boundary_condition == periodic) {
       for (int g = 0; g < nghost; g++) {
            state[0][g] = state[0][nx - 2*nghost + g];
            state[1][g] = state[1][nx - 2*nghost + g];
        }

        // Right ghosts ← Left interior
        for (int g = 0; g < nghost; g++) {
            state[0][nx - nghost + g] = state[0][nghost + g];
            state[1][nx - nghost + g] = state[1][nghost + g];
        }
    }
    else if constexpr (boundary_condition == reflecting) {
      // Left boundary
        for (int g = 0; g < nghost; g++) {
            DensityMatrix temp_plus = state[0][2*nghost - 1 - g];
            DensityMatrix temp_minus = state[1][2*nghost - 1 - g];
            
            state[0][g] = temp_minus;  // Right-moving ← Left-moving
            state[1][g] = temp_plus;   // Left-moving ← Right-moving
        }
        // Right boundary
        for (int g = 0; g < nghost; g++) {
            DensityMatrix temp_plus = state[0][nx - 2*nghost + g];
            DensityMatrix temp_minus = state[1][nx - 2*nghost + g];
            
            state[0][nx - nghost + g] = temp_minus;  // Right-moving ← Left-moving
            state[1][nx - nghost + g] = temp_plus;   // Left-moving ← Right-moving
        }   


    }
    else if constexpr (boundary_condition == outflow) {
        // Copy nearest interior to ghosts
        for (int g = 0; g < nghost; g++) {
            state[0][g] = state[0][nghost];             // Left ghost ← First interior
            state[1][g] = state[1][nghost];
            
            state[0][nx - nghost + g] = state[0][nx - nghost - 1];  // Right ghost ← Last interior
            state[1][nx - nghost + g] = state[1][nx - nghost - 1];
        }
    }
}

// -----------------------------------------------------------------------------
// Compute Hamiltonian Commutator: -i[H, ρ]
// This is the quantum mechanical term responsible for flavor oscillations
// -----------------------------------------------------------------------------
inline DensityMatrix compute_commutator(const DensityMatrix& rho, const VacuumHamiltonian& H) {
    // For 2x2 Hermitian matrices:
    // f = [ f_ee        f_ex - i*f_ex_conj ]
    //     [ f_ex + i*f_ex_conj    f_xx     ]
    //
    // H = [ H_ee    H_ex ]
    //     [ H_ex    H_xx ]

    DensityMatrix result;
    result.f_ee = +2.0 * H.H_ex * rho.f_ex_conj;
    result.f_xx = -result.f_ee;  // Opposite sign for conservation
    result.f_ex = H.H_ex * (rho.f_xx - rho.f_ee);
    result.f_ex_conj = -(H.H_ee - H.H_xx) * rho.f_ex;
    
    return result;
}

// -----------------------------------------------------------------------------
// Calculate Vacuum Oscillation Component to Change in Flavor
// 
// Compute RHS of Discrete, Two-Angle Boltzmann Equation with Oscillations
// -----------------------------------------------------------------------------
template<int nx>
array<array<DensityMatrix, nx>, 2> compute_rhs_with_oscillations(
    const array<array<DensityMatrix, nx>, 2>& state,
    const array<array<array<double, nx>, 3>, 2>& moments,
    double dx, double dt, double* max_speed) {
    
    array<array<DensityMatrix, nx>, 2> dstate_dt{};

    // Compute vacuum Hamiltonian (energy-dependent in full code)
    double omega_vac = compute_omega_vac(E_neutrino_MeV);
    VacuumHamiltonian H_vac(omega_vac);

    *max_speed = abs(c_light * mu0);

    // Loop over interior zones (exclude ghosts)
    for (int i = nghost; i < nx - nghost; i++) {
        
        // ===== DIRECTION 0: f_plus (μ > 0, moving right) =====
        
        // Upwind derivative (backward difference for rightward motion)
        DensityMatrix dfdx_plus;
        dfdx_plus.f_ee = (state[0][i].f_ee - state[0][i-1].f_ee) / dx;
        dfdx_plus.f_ex = (state[0][i].f_ex - state[0][i-1].f_ex) / dx;
        dfdx_plus.f_ex_conj = (state[0][i].f_ex_conj - state[0][i-1].f_ex_conj) / dx;
        dfdx_plus.f_xx = (state[0][i].f_xx - state[0][i-1].f_xx) / dx;
        
        // Streaming term: -μ₀ c ∂f/∂x
        DensityMatrix streaming_plus = (-c_light * mu0) * dfdx_plus;
        
        // Oscillation term: -i[H, ρ]
        DensityMatrix oscillation_plus = compute_commutator(state[0][i], H_vac);
        
        // Collision term (placeholder - set to zero for pure oscillation test)
        DensityMatrix collision_plus(0.0, 0.0, 0.0, 0.0);
        
        // Total RHS for f_plus
        dstate_dt[0][i] = streaming_plus + oscillation_plus + collision_plus;
        
        
        // ===== DIRECTION 1: f_minus (μ < 0, moving left) =====
        
        // Upwind derivative (forward difference for leftward motion)
        DensityMatrix dfdx_minus;
        dfdx_minus.f_ee = (state[1][i+1].f_ee - state[1][i].f_ee) / dx;
        dfdx_minus.f_ex = (state[1][i+1].f_ex - state[1][i].f_ex) / dx;
        dfdx_minus.f_ex_conj = (state[1][i+1].f_ex_conj - state[1][i].f_ex_conj) / dx;
        dfdx_minus.f_xx = (state[1][i+1].f_xx - state[1][i].f_xx) / dx;
        
        // Streaming term: -μ₀ c ∂f/∂x  (note: μ₀ < 0 for this direction)
        DensityMatrix streaming_minus = (-c_light * (-mu0)) * dfdx_minus;
        
        // Oscillation term: -i[H, ρ]
        DensityMatrix oscillation_minus = compute_commutator(state[1][i], H_vac);
        
        // Collision term (placeholder)
        DensityMatrix collision_minus(0.0, 0.0, 0.0, 0.0);
        
        // Total RHS for f_minus
        dstate_dt[1][i] = streaming_minus + oscillation_minus + collision_minus;
    }
    
    return dstate_dt;
}

// -----------------------------------------------------------------------------
// Wrapper to apply integrator based on compile-time choice
// -----------------------------------------------------------------------------
template<typename RHS>
void apply_integrator_choice(DensityMatrix& state_elem, double dt, RHS rhs) {
    if constexpr (integrator == rk2) {
        Integrators::rk2(state_elem, dt, rhs);
    }
    else if constexpr (integrator == rk4) { 
        Integrators::rk4(state_elem, dt, rhs);
    }
    else if constexpr (integrator == forward_euler) {
        Integrators::forward_euler(state_elem, dt, rhs);
    }
}
