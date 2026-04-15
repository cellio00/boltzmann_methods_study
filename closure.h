#pragma once
#include <array>
#include <cmath>
#include "state.h"
#include "parameters.h"
using namespace std;

// Compute angular moments from distribution functions
// For two-flavor system, we compute flavor-diagonal moments
//
// J_α = ∫ f_α dΩ  (number density for flavor α)
// H_α = ∫ μ f_α dΩ  (number flux for flavor α) 
// K_α = ∫ μ² f_α dΩ  (pressure/stress for flavor α)
//
// In discrete ordinates with two angles (±μ₀):
// J_α ≈ 0.5 * (f_α^+ + f_α^-)  
// H_α ≈ 0.5 * μ₀ * (f_α^+ - f_α^-)
// K_α ≈ 0.5 * μ₀² * (f_α^+ + f_α^-)

template<int N>
std::array<std::array<std::array<double, N>, 3>, 2>
compute_moments(const std::array<std::array<DensityMatrix, N>, 2>& state)
{
    // moments[flavor][moment_type][zone]
    // flavor: 0 = electron, 1 = x-flavor (mu/tau)
    // moment_type: 0 = J (density), 1 = H (flux), 2 = K (pressure)
    std::array<std::array<std::array<double, N>, 3>, 2> moments{};

    for (int i = 0; i < N; i++) {
        // Extract density matrix elements for both directions
        const double f_ee_plus = state[0][i].f_ee;
        const double f_xx_plus = state[0][i].f_xx;
        
        const double f_ee_minus = state[1][i].f_ee;
        const double f_xx_minus = state[1][i].f_xx;
        
        // Electron neutrino moments
        const double J_e = 0.5 * (f_ee_plus + f_ee_minus);
        const double H_e = 0.5 * mu0 * (f_ee_plus - f_ee_minus);
        const double K_e = 0.5 * mu0 * mu0 * (f_ee_plus + f_ee_minus);
        
        // X-flavor neutrino moments  
        const double J_x = 0.5 * (f_xx_plus + f_xx_minus);
        const double H_x = 0.5 * mu0 * (f_xx_plus - f_xx_minus);
        const double K_x = 0.5 * mu0 * mu0 * (f_xx_plus + f_xx_minus);
        
        // Store electron neutrino moments
        moments[0][0][i] = J_e;
        moments[0][1][i] = H_e;
        moments[0][2][i] = K_e;
        
        // Store x-flavor neutrino moments
        moments[1][0][i] = J_x;
        moments[1][1][i] = H_x;
        moments[1][2][i] = K_x;
    }

    return moments;
}

// Compute total (flavor-summed) moments
// Useful for checking conservation
template<int N>
std::array<std::array<double, N>, 3>
compute_total_moments(const std::array<std::array<DensityMatrix, N>, 2>& state)
{
    std::array<std::array<double, N>, 3> total_moments{};
    
    auto flavor_moments = compute_moments<N>(state);
    
    for (int i = 0; i < N; i++) {
        // Sum over flavors
        total_moments[0][i] = flavor_moments[0][0][i] + flavor_moments[1][0][i];  // J_total
        total_moments[1][i] = flavor_moments[0][1][i] + flavor_moments[1][1][i];  // H_total
        total_moments[2][i] = flavor_moments[0][2][i] + flavor_moments[1][2][i];  // K_total
    }
    
    return total_moments;
}

// Compute flavor-mixing diagnostics
// Useful for visualizing oscillations
template<int N>
std::array<std::array<double, N>, 2>
compute_flavor_coherence(const std::array<std::array<DensityMatrix, N>, 2>& state)
{
    // Returns: [0] = |ρ_ex| (magnitude of coherence)
    //          [1] = arg(ρ_ex) (phase of coherence)
    std::array<std::array<double, N>, 2> coherence{};
    
    for (int i = 0; i < N; i++) {
        // Average over directions (could also look at each separately)
        double rho_ex_real = 0.5 * (state[0][i].f_ex + state[1][i].f_ex);
        double rho_ex_imag = 0.5 * (state[0][i].f_ex_conj + state[1][i].f_ex_conj);
        
        coherence[0][i] = std::sqrt(rho_ex_real*rho_ex_real + rho_ex_imag*rho_ex_imag);
        coherence[1][i] = std::atan2(rho_ex_imag, rho_ex_real);
    }
    
    return coherence;
}
