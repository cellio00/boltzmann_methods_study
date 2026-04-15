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
constexpr int nx_valid = nx - 2*nghost;  // = 98
constexpr double dx = 1.0 / nx_valid;

// Time integration
constexpr double courant_factor = 0.1;
constexpr double tend = 1.0;

enum Integrator_Type {
    forward_euler = 0,
    rk2 = 1,
    rk4 = 2
};

// Transport parameters
constexpr double mu0 = 1.0;         // discrete propagation cosine(theta)
constexpr double c_light = 1.0;    // massless relativistic neutrino limit, set c=1
constexpr double kappa = 0.25;      // absorption / emissivity
constexpr double sigma_s = 0.3;    // isotropic scattering opacity (set to 0 for pure oscillations)

// Neutrino oscillation parameters (two-flavor approximation)
// Using solar neutrino parameters as example
constexpr double theta12_deg = 33.82;  // vacuum mixing angle in degrees [PDG 2020: 33.82°]
constexpr double theta12_rad = theta12_deg * M_PI / 180.0;  // in radians

// Mass-squared differences
constexpr double delta_m2_21_eV2 = 7.53e-5;  // Δm²₂₁ in eV² [PDG 2020: 7.53×10⁻⁵ eV²]

// Neutrino energy (this should ideally be per energy bin in a real code)
constexpr double E_neutrino_MeV = 10.0;  // neutrino energy in MeV

// Physical constants
constexpr double eV_to_erg = 1.60218e-12;       // eV to erg conversion
constexpr double MeV_to_eV = 1.0e6;             // MeV to eV
constexpr double c_cgs = 2.99792458e10;         // speed of light in cm/s
constexpr double hbar_cgs = 1.05457266e-27;     // reduced Planck constant in erg s
constexpr double hbar_eV_s = 6.582119569e-16;   // hbar in eV s

// Vacuum oscillation frequency
// wvac = Δm2 / (2E) in natural units
// For proper units: need to convert carefully
inline double compute_omega_vac(double E_MeV) {
    double E_eV = E_MeV * MeV_to_eV;
    double omega_eV = delta_m2_21_eV2 / (2.0 * E_eV);  // in eV

    double hbar_c_eV_cm = hbar_eV_s * c_cgs;  // ℏc in eV·cm
    return omega_eV / hbar_c_eV_cm;
}

// Oscillation length in cm
inline double compute_L_osc(double E_MeV) {
    // L_osc = 2π/w = 4piE/Δm2
    // Simplified formula: L_osc[km] = 2.48 × E[MeV] / Δm²[eV²]
    return 2.48e5 * E_MeV / delta_m2_21_eV2;  // in cm (×1e5 for km to cm)
}

// Vacuum Hamiltonian in flavor basis (2x2 matrix elements)
// H = ω_vac × [ -cos(2θ)    sin(2θ)  ]
//             [  sin(2θ)    cos(2θ)  ]
struct VacuumHamiltonian {
    double H_ee;
    double H_ex;  // real part (Hermitian matrix)
    double H_xx;
    
    VacuumHamiltonian(double omega) {
        double cos2theta = std::cos(2.0 * theta12_rad);
        double sin2theta = std::sin(2.0 * theta12_rad);
        
        H_ee = - cos2theta;
        H_ex =  sin2theta;
        H_xx =  cos2theta;
    }
};

// Weak interaction parameters (for future use with matter effects)
constexpr double GF_GeV2 = 1.1663787e-5;  // Fermi constant in GeV⁻²

// Initial conditions - Two-beam collision model
// Electron neutrino beam from left moving RIGHT
// X-flavor neutrino beam from right moving LEFT
// They oscillate as they pass through each other
constexpr double fplusL  = 1.0;   // Left beam: electron neutrinos moving RIGHT 
constexpr double fplusR  = 0.0;   // (no right-moving on right side)
constexpr double fminusL = 0.0;   // (no left-moving on left side)
constexpr double fminusR = 1.0;   // Right beam: x-flavor neutrinos moving LEFT 

// Boundary conditions
enum BoundaryCondition {
    periodic = 0,
    reflecting = 1,
    outflow = 2
};

// Manually set boundary condition type and integration type at compile time
constexpr BoundaryCondition boundary_condition = reflecting;
constexpr Integrator_Type integrator = forward_euler;  // Use RK2 for better accuracy with oscillations
