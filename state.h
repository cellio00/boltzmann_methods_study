#pragma once
#include <vector>
#include <cassert>
#include <array>
#include <cmath>
#include "parameters.h"

// Density matrix for two-flavor neutrino system
// f = [ f_ee    f_ex         ]
//     [ f_ex*   f_xx         ]
// where f_ex is complex, stored as real and imaginary parts

struct DensityMatrix {
    double f_ee;       // electron neutrino diagonal element
    double f_ex;       // off-diagonal real part (e-x mixing)
    double f_ex_conj;  // off-diagonal imaginary part 
    double f_xx;       // heavy flavor (mu/tau) diagonal element
    
    // Constructor for easy initialization
    DensityMatrix(double ee = 0.0, double ex_real = 0.0, double ex_imag = 0.0, double xx = 0.0) 
        : f_ee(ee), f_ex(ex_real), f_ex_conj(ex_imag), f_xx(xx) {}
    
    // Hermiticity check (for debugging)
    bool is_hermitian(double tol = 1e-10) const {
        return (f_ee >= 0 && f_xx >= 0);  // diagonal elements should be non-negative
    }
    
    // Trace (should be conserved in vacuum oscillations)
    double trace() const {
        return f_ee + f_xx;
    }
};

// Operators for DensityMatrix
inline DensityMatrix operator+(const DensityMatrix& a, const DensityMatrix& b) {
    return DensityMatrix(
        a.f_ee + b.f_ee,
        a.f_ex + b.f_ex,
        a.f_ex_conj + b.f_ex_conj,
        a.f_xx + b.f_xx
    );
}

inline DensityMatrix operator-(const DensityMatrix& a, const DensityMatrix& b) {
    return DensityMatrix(
        a.f_ee - b.f_ee,
        a.f_ex - b.f_ex,
        a.f_ex_conj - b.f_ex_conj,
        a.f_xx - b.f_xx
    );
}

inline DensityMatrix operator*(double c, const DensityMatrix& a) {
    return DensityMatrix(
        c * a.f_ee,
        c * a.f_ex,
        c * a.f_ex_conj,
        c * a.f_xx
    );
}

inline DensityMatrix operator*(const DensityMatrix& a, double c) {
    return c * a;
}

// State type: array of density matrices for + and - directions at each zone
// state[direction][zone] where direction: 0 = +, 1 = -
using State = std::array<std::array<DensityMatrix, nx>, 2>;

// Operators for State
inline State operator+(const State& a, const State& b) {
    State out{};
    for (int v = 0; v < 2; ++v) {
        for (int i = 0; i < nx; ++i) {
            out[v][i] = a[v][i] + b[v][i];
        }
    }
    return out;
}

inline State operator-(const State& a, const State& b) {
    State out{};
    for (int v = 0; v < 2; ++v) {
        for (int i = 0; i < nx; ++i) {
            out[v][i] = a[v][i] - b[v][i];
        }
    }
    return out;
}

inline State operator*(double c, const State& a) {
    State out{};
    for (int v = 0; v < 2; ++v) {
        for (int i = 0; i < nx; ++i) {
            out[v][i] = c * a[v][i];
        }
    }
    return out;
}

inline State operator*(const State& a, double c) {
    return c * a;
}

// Initialize state with two-beam collision setup
// Left side: right-moving electron neutrinos (f_plus = f_ee)
// Right side: left-moving x-flavor neutrinos (f_minus = f_xx)
template<int N>
State set_initial_conditions(double fplusL, double fplusR,
                             double fminusL, double fminusR) {
    State state{};
    const int i_mid = N / 2;

    for (int i = 0; i < N; ++i) {
        if (i < i_mid) {
            state[0][i] = DensityMatrix(fplusL, 0.0, 0.0, 0.0);  // f_ee = fplusL
            state[1][i] = DensityMatrix(fminusL, 0.0, 0.0, 0.0);  // all zeros if fminusL=0
        } else {
            state[0][i] = DensityMatrix(fplusR, 0.0, 0.0, 0.0);   // all zeros if fplusR=0
            // f_minus (←): x-flavor neutrinos
            state[1][i] = DensityMatrix(0.0, 0.0, 0.0, fminusR);  // f_xx = fminusR
        }
    }

    return state;
}
