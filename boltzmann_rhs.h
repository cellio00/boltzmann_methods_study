// Do the physics here! 
#pragma once
#include <array>
#include <algorithm>
#include <cmath>
#include "integrators.h"
#include "state.h"
#include "grid.h"
using namespace std;

// -----------------------------------------------------------------------------
// Setup Boundary Conditions - Need Ghost Cells for Periodic!
// -----------------------------------------------------------------------------
template<int N, int NG>
void apply_boundary_conditions(std::array<std::array<double, N>, 2>& state)
{
    static_assert(NG >= 1, "Need at least one ghost zone for this implementation.");

    if constexpr (boundary_condition == periodic) {
        // Left ghost(s) from right interior
        for (int g = 0; g < NG; g++) {
            state[0][g] = state[0][N - 2*NG + g];
            state[1][g] = state[1][N - 2*NG + g];
        }

        // Right ghost(s) from left interior
        for (int g = 0; g < NG; g++) {
            state[0][N - NG + g] = state[0][NG + g];
            state[1][N - NG + g] = state[1][NG + g];
        }
    }
    else if constexpr (boundary_condition == reflecting) {
        // At a reflecting wall, f+ <-> f- swap at the boundary
        for (int g = 0; g < NG; g++) {
            state[0][g]         = state[1][2*NG - 1 - g];
            state[1][g]         = state[0][2*NG - 1 - g];

            state[0][N - NG + g] = state[1][N - NG - 1 - g];
            state[1][N - NG + g] = state[0][N - NG - 1 - g];
        }
    }
    else if constexpr (boundary_condition == outflow) {
        for (int g = 0; g < NG; g++) {
            state[0][g]          = state[0][NG];
            state[1][g]          = state[1][NG];
            state[0][N - NG + g] = state[0][N - NG - 1];
            state[1][N - NG + g] = state[1][N - NG - 1];
        }
    }
}

template<typename RHS>
void apply_integrator(State& state, double dt, RHS rhs)
{
    if constexpr (integrator == rk2) {
        Integrators::rk2(state, dt, rhs);
    }
    else if constexpr (integrator == rk4) {
        Integrators::rk4(state, dt, rhs);
    }
    else if constexpr (integrator == forward_euler) {
        Integrators::forward_euler(state, dt, rhs);
    }
}



// -----------------------------------------------------------------------------
// Compute RHS of Discrete, Two-Angle Boltzmann Equation
// -----------------------------------------------------------------------------
template<int N>
std::array<std::array<double, N>, 2>
compute_rhs(const std::array<std::array<double, N>, 2>& state,
            const std::array<std::array<double, N>, 3>& moments,
            double dx,
            double* max_signal_speed)
{
    std::array<std::array<double, N>, 2> dstate_dt{};

    *max_signal_speed = std::abs(c_light * mu0);

    // only update valid cells, not ghosts
    for (int i = nghost; i < N - nghost; i++) {
        const double fplus  = state[0][i];
        const double fminus = state[1][i];
        const double J      = moments[0][i];

        // up-wind streaming terms; approximate advection 
        // f+ travels right, so backward difference
        const double dfdx_plus =
            (state[0][i] - state[0][i - 1]) / dx;

        // f- travels left, so forward difference
        const double dfdx_minus =
            (state[1][i + 1] - state[1][i]) / dx;

        // collision/source pieces
        const double collision_plus =
            kappa * (f_eq - fplus) + sigma_s * (J - fplus);

        const double collision_minus =
            kappa * (f_eq - fminus) + sigma_s * (J - fminus);

        dstate_dt[0][i] =
            -c_light * mu0 * dfdx_plus + collision_plus;

        dstate_dt[1][i] =
            +c_light * mu0 * dfdx_minus + collision_minus;
    }

    // keep ghost-zone RHS zero! No sources there.
    for (int i = 0; i < nghost; i++) {
        dstate_dt[0][i] = 0.0;
        dstate_dt[1][i] = 0.0;
        dstate_dt[0][N - 1 - i] = 0.0;
        dstate_dt[1][N - 1 - i] = 0.0;
    }

    return dstate_dt;
}