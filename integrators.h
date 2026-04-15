#pragma once
#include <cmath>
#include "state.h"

namespace Integrators {
    
    // Forward Euler: u^{n+1} = u^n + dt * f(u^n)
    template <typename RHS>
    void forward_euler(DensityMatrix& u, double dt, RHS rhs) {
        DensityMatrix k1 = rhs(u);
        u = u + dt * k1;
    }

    // RK2 (Heun's method / Improved Euler)
    // u^{n+1} = u^n + dt * f(u^n + dt/2 * f(u^n))
    template <typename RHS>
    void rk2(DensityMatrix& u, double dt, RHS rhs) {
        DensityMatrix k1 = rhs(u);
        DensityMatrix u_mid = u + 0.5 * dt * k1;
        DensityMatrix k2 = rhs(u_mid);
        u = u + dt * k2;
    }

    // RK4 (Classic 4th order Runge-Kutta)
    template <typename RHS>
    void rk4(DensityMatrix& u, double dt, RHS rhs) {
        DensityMatrix k1 = rhs(u);
        DensityMatrix k2 = rhs(u + 0.5 * dt * k1);
        DensityMatrix k3 = rhs(u + 0.5 * dt * k2);
        DensityMatrix k4 = rhs(u + dt * k3);

        u = u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    
    // State-level integrators (for full state array)
    template <typename RHS>
    void forward_euler(State& u, double dt, RHS rhs) {
        State k1 = rhs(u);
        u = u + dt * k1;
    }
    
    template <typename RHS>
    void rk2(State& u, double dt, RHS rhs) {
        State k1 = rhs(u);
        State u_mid = u + 0.5 * dt * k1;
        State k2 = rhs(u_mid);
        u = u + dt * k2;
    }
    
    template <typename RHS>
    void rk4(State& u, double dt, RHS rhs) {
        State k1 = rhs(u);
        State k2 = rhs(u + 0.5 * dt * k1);
        State k3 = rhs(u + 0.5 * dt * k2);
        State k4 = rhs(u + dt * k3);

        u = u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
}
