//implement different integration methods here 
#pragma once
#include <cmath>
#include "state.h"

namespace Integrators {
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
