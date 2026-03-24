#pragma once
#include <vector>
#include <cassert>
#include <array>
#include "parameters.h"

using State = std::array<std::array<double, nx>, 2>;
// state[0][i] = f_plus
// state[1][i] = f_minus

inline State operator+(const State& a, const State& b)
{
    State out{};
    for (int v = 0; v < 2; ++v) {
        for (int i = 0; i < nx; ++i) {
            out[v][i] = a[v][i] + b[v][i];
        }
    }
    return out;
}

inline State operator-(const State& a, const State& b)
{
    State out{};
    for (int v = 0; v < 2; ++v) {
        for (int i = 0; i < nx; ++i) {
            out[v][i] = a[v][i] - b[v][i];
        }
    }
    return out;
}

inline State operator*(double c, const State& a)
{
    State out{};
    for (int v = 0; v < 2; ++v) {
        for (int i = 0; i < nx; ++i) {
            out[v][i] = c * a[v][i];
        }
    }
    return out;
}

inline State operator*(const State& a, double c)
{
    return c * a;
}

template<int N>
State set_initial_conditions(double fplusL, double fplusR,
                             double fminusL, double fminusR)
{
    State state{};

    const int i_mid = N / 2;

    for (int i = 0; i < N; ++i) {
        if (i < i_mid) {
            state[0][i] = fplusL;
            state[1][i] = fminusL;
        } else {
            state[0][i] = fplusR;
            state[1][i] = fminusR;
        }
    }

    return state;
}