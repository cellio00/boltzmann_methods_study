// -----------------------------------------------------------------------------
// Define a moment-closure technique to make this problem solveable lol
// ----------------------------------------------------------------------------- 
#define CLOSURE_H
#pragma once
#include <array>
#include <cmath>
using namespace std;


// moments[0][i] = J
// moments[1][i] = H
// moments[2][i] = K

template<int N>
std::array<std::array<double, N>, 3>
compute_moments(const std::array<std::array<double, N>, 2>& state)
{
    std::array<std::array<double, N>, 3> moments{};

    for (int i = 0; i < N; i++) {
        const double fplus  = state[0][i];
        const double fminus = state[1][i];

        const double J = 0.5 * (fplus + fminus);
        const double H = 0.5 * mu0 * (fplus - fminus);
        const double K = mu0 * mu0 * J;  // exact 2-angle closure

        moments[0][i] = J;
        moments[1][i] = H;
        moments[2][i] = K;
    }

    return moments;
}