#pragma once
#include <vector>
#include <fstream>
#include <array>

// cell-centered coordinates
inline double x_coordinate(int i, double dx) {
    return (i + 0.5) * dx;
}

template<int N>
void print(std::ofstream& output,
           int it,
           double t,
           double dx,
           const std::array<std::array<double, N>, 2>& state,
           const std::array<std::array<double, N>, 3>& moments)
{
    for (int i = 0; i < N; i++) {
        output << it << " "
               << t << " "
               << i << " "
               << x_coordinate(i, dx) << " "
               << state[0][i] << " "
               << state[1][i] << " "
               << moments[0][i] << " "
               << moments[1][i] << " "
               << moments[2][i] << "\n";
    }
    output << "\n";
}