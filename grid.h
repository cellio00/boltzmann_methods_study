#pragma once
#include <vector>
#include <fstream>
#include <array>
#include <iomanip>
#include "state.h"

// cell-centered coordinates
inline double x_coordinate(int i, double dx) {
    return (i + 0.5) * dx;
}

// Print function for outputting state and moments
template<int N>
void print(std::ofstream& output,
           int it,
           double t,
           double dx,
           const std::array<std::array<DensityMatrix, N>, 2>& state,
           const std::array<std::array<std::array<double, N>, 3>, 2>& moments)
{
    // Output format:
    // index(t) t index(x) x f_ee_+ f_ex_+ f_ex_conj_+ f_xx_+ f_ee_- f_ex_- f_ex_conj_- f_xx_- J_e H_e K_e J_x H_x K_x
    
    output << std::scientific << std::setprecision(6);
    
    for (int i = 0; i < N; i++) {
        output << it << " "
               << t << " "
               << i << " "
               << x_coordinate(i, dx) << " ";
        
        // f_plus components
        output << state[0][i].f_ee << " "
               << state[0][i].f_ex << " "
               << state[0][i].f_ex_conj << " "
               << state[0][i].f_xx << " ";
        
        // f_minus components
        output << state[1][i].f_ee << " "
               << state[1][i].f_ex << " "
               << state[1][i].f_ex_conj << " "
               << state[1][i].f_xx << " ";
        
        // Electron neutrino moments
        output << moments[0][0][i] << " "  // J_e
               << moments[0][1][i] << " "  // H_e
               << moments[0][2][i] << " "; // K_e
        
        // X-flavor neutrino moments
        output << moments[1][0][i] << " "  // J_x
               << moments[1][1][i] << " "  // H_x
               << moments[1][2][i] << "\n"; // K_x
    }
    output << "\n";  // Extra newline between timesteps for gnuplot
    output << std::flush;
}

// Header for output file
inline void print_header(std::ofstream& output) {
    output << "# Columns:\n";
    output << "# 1:it 2:t 3:ix 4:x ";
    output << "5:f_ee+ 6:f_ex+ 7:f_ex_conj+ 8:f_xx+ ";
    output << "9:f_ee- 10:f_ex- 11:f_ex_conj- 12:f_xx- ";
    output << "13:J_e 14:H_e 15:K_e ";
    output << "16:J_x 17:H_x 18:K_x\n";
}
