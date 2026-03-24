all:
	clang++ -std=c++17 boltzmann_solver.cpp -o boltzmann_solver.x
	./boltzmann_solver.x
