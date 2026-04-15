CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -march=native
LDFLAGS = 

TARGET = boltzmann_solver.x

SOURCES = boltzmann_solver.cpp
HEADERS = state.h \
          parameters.h \
          boltzmann_rhs.h \
          closure.h \
          integrators.h \
          grid.h


all: $(TARGET)

# Build the executable
$(TARGET): $(SOURCES) $(HEADERS)
	@echo "Compiling $(TARGET)..."
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LDFLAGS)
	@echo "Done"

# Clean build artifacts
clean:
	@echo "Cleaning..."
	rm -f $(TARGET) output.dat *.o

debug: CXXFLAGS = -std=c++17 -O0 -g -Wall -Wextra
debug: $(TARGET)


.PHONY: all clean run plot debug help
