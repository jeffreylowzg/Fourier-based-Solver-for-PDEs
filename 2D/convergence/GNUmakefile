# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -O3

# Libraries (link FFTW)
LIBS = -lfftw3
LIBS += -I/home/jrohozen/fftw3/include -L/home/jrohozen/fftw3/lib

# Executable names
TARGET = convergence.exe

# Source files for the 2D solver
SOURCES = main.cpp heat_solver2d.cpp

all: $(TARGET)

# Build the 2D solver executable
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LIBS)

# Run the convergence test
run: $(TARGET) 
	./$(TARGET) be vary_n vary_dt

# Clean the executables and output data
clean:
	rm -f $(TARGET)
	rm -rf data

.PHONY: all run clean
