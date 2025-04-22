# Compiler
CXX = g++

# Project Name
PROJ = heat_solver

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -O3

# Libraries (link FFTW)
LIBS = -lfftw3
LIBS += -I/home/jrohozen/fftw3/include -L/home/jrohozen/fftw3/lib

# Executable name
TARGET = $(PROJ).exe

# Source files
SOURCES = main.cpp heat_solver.cpp

# Default target
all: $(TARGET)

# Linking and compiling
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LIBS)

# Run the simulation
# Optionally run with the "spectral" argument to use the FFTW spectral method:
#   ./heat_solver.exe spectral
run: $(TARGET)
	./$(TARGET)

# Clean the executable and output data
clean:
	rm -f $(TARGET)
	rm -rf data

.PHONY: all run clean
