# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -O3

# Libraries (link FFTW)
LIBS = -lfftw3

# Executable names
TARGET = reaction2D.exe

# Source files for the 2D solver
SOURCES = main.cpp reaction2D.cpp

# Default target: build both executables
all: $(TARGET)

# Build the 2D solver executable
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LIBS)

# Run the 2D solver simulation.
# Use the initial condition argument: sine or gaussian.
run: $(TARGET) 
	./$(TARGET) sine 

# Clean the executables and output data
clean:
	rm -f $(TARGET)
	rm -rf data

.PHONY: all run clean
