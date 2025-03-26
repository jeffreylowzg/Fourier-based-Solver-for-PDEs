# Compiler
CXX = g++

# Name
PROJ = heat_solver

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra -O3

# Executable name
TARGET = $(PROJ).exe

# Source files
SOURCES = heat_solver.cpp

# Default target
all: $(TARGET)

# Linking and compiling
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES)

# Run the simulation
run: $(TARGET)
	./$(TARGET)

# Clean the executable and output data
clean:
	rm -f $(TARGET) data/

.PHONY: all run clean
