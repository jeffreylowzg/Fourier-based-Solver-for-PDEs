include ../../Make.include

# Project Name
PROJ = convergence_test

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
	./$(TARGET) be vary_n

# Clean the executable and output data
clean:
	rm -rf $(TARGET) data

.PHONY: all run clean
