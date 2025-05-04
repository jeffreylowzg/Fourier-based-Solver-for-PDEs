include ../../Make.include

# Executable names
TARGET = heat_solver2d.exe

# Source files for the 2D solver
SOURCES = main.cpp heat_solver2d.cpp

# Default target: build both executables
all: $(TARGET)

# Build the 2D solver executable
$(TARGET): $(SOURCES)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCES) $(LIBS)

# Run the 2D solver simulation.
# Use the initial condition argument: sine or gaussian.
run: $(TARGET) 
	./$(TARGET) spectral_be sine none

# Clean the executables and output data
clean:
	rm -rf $(TARGET) data

.PHONY: all run clean
