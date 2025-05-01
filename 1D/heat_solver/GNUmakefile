include ../../Make.include

# Project Name
PROJ = heat_solver

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
run: $(TARGET)
	./$(TARGET) spectral sine

# Clean the executable and output data
clean:
	rm -rf $(TARGET) data

.PHONY: all run clean
