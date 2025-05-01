include ../../Make.include

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
	rm -rf $(TARGET) data

.PHONY: all run clean
