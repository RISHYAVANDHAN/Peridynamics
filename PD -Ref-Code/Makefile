CXX = g++
CXXFLAGS = -fopenmp -O3
SOURCE = main.cpp

OUTPUT = main

all: $(OUTPUT)

$(OUTPUT): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(SOURCE)

clean:
	rm -f $(OUTPUT)
