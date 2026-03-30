CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -pedantic

SRC = src/main.cpp src/geometry.cpp src/polygon.cpp src/simplify.cpp
OUT = simplify

all: $(OUT)

$(OUT): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(OUT)

clean:
	rm -f $(OUT)