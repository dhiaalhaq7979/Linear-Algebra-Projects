CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Iinclude

build/linear_algebra: src/main.cpp src/Matrix.cpp include/Matrix.h
	mkdir -p build
	$(CXX) $(CXXFLAGS) src/main.cpp src/Matrix.cpp -o $@

run: build/linear_algebra
	./build/linear_algebra

clean:
	rm -rf build
