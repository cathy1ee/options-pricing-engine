CXX      ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra
REQUIRED := -std=c++17 -Iinclude
SRC      := src/black_scholes.cpp src/monte_carlo.cpp

all: test bench

build:
	mkdir -p build

build/test_pricing: tests/test_pricing.cpp $(SRC) include/pricing.hpp | build
	$(CXX) $(REQUIRED) $(CXXFLAGS) tests/test_pricing.cpp $(SRC) -o $@

build/convergence: bench/convergence.cpp $(SRC) include/pricing.hpp | build
	$(CXX) $(REQUIRED) $(CXXFLAGS) bench/convergence.cpp $(SRC) -o $@

test: build/test_pricing
	./build/test_pricing

bench: build/convergence
	./build/convergence

clean:
	rm -rf build

.PHONY: all test bench clean
