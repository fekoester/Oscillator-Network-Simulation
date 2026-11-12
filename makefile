CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
CPPFLAGS += -Iinclude
HEADERS = $(wildcard include/oscnet/*.hpp)

all: oscnet

oscnet: src/main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

oscnet_tests: tests/test_oscnet.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

test: oscnet_tests
	./oscnet_tests

clean:
	rm -f oscnet oscnet_tests

.PHONY: all test clean
