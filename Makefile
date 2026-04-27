CXX ?= c++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic
CPPFLAGS ?= -Iinclude
BUILD_DIR ?= build

.PHONY: all run test clean

CORE_SRCS := src/solver.cpp src/phase.cpp src/optimizer.cpp src/gf2_dense.cpp src/density_kernel.cpp src/clifford_tableau.cpp src/laph.cpp
CORE_OBJS := $(CORE_SRCS:src/%.cpp=$(BUILD_DIR)/%.o)

all: $(BUILD_DIR)/laph_demo $(BUILD_DIR)/laph_tests

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: src/%.cpp include/laph/*.hpp include/laph/math/*.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/liblaph_core.a: $(CORE_OBJS) | $(BUILD_DIR)
	ar rcs $@ $(CORE_OBJS)

$(BUILD_DIR)/laph_demo: examples/demo.cpp $(BUILD_DIR)/liblaph_core.a include/laph/*.hpp include/laph/math/*.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD_DIR)/liblaph_core.a -o $@

$(BUILD_DIR)/laph_tests: tests/unit_tests.cpp $(BUILD_DIR)/liblaph_core.a include/laph/*.hpp include/laph/math/*.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(BUILD_DIR)/liblaph_core.a -o $@

run: $(BUILD_DIR)/laph_demo
	./$(BUILD_DIR)/laph_demo

test: $(BUILD_DIR)/laph_tests
	./$(BUILD_DIR)/laph_tests

clean:
	rm -rf $(BUILD_DIR)
