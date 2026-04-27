CXX ?= c++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -pedantic
CPPFLAGS ?= -Iinclude
BUILD_DIR ?= build

.PHONY: all run test clean

all: $(BUILD_DIR)/laph_demo $(BUILD_DIR)/laph_tests

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/laph_demo: examples/demo.cpp include/laph/laph.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(BUILD_DIR)/laph_tests: tests/unit_tests.cpp include/laph/laph.hpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

run: $(BUILD_DIR)/laph_demo
	./$(BUILD_DIR)/laph_demo

test: $(BUILD_DIR)/laph_tests
	./$(BUILD_DIR)/laph_tests

clean:
	rm -rf $(BUILD_DIR)
