CXX      ?= g++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -Wpedantic -Iinclude
DEPFLAGS := -MMD -MP
BUILD    := build
OBJ      := $(BUILD)/obj
BIN      := $(BUILD)/bin

LIB_SRC  := src/black_scholes.cpp src/implied_vol.cpp src/binomial.cpp src/normal_model.cpp
LIB_OBJ  := $(LIB_SRC:%.cpp=$(OBJ)/%.o)
ALL_OBJ  := $(LIB_OBJ) $(OBJ)/src/main.o $(OBJ)/tests/tests.o

.PHONY: all test clean run
all: $(BIN)/bachelier

# -MMD -MP emits a .d file per object listing the headers it included, so a
# header edit rebuilds every translation unit that sees it. Without this, a
# struct layout change silently leaves stale objects linked against each other.
$(OBJ)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

$(BIN)/bachelier: $(LIB_OBJ) $(OBJ)/src/main.o
	@mkdir -p $(BIN)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN)/bachelier_tests: $(LIB_OBJ) $(OBJ)/tests/tests.o
	@mkdir -p $(BIN)
	$(CXX) $(CXXFLAGS) $^ -o $@

test: $(BIN)/bachelier_tests
	@./$(BIN)/bachelier_tests

run: $(BIN)/bachelier
	@./$(BIN)/bachelier price -s 100 -k 100 -t 1 -r 0.05 -v 0.2

clean:
	@rm -rf $(BUILD)

-include $(ALL_OBJ:.o=.d)
