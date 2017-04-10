
all: $(BIN)/MultGPU $(BIN)/Generator

$(BIN)/MultGPU: Main.cpp
	g++ --std=c++11 Main.cpp Integer.o Kernel.o -o $(BIN)/MultiGPU

$(BIN)/Generator: Generator.cpp
	g++ --std=c++11 Generator.cpp -p $(BIN)/Generator
