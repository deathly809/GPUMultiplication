
$(BIN)/Kernel.o: Kernel.cu 
	nvcc -c --std=c++11 --use_fast_math -m 64 --gpu-architecture=compute_30 Kernel.cu -o $(BIN)/Kernel.o
