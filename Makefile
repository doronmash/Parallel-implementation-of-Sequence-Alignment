build:
	mpicxx -fopenmp -c main.c -o main.o
	mpicxx -fopenmp -c fileUtil.c -o fileUtil.o
	mpicxx -fopenmp -c auxiliaryfunctions.c -o auxiliaryfunctions.o
	mpicxx -fopenmp -c cFunctions.c -o cFunctions.o
	nvcc -I./common/inc -c cudaFunctions.cu -o cudaFunctions.o
	mpicxx -fopenmp -o program main.o fileUtil.o auxiliaryfunctions.o cFunctions.o 		cudaFunctions.o /usr/local/cuda/lib64/libcudart_static.a -ldl -lrt
clean:
	rm -f *.o ./program
	
run:
	mpiexec -np 2 ./program
