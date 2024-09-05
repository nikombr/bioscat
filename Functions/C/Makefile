# Makefile 

TARGET = gaussian_process

SRC = src/

SOURCES	= $(SRC)cholesky.cu $(SRC)cudaMalloc2d.cu $(SRC)matrix_vector_multiplication.cu $(SRC)normal_random_vector.cu $(SRC)squared_exponential.cu
OBJECTS	= $(SRC)cholesky.o $(SRC)cudaMalloc2d.o $(SRC)matrix_vector_multiplication.o $(SRC)normal_random_vector.o $(SRC)squared_exponential.o

# options and settings for the GCC compilers

CC		= gcc
CCC		= g++
CXX		= nvcc
DEFS	= 
OPT		= -g -O3
IPO		= 
ISA		= 
CHIP	= 
ARCH	= 
PARA	= 
XOPTS 	= -Xptxas=-v -arch=sm_90 -lineinfo -diag-suppress 2464 -diag-suppress 68  -diag-suppress 177 # Suppressing error regarding string conversion
CFLAGS	= $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) -lstdc++
CXXFLAGS= --compiler-options "$(OPT) $(PARA)" $(XOPTS)
LDFLAGS = -lm $(INCLUDES) $(SOFLAGS) $(XLIBS)
LIBS	= 
CUDA_PATH ?= /appl/cuda/12.2.0
INCLUDES = -I$(CUDA_PATH)/include -I$(CUDA_PATH)/samples/Common 
SOFLAGS = -L$(CUDA_PATH)/lib64 
XLIBS =  -lcudart

gaussian_process.so: src/gaussian_process.c $(OBJECTS) src/setup_covariance_matrix.o
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

src/setup_covariance_matrix.o: src/setup_covariance_matrix.cu src/squared_exponential.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

%.o: %.cu
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $< -lcudart

## %.o: %.c
## 	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	@/bin/rm -f core *.o *~
	@/bin/rm -f core src/*.o
	

realclean: clean
	@/bin/rm -f 
	@/bin/rm -f core gaussian_process.so
