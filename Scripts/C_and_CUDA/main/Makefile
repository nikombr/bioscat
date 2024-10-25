# Makefile 
SOURCES      	= src/BioScat.cpp src/forward.cpp
OBJECTS      	= obj/BioScat.o   obj/forward.o   
SOURCES_MATRIX 	= src/ComplexMatrix.cpp src/RealMatrix.cpp
OBJECTS_MATRIX 	= obj/ComplexMatrix.o   obj/RealMatrix.o
SOURCES_GP 		= src/GP/cholesky.cu src/GP/cudaMalloc2d.cu src/GP/covariance_matrix.cu  src/GP/realisation.cu src/GP/gaussian_process_inner.cpp src/GP/GaussianProcess.cpp
OBJECTS_GP 		= obj/GP/cholesky.o  obj/GP/cudaMalloc2d.o  obj/GP/covariance_matrix.o   obj/GP/realisation.o  obj/GP/gaussian_process_inner.o   obj/GP/GaussianProcess.o
SOURCES_2D 		= src/2D/getNanostructure.cpp src/2D/getSegments.cu
OBJECTS_2D 		= obj/2D/getNanostructure.o   obj/2D/getSegments.o

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
PARA	= -fopenmp -llapack -lblas -lopenblas -lcusolver -lcublas -lcurand
XOPTS 	= -Xptxas=-v -arch=sm_90 -lineinfo -diag-suppress 2464 -diag-suppress 68  -diag-suppress 177 # Suppressing error regarding string conversion
CFLAGS	= $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) -lstdc++
CXXFLAGS= --compiler-options "$(OPT) $(PARA)" $(XOPTS)
LDFLAGS = -lm $(INCLUDES) $(SOFLAGS) $(XLIBS)
LIBS	= 
CUDA_PATH ?= /appl/cuda/12.2.0
INCLUDES = -I$(CUDA_PATH)/include -I$(CUDA_PATH)/samples/Common 
SOFLAGS = -L$(CUDA_PATH)/lib64 
XLIBS =  -lcudart

forward.so: src/executeForward.c $(OBJECTS_MATRIX) $(OBJECTS_2D) $(OBJECTS_GP) $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

obj/%.o: src/%.cu $(OBJECTS_MATRIX) $(OBJECTS_2D)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/%.o: src/%.cpp $(OBJECTS_MATRIX)  $(OBJECTS_2D)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/2D/%.o: src/2D/%.cu $(OBJECTS_MATRIX) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/2D/%.o: src/2D/%.cpp $(OBJECTS_MATRIX) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/3D/%.o: src/3D/%.cu
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/3D/%.o: src/3D/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/GP/%.o: src/GP/%.cu
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/GP/%.o: src/GP/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

## %.o: %.c
## 	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	@/bin/rm -f core *.o *~
	@/bin/rm -f core obj/*.o
	@/bin/rm -f core obj/GP/*.o
	@/bin/rm -f core obj/2D/*.o
	@/bin/rm -f core obj/3D/*.o
	

realclean: clean
	@/bin/rm -f 
	@/bin/rm -f core forward.so
	@/bin/rm -f core backward.so
