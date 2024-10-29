# Makefile 
#SOURCES      	= src/BioScat.cpp src/forward.cpp
#OBJECTS      	= obj/BioScat.o   obj/forward.o
#SOURCES_SEG     = src/Segment.cu
#OBJECTS_SEG     = obj/Segment.o
#SOURCES_MATRIX 	= src/ComplexMatrix.cu src/RealMatrix.cu
#OBJECTS_MATRIX 	= obj/ComplexMatrix.o  obj/RealMatrix.o
#SOURCES_GP 		= src/GP/cholesky.cu src/GP/cudaMalloc2d.cu src/GP/covariance_matrix.cu  src/GP/realisation.cu src/GP/gaussian_process_inner.cpp src/GP/GaussianProcess.cpp
#OBJECTS_GP 		= obj/GP/cholesky.o  obj/GP/cudaMalloc2d.o  obj/GP/covariance_matrix.o   obj/GP/realisation.o  obj/GP/gaussian_process_inner.o   obj/GP/GaussianProcess.o
#SOURCES_2D 		= src/2D/getNanostructure.cpp src/2D/Segment.cu
#OBJECTS_2D 		= obj/2D/getNanostructure.o   obj/2D/Segment.o

# Files and dependencies

SOURCE_REALMATRIX = src/RealMatrix.cu
OBJECT_REALMATRIX = obj/RealMatrix.o
DEPEND_REALMATRIX = 

SOURCE_COMPLEXMATRIX = src/ComplexMatrix.cu
OBJECT_COMPLEXMATRIX = obj/ComplexMatrix.o
DEPEND_COMPLEXMATRIX = 

SOURCE_GET_NANOSTRUCTURE = src/getNanostructure.cpp
OBJECT_GET_NANOSTRUCTURE = obj/getNanostructure.o
DEPEND_GET_NANOSTRUCTURE = $(OBJECT_REALMATRIX)

SOURCE_SEGMENT_COMPUTE_FIELD_MATRICES = src/segment/computeFieldMatrices.cu
OBJECT_SEGMENT_COMPUTE_FIELD_MATRICES = obj/segment/computeFieldMatrices.o 
DEPEND_SEGMENT_COMPUTE_FIELD_MATRICES = $(OBJECT_REALMATRIX)

SOURCE_SEGMENT_COMPUTE_FIELD_VECTORS = src/segment/computeFieldVectors.cu
OBJECT_SEGMENT_COMPUTE_FIELD_VECTORS = obj/segment/computeFieldVectors.o
DEPEND_SEGMENT_COMPUTE_FIELD_VECTORS = $(OBJECT_REALMATRIX)

SOURCE_SEGMENT_COMPUTE_FIELDS_FOR_LINEAR_SYSTEM = src/segment/computeFieldsForLinearSystem.cu
OBJECT_SEGMENT_COMPUTE_FIELDS_FOR_LINEAR_SYSTEM = obj/segment/computeFieldsForLinearSystem.o
DEPEND_SEGMENT_COMPUTE_FIELDS_FOR_LINEAR_SYSTEM = $(OBJECT_REALMATRIX) $(OBJECT_SEGMENT_COMPUTE_FIELD_MATRICES) $(OBJECT_SEGMENT_COMPUTE_FIELD_VECTORS)

SOURCE_SEGMENT_SEGMENT = src/segment/Segment.cpp
OBJECT_SEGMENT_SEGMENT = obj/segment/Segment.o
DEPEND_SEGMENT_SEGMENT = $(OBJECT_REALMATRIX)

SOURCE_SEGMENT_SETUP = src/segment/setup.cu
OBJECT_SEGMENT_SETUP = obj/segment/setup.o
DEPEND_SEGMENT_SETUP = $(OBJECT_REALMATRIX)

SOURCE_SEGMENT_SETUP_RIGHT_HAND_SIDE = src/segment/setupRightHandSide.cu
OBJECT_SEGMENT_SETUP_RIGHT_HAND_SIDE = obj/segment/setupRightHandSide.o
DEPEND_SEGMENT_SETUP_RIGHT_HAND_SIDE = $(OBJECT_REALMATRIX)

SOURCE_SEGMENT_SETUP_SYSTEM_MATRIX = src/segment/setupSystemMatrix.cu
OBJECT_SEGMENT_SETUP_SYSTEM_MATRIX = obj/segment/setupSystemMatrix.o
DEPEND_SEGMENT_SETUP_SYSTEM_MATRIX = $(OBJECT_REALMATRIX)

SOURCE_SEGMENT_SOLVE_LINEAR_SYSTEM = src/segment/solveLinearSystem.cu
OBJECT_SEGMENT_SOLVE_LINEAR_SYSTEM = obj/segment/solveLinearSystem.o
DEPEND_SEGMENT_SOLVE_LINEAR_SYSTEM = $(OBJECT_REALMATRIX)

OBJECT_SEGMENT = $(OBJECT_SEGMENT_COMPUTE_FIELD_MATRICES) $(OBJECT_SEGMENT_COMPUTE_FIELD_VECTORS) $(OBJECT_SEGMENT_COMPUTE_FIELDS_FOR_LINEAR_SYSTEM) $(OBJECT_SEGMENT_SEGMENT) $(OBJECT_SEGMENT_SETUP) $(OBJECT_SEGMENT_SETUP_RIGHT_HAND_SIDE) $(OBJECT_SEGMENT_SETUP_SYSTEM_MATRIX) $(OBJECT_SEGMENT_SOLVE_LINEAR_SYSTEM)

SOURCE_BIOSCAT = src/BioScat.cpp
OBJECT_BIOSCAT = obj/BioScat.o
DEPEND_BIOSCAT = $(OBJECT_SEGMENT)

SOURCE_FORWARD = src/forward.cpp
OBJECT_FORWARD = obj/forward.o
DEPEND_FORWARD = $(OBJECT_COMPLEXMATRIX) $(OBJECT_BIOSCAT)

OBJECTS = $(OBJECT_REALMATRIX) $(OBJECT_COMPLEXMATRIX) $(OBJECT_SEGMENT) $(OBJECT_GET_NANOSTRUCTURE) $(OBJECT_BIOSCAT) $(OBJECT_FORWARD)

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
PARA	= -fopenmp -llapacke -llapack -lblas -lopenblas -lcusolver -lcublas -lcurand
XOPTS 	= -Xptxas=-v -arch=sm_90 -lineinfo -diag-suppress 2464 -diag-suppress 68  -diag-suppress 177 # Suppressing error regarding string conversion
CFLAGS	= $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) -lstdc++
CXXFLAGS= --compiler-options "$(OPT) $(PARA)" $(XOPTS)
LDFLAGS = -lm $(INCLUDES) $(SOFLAGS) $(XLIBS)
LIBS	= 
CUDA_PATH ?= /appl/cuda/12.2.0
INCLUDES = -I$(CUDA_PATH)/include -I$(CUDA_PATH)/samples/Common  -D_MAKE_2D
SOFLAGS = -L$(CUDA_PATH)/lib64 
XLIBS =  -lcudart

forward.so: src/executeForward.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

$(OBJECT_REALMATRIX): $(SOURCE_REALMATRIX) $(DEPEND_REALMATRIX) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_COMPLEXMATRIX): $(SOURCE_COMPLEXMATRIX) $(DEPEND_COMPLEXMATRIX) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_GET_NANOSTRUCTURE): $(SOURCE_GET_NANOSTRUCTURE) $(DEPEND_GET_NANOSTRUCTURE) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_SEGMENT_COMPUTE_FIELD_MATRICES): $(SOURCE_SEGMENT_COMPUTE_FIELD_MATRICES) $(DEPEND_SEGMENT_COMPUTE_FIELD_MATRICES) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_SEGMENT_COMPUTE_FIELD_VECTORS): $(SOURCE_SEGMENT_COMPUTE_FIELD_VECTORS) $(DEPEND_SEGMENT_COMPUTE_FIELD_VECTORS) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_SEGMENT_COMPUTE_FIELDS_FOR_LINEAR_SYSTEM): $(SOURCE_SEGMENT_COMPUTE_FIELDS_FOR_LINEAR_SYSTEM) $(DEPEND_SEGMENT_COMPUTE_FIELDS_FOR_LINEAR_SYSTEM) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_SEGMENT_SEGMENT): $(SOURCE_SEGMENT_SEGMENT) $(DEPEND_SEGMENT_SEGMENT) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_SEGMENT_SETUP): $(SOURCE_SEGMENT_SETUP) $(DEPEND_SEGMENT_SETUP) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_SEGMENT_SETUP_RIGHT_HAND_SIDE): $(SOURCE_SEGMENT_SETUP_RIGHT_HAND_SIDE) $(DEPEND_SEGMENT_SETUP_RIGHT_HAND_SIDE) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_SEGMENT_SETUP_SYSTEM_MATRIX): $(SOURCE_SEGMENT_SETUP_SYSTEM_MATRIX) $(DEPEND_SEGMENT_SETUP_SYSTEM_MATRIX) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_SEGMENT_SOLVE_LINEAR_SYSTEM): $(SOURCE_SEGMENT_SOLVE_LINEAR_SYSTEM) $(DEPEND_SEGMENT_SOLVE_LINEAR_SYSTEM) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_BIOSCAT): $(SOURCE_BIOSCAT) $(DEPEND_BIOSCAT) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

$(OBJECT_FORWARD): $(SOURCE_FORWARD) $(DEPEND_FORWARD) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

#obj/segment/%.o: src/segment/%.cu $(DEPEND_SEGMENT)
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

#obj/segment/%.o: src/segment/%.cpp $(DEPEND_SEGMENT)
#	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

##obj/2D/%.o: src/2D/%.cu $(OBJECTS_MATRIX) 
##	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

##obj/2D/%.o: src/2D/%.cpp $(OBJECTS_MATRIX) 
##	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

##obj/3D/%.o: src/3D/%.cu
##	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

##obj/3D/%.o: src/3D/%.cpp
##	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

##obj/GP/%.o: src/GP/%.cu
##	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

##obj/GP/%.o: src/GP/%.cpp
##	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

## %.o: %.c
## 	$(CC) $(CFLAGS) -o $@ -c $<


clean:
	@/bin/rm -f core *.o *~
	@/bin/rm -f core obj/*.o
	@/bin/rm -f core obj/GP/*.o
	@/bin/rm -f core obj/2D/*.o
	@/bin/rm -f core obj/3D/*.o
	@/bin/rm -f core obj/segment/*.o
	

realclean: clean
	@/bin/rm -f 
	@/bin/rm -f core forward.so
	@/bin/rm -f core backward.so
