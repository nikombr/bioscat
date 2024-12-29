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

SOURCES_UTILS = src/utils/cudaMalloc2d.cu src/utils/ComplexMatrix.cu src/utils/RealMatrix.cu
OBJECTS_UTILS = obj/utils/cudaMalloc2d.o  obj/utils/ComplexMatrix.o  obj/utils/RealMatrix.o

SOURCES_GAUSSIAN_PROCESS = src/GaussianProcess/cholesky.cu src/GaussianProcess/compute_inverse.cu src/GaussianProcess/compute_prior.cu  src/GaussianProcess/covariance_matrix.cu  src/GaussianProcess/GaussianProcess.cpp  src/GaussianProcess/realisation.cu 
OBJECTS_GAUSSIAN_PROCESS = obj/GaussianProcess/cholesky.o  obj/GaussianProcess/compute_inverse.o  obj/GaussianProcess/compute_prior.o   obj/GaussianProcess/covariance_matrix.o   obj/GaussianProcess/GaussianProcess.o    obj/GaussianProcess/realisation.o 
DEPENDS_GAUSSIAN_PROCESS = $(OBJECTS_UTILS)

SOURCES_SEGMENT_DEPEND = src/Segment/depend/computeFieldMatricesKernel.cu src/Segment/depend/computeFieldsFromMatricesKernel.cu src/Segment/depend/kernels.cu
OBJECTS_SEGMENT_DEPEND = obj/Segment/depend/computeFieldMatricesKernel.o  obj/Segment/depend/computeFieldsFromMatricesKernel.o  obj/Segment/depend/kernels.o

SOURCES_SEGMENT_SETUP = src/Segment/setup/computeExteriorPointsAndNormalVectors.cu src/Segment/setup/computeInteriorPoints.cu src/Segment/setup/computeTestPoints.cu
OBJECTS_SEGMENT_SETUP = obj/Segment/setup/computeExteriorPointsAndNormalVectors.o  obj/Segment/setup/computeInteriorPoints.o  obj/Segment/setup/computeTestPoints.o
DEPENDS_SEGMENT_SETUP = $(OBJECTS_SEGMENT_DEPEND) $(OBJECTS_UTILS)

SOURCES_SEGMENT = src/Segment/allocate.cpp src/Segment/computeFarFieldPattern.cu src/Segment/computeIncidentFields.cu src/Segment/computeInteriorFields.cu src/Segment/computeInteriorFieldMatrices.cu src/Segment/computeScatteredFields.cu src/Segment/computeScatteredFieldMatrices.cu src/Segment/free.cpp src/Segment/Segment.cpp src/Segment/setup.cpp src/Segment/setupRightHandSide.cu src/Segment/setupSystemMatrix.cu src/Segment/solveLinearSystem.cu
OBJECTS_SEGMENT = obj/Segment/allocate.o   obj/Segment/computeFarFieldPattern.o  obj/Segment/computeIncidentFields.o  obj/Segment/computeInteriorFields.o  obj/Segment/computeInteriorFieldMatrices.o  obj/Segment/computeScatteredFields.o  obj/Segment/computeScatteredFieldMatrices.o  obj/Segment/free.o   obj/Segment/Segment.o   obj/Segment/setup.o   obj/Segment/setupRightHandSide.o  obj/Segment/setupSystemMatrix.o  obj/Segment/solveLinearSystem.o
DEPENDS_SEGMENT = $(OBJECTS_SEGMENT_SETUP) $(OBJECTS_SEGMENT_DEPEND) $(OBJECTS_UTILS)

SOURCES_SEGMENT_ALL_DEPEND = src/Segment/alldepend/computeFieldsForLinearSystem.cpp
OBJECTS_SEGMENT_ALL_DEPEND = obj/Segment/alldepend/computeFieldsForLinearSystem.o
DEPENDS_SEGMENT_ALL_DEPEND = $(OBJECTS_SEGMENT_SETUP) $(OBJECTS_SEGMENT_DEPEND) $(OBJECTS_SEGMENT) $(OBJECTS_UTILS)

SOURCES_BIOSCAT_DEPEND = src/BioScat/depend/combinePolarisation.cu
OBJECTS_BIOSCAT_DEPEND = obj/BioScat/depend/combinePolarisation.o
DEPENDS_BIOSCAT_DEPEND = $(OBJECTS_GAUSSIAN_PROCESS) $(OBJECTS_UTILS) $(OBJECTS_SEGMENT_DEPEND) $(OBJECTS_SEGMENT_SETUP) $(OBJECTS_SEGMENT) $(OBJECTS_SEGMENT_ALL_DEPEND)

SOURCES_BIOSCAT = src/BioScat/allocateSegments.cpp src/BioScat/BioScat.cpp src/BioScat/computeFields.cpp src/BioScat/computeSubFields.cu src/BioScat/dumpFarFields.cpp src/BioScat/dumpFields.cpp src/BioScat/fowardSolver.cpp src/BioScat/free.cpp src/BioScat/prepareForward.cpp src/BioScat/reset.cpp src/BioScat/setupObservationPoints.cpp src/BioScat/setupSegments.cpp src/BioScat/computeReflectance.cpp src/BioScat/getNanostructure.cpp src/BioScat/computeFarFieldPattern.cpp
OBJECTS_BIOSCAT = obj/BioScat/allocateSegments.o   obj/BioScat/BioScat.o   obj/BioScat/computeFields.o   obj/BioScat/computeSubFields.o  obj/BioScat/dumpFarFields.o   obj/BioScat/dumpFields.o   obj/BioScat/fowardSolver.o   obj/BioScat/free.o   obj/BioScat/prepareForward.o   obj/BioScat/reset.o   obj/BioScat/setupObservationPoints.o   obj/BioScat/setupSegments.o   obj/BioScat/computeReflectance.o   obj/BioScat/getNanostructure.o   obj/BioScat/computeFarFieldPattern.o
DEPENDS_BIOSCAT = $(OBJECTS_GAUSSIAN_PROCESS) $(OBJECTS_UTILS) $(OBJECTS_SEGMENT_DEPEND) $(OBJECTS_SEGMENT_SETUP) $(OBJECTS_SEGMENT) $(OBJECTS_SEGMENT_ALL_DEPEND) $(OBJECTS_BIOSCAT_DEPEND)

SOURCES_REFLECTANCE = src/reflectance/computeReflectanceMatrix.cu
OBJECTS_REFLECTANCE = obj/reflectance/computeReflectanceMatrix.o
DEPENDS_REFLECTANCE = $(OBJECTS_GAUSSIAN_PROCESS) $(OBJECTS_UTILS) $(OBJECTS_SEGMENT_DEPEND) $(OBJECTS_SEGMENT_SETUP) $(OBJECTS_SEGMENT) $(OBJECTS_SEGMENT_ALL_DEPEND) $(OBJECTS_BIOSCAT) $(OBJECTS_BIOSCAT_DEPEND)

SOURCES_EXECPPCU = src/execppcu/computeFarFieldPattern.cpp src/execppcu/computeGaussianProcess.cpp src/execppcu/forward.cpp src/execppcu/generateArtificialData.cpp src/execppcu/inverse.cpp src/execppcu/profileGaussianProcess.cpp src/execppcu/profileForward.cpp
OBJECTS_EXECPPCU = obj/execppcu/computeFarFieldPattern.o   obj/execppcu/computeGaussianProcess.o   obj/execppcu/forward.o   obj/execppcu/generateArtificialData.o   obj/execppcu/inverse.o   obj/execppcu/profileGaussianProcess.o   obj/execppcu/profileForward.o
DEPENDS_EXECPPCU = $(OBJECTS_GAUSSIAN_PROCESS) $(OBJECTS_UTILS) $(OBJECTS_SEGMENT_DEPEND) $(OBJECTS_SEGMENT_SETUP) $(OBJECTS_SEGMENT) $(OBJECTS_SEGMENT_ALL_DEPEND) $(OBJECTS_REFLECTANCE) $(OBJECTS_BIOSCAT_DEPEND)

OBJECTS = $(OBJECTS_GAUSSIAN_PROCESS) $(OBJECTS_UTILS) $(OBJECTS_SEGMENT_DEPEND) $(OBJECTS_SEGMENT_SETUP) $(OBJECTS_SEGMENT) $(OBJECTS_SEGMENT_ALL_DEPEND) $(OBJECTS_EXECPPCU) $(OBJECTS_REFLECTANCE) $(OBJECTS_BIOSCAT) $(OBJECTS_BIOSCAT_DEPEND)

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
PARA	= -fopenmp -llapacke -llapack -lblas -lopenblas -lcusolver -lcublas -lcurand -Ilib -Ilib/utils -Ilib/reflectance -Ilib/Segment 
XOPTS 	= -Xptxas=-v -arch=sm_90 -lineinfo -diag-suppress 2464 -diag-suppress 68  -diag-suppress 177 # Suppressing error regarding string conversion
CFLAGS	= $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) -lstdc++
CXXFLAGS= --compiler-options "$(OPT) $(PARA)" $(XOPTS)
LDFLAGS = -lm $(INCLUDES) $(SOFLAGS) $(XLIBS)
LIBS	= 
CUDA_PATH ?= /appl/cuda/12.2.0
INCLUDES = -I$(CUDA_PATH)/include -I$(CUDA_PATH)/samples/Common  -D_MAKE_2D
SOFLAGS = -L$(CUDA_PATH)/lib64 
XLIBS =  -lcudart

#all: $(OBJECTS_GAUSSIAN_PROCESS) $(OBJECTS_UTILS) $(OBJECTS_SEGMENT)
all: so/forward.so so/computeGaussianProcess.so so/generateArtificialData.so so/inverse.so so/computeFarFieldPattern.so so/profileGaussianProcess.so so/profileForward.so

so/forward.so: src/exec/executeForward.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

so/profileForward.so: src/exec/executeProfileForward.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

so/computeFarFieldPattern.so: src/exec/executeComputeFarFieldPattern.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

so/inverse.so: src/exec/executeInverse.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

so/computeGaussianProcess.so: src/exec/executeComputeGaussianProcess.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

so/profileGaussianProcess.so: src/exec/executeProfileGaussianProcess.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

so/generateArtificialData.so: src/exec/executeGenerateArtificialData.c $(OBJECTS)
	$(CC) $(CFLAGS) $(INCLUDES)  -fPIC -shared -o $@ $^ $(SOFLAGS) $(LIBS)  -lcudart

obj/GaussianProcess/%.o: src/GaussianProcess/%.cu $(DEPENDS_GAUSSIAN_PROCESS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/GaussianProcess/%.o: src/GaussianProcess/%.cpp $(DEPENDS_GAUSSIAN_PROCESS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/utils/%.o: src/utils/%.cu 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/utils/%.o: src/utils/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/Segment/depend/%.o: src/Segment/depend/%.cu 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/Segment/depend/%.o: src/Segment/depend/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/Segment/setup/%.o: src/Segment/setup/%.cu $(DEPENDS_SEGMENT_SETUP)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/Segment/setup/%.o: src/Segment/setup/%.cpp $(DEPENDS_SEGMENT_SETUP)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/Segment/%.o: src/Segment/%.cu $(DEPENDS_SEGMENT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/Segment/%.o: src/Segment/%.cpp $(DEPENDS_SEGMENT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/Segment/alldepend/%.o: src/Segment/alldepend/%.cu $(DEPENDS_SEGMENT_ALL_DEPEND)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/Segment/alldepend/%.o: src/Segment/alldepend/%.cpp $(DEPENDS_SEGMENT_ALL_DEPEND)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/BioScat/%.o: src/BioScat/%.cu $(DEPENDS_BIOSCAT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/BioScat/%.o: src/BioScat/%.cpp $(DEPENDS_BIOSCAT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/BioScat/depend/%.o: src/BioScat/depend/%.cu $(DEPENDS_BIOSCAT_DEPEND)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/BioScat/depend/%.o: src/BioScat/depend/%.cpp $(DEPENDS_BIOSCAT_DEPEND)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/reflectance/%.o: src/reflectance/%.cu $(DEPENDS_REFLECTANCE)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/reflectance/%.o: src/reflectance/%.cpp $(DEPENDS_REFLECTANCE)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/execppcu/%.o: src/execppcu/%.cu $(DEPENDS_EXECPPCU)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart

obj/execppcu/%.o: src/execppcu/%.cpp $(DEPENDS_EXECPPCU)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -Xcompiler -fPIC -shared -o $@  $^ -lcudart
	

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
	@/bin/rm -f core obj/GaussianProcess/*.o
	@/bin/rm -f core obj/Segment/*.o
	@/bin/rm -f core obj/Segment/alldepend/*.o
	@/bin/rm -f core obj/Segment/setup/*.o
	@/bin/rm -f core obj/Segment/depend/*.o
	@/bin/rm -f core obj/BioScat/*.o
	@/bin/rm -f core obj/BioScat/depend/*.o
	@/bin/rm -f core obj/utils/*.o
	@/bin/rm -f core obj/execppcu/*.o
	

realclean: clean
	@/bin/rm -f 
	@/bin/rm -f core so/*.o
