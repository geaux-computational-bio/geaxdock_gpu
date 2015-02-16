
#COMPUTE_DEVICE := USE_GPU
COMPUTE_DEVICE := USE_CPU
#COMPUTE_DEVICE := USE_MIC

CPU_COMPILER := GCC
#CPU_COMPILER := ICC







EXE := dock

ifeq ($(COMPUTE_DEVICE), USE_GPU)
	OBJ_CPU := dock.o load.o data.o rmsd.o util.o hdf5io.o
	OBJ_GPU := run.o

	SH := sh
	CPP_HOST := g++
	CPP_DEV := nvcc
	FC := gfortran

	HOSTFLAGS :=
	DEVFLAGS := -gencode arch=compute_20,code=sm_20
#	DEVFLAGS := -gencode=arch=compute_20,code=\"sm_20,compute_20\"
#	DEVFLAGS := -gencode arch=compute_30,code=sm_30
#	DEVFLAGS := -gencode arch=compute_35,code=\"sm_35,compute_35\"
#	DEVFLAGS := -arch=sm_35
# -ccbin=g++ # --ptxas-options=-v #-keep #-rdc=true
	LINKFLAGS := -lcudart -lhdf5 -lm
else ifeq ($(COMPUTE_DEVICE), USE_CPU)
	OBJ_CPU := dock.o load.o data.o rmsd.o util_cpu.o run_cpu.o
	OBJ_GPU :=

	ifeq ($(CPU_COMPILER), GCC)
		SH := sh
		CPP_HOST := g++
		FC := gfortran

		HOSTFLAGS := -fopenmp
		LINKFLAGS := -lm
	else ifeq ($(CPU_COMPILER), ICC)
		SH := sh
		CPP_HOST := icc
		FC := ifort

		HOSTFLAGS := -openmp #-no-vec #-vec-report3
		LINKFLAGS := -lm
	endif
else ifeq ($(COMPUTE_DEVICE), USE_MIC)
	OBJ_CPU := dock.o load.o data.o rmsd.o util_cpu.o run_cpu.o
	OBJ_GPU :=

	SH := sh
	CPP_HOST := icc
	FC := ifort

	HOSTFLAGS := -openmp -mmic #-no-vec #-vec-report3
	LINKFLAGS := -lm
endif




OPTFLAGS := -O3

HEADPATH := -I./include -I./include/cuda

#LIBPATH := -L/usr/local/cuda-4.2.9/lib64
#LIBPATH := -L/usr/local/cuda-5.5/lib64
#LIBPATH := -L/usr/local/packages/cuda/4.2.9/cuda/lib64 -L/usr/local/packages/hdf5/1.8.10/gcc-4.4.6-openmpi-1.6.2/lib
#LIBPATH := -L/usr/local/packages/cuda/5.0/lib64 -L/usr/local/packages/hdf5/1.8.10/gcc-4.4.6-openmpi-1.6.2/lib
#LIBPATH := -L/usr/local/packages/cuda/5.0/lib64 -L/usr/local/packages/hdf5/1.8.9/Intel-13.0.0-openmpi-1.6.2/lib
#LIBPATH := -L/opt/intel/lib/intel64



HOSTFLAGS += -Wall $(OPTFLAGS) $(HEADPATH)
DEVFLAGS += $(OPTFLAGS) $(HEADPATH)
FFLAGS := $(OPTFLAGS)

# -D COMPUTE_DEVICE=$(COMPUTE_DEVICE)







default: $(EXE) #analysis

$(EXE): $(OBJ_CPU) $(OBJ_GPU) 
	$(CPP_HOST) $(HOSTFLAGS) $(LIBPATH) $(OBJ_CPU) $(OBJ_GPU) -o $@ $(LINKFLAGS) 
	@cp $(EXE) ../bin/

hdf5io.o: hdf5io.C
	h5c++ -c $<

%.o: %.C
	$(CPP_HOST) $(HOSTFLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.cu
	$(CPP_DEV) $(DEVFLAGS) -c $<


count: count.o load.o data.o rmsd.o
	$(CPP_HOST) -O0 -g $^ -o $@ $(LINKFLAGS)
	@cp count ../test/
	@echo "copy count to ../test/"

analysis: analysis.o util.o hdf5io.o
	$(CPP_HOST) $(LIBPATH) $^ -o $@ $(LINKFLAGS)

clean:
	@(rm -f ${EXE} ${OBJ_CPU} ${OBJ_GPU} analysis count *i *cubin* *.ptx *cudafe* *.fatbin* *.hash count count.o ../bin/dock)


