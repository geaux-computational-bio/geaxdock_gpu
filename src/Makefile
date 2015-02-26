# MAKEFILE CONSTANTS

TARGET_GPU := 0
TARGET_CPU := 1
TARGET_MIC := 2

# GPUs per compute node
NGPU := 1










EXE := dock
OBJ_CPU := dock.o load.o data.o rmsd.o util.o hdf5io.o
OBJ_GPU := run.o
SH := sh
CPP_HOST := g++
CPP_DEV := nvcc
FC := gfortran


DMACRO := -DTARGET_GPU=$(TARGET_GPU) -DTARGET_CPU=$(TARGET_CPU) -DTARGET_MIC=$(TARGET_MIC)
GPU := GPU_NA







# HOST: mike, shelob, qb, smic las07
HOST := qb

ifeq ($(HOST), mike)
# edit include/yeah/cudautil.h, locate "CUDA_LDG_P", enable the second branch
# edit ~/.soft
# +cuda-5.0
# +hdf5-1.8.10-gcc-4.4.6-openmpi-1.6.2
	LIBPATH := -L/usr/local/packages/cuda/5.0/lib64 -L/usr/local/packages/hdf5/1.8.10/gcc-4.4.6-openmpi-1.6.2/lib
	DMARCRO += -DTARGET_DEVICE=$(TARGET_GPU)
	GPU := GF110

else ifeq ($(HOST), shelob)
# edit ~/.soft
# +cuda-5.5.22
# +hdf5-1.8.10-gcc-4.4.6-openmpi-1.6.2
	LIBPATH := -L/usr/local/packages/cuda/5.5.22/lib64 -L/usr/local/packages/hdf5/1.8.10/gcc-4.4.6/lib
	DMARCRO += -DTARGET_DEVICE=$(TARGET_GPU)
	GPU := GK110

else ifeq ($(HOST), qb)
	LIBPATH := -L/usr/local/packages/cuda/5.5.22/lib64 -L/usr/local/packages/hdf5/1.8.10/gcc-4.4.6/lib
	DMARCRO += -DTARGET_DEVICE=$(TARGET_GPU)
	GPU := GK110

else ifeq ($(HOST), smic)
# module load intel/14.0.2
# module load cuda/6.5
# module load hdf5/1.8.12/INTEL-140-MVAPICH2-2.0
# module unload gcc/4.9.0
# module unload GCC/4.9.0
	DMARCRO += -DTARGET_DEVICE=$(TARGET_MIC)

else ifeq ($(HOST), las07)
	LIBPATH := -L/opt/cuda/cuda-6.5/lib64 -L/opt/hdf5/1.8.13_gcc/lib
	DMARCRO += -DTARGET_DEVICE=$(TARGET_GPU)
	GPU := GM204
endif



# best CUDA kernel configurations
# GPU		Fermi M2090	Kepler K20Xm	Maxwell GTX980
#		GF110		GK110		GM204
# ---------------------------------------------------------------
# GD		16		56		64				GridDim.x
# BDy		8		8		8				BlockDim.y
# BDx		64		96		16				BlockDim.x
# TperB		512		768		128				Threads per Block

ifeq ($(GPU), GF110)
	DEVFLAGS := -gencode arch=compute_20,code=sm_20
	GD := 16
	BDy := 8
	BDx := 64
	TperB := 512
else ifeq ($(GPU), GK110)
	DEVFLAGS := -gencode arch=compute_35,code=sm_35 -Xptxas -dlcm=ca
	GD := 15
	BDy := 8
	BDx := 96
	TperB := 768
else ifeq ($(GPU), GM204)
	DEVFLAGS := -gencode arch=compute_52,code=compute_52 -Xptxas -dlcm=ca
	GD := 64
	BDy := 8
	BDx := 16
	TperB := 128
else
	GD := 16
	BDy := 8
	BDx := 32
	TperB := 256
endif

DMARCRO_GPU := -DGD=$(GD) -DBDy=$(BDy) -DBDx=$(BDx) -DTperB=$(TperB) -DNGPU=$(NGPU)





HOSTFLAGS :=
HEADPATH := -I./modules -I./modules/nvidia_gpucomputingsdk_4.2.9_c_common_inc
OPTFLAGS := -O3
LINKFLAGS := -lcudart -lhdf5 -lm


HOSTFLAGS += -Wall $(OPTFLAGS) $(HEADPATH) $(DMARCRO)
DEVFLAGS += $(OPTFLAGS) $(HEADPATH) $(DMARCRO) $(DMARCRO_GPU)
FFLAGS := $(OPTFLAGS)







default: $(EXE)

$(EXE): $(OBJ_CPU) $(OBJ_GPU) 
	$(CPP_HOST) $(HOSTFLAGS) $(LIBPATH) $(OBJ_CPU) $(OBJ_GPU) -o $@ $(LINKFLAGS)

hdf5io.o: hdf5io.C
	h5c++ -c $<

%.o: %.C
	$(CPP_HOST) $(HOSTFLAGS) -c $<

%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.cu
	$(CPP_DEV) $(DEVFLAGS) -c $<

install:
	@cp $(EXE) ../work/

clean:
	@(rm -f ${EXE} ${OBJ_CPU} ${OBJ_GPU} count *i *cubin* *.ptx *cudafe* *.fatbin* *.hash count count.o ../bin/dock)


