
#--------------------------------------General------------------------------------------
FC      := mpif90
#SRCDIR  := .
BUILDDIR := build
TARGET  := Selfk_LU_parallel_3D.x


#--------------------------------------Sources and header files------------------------------------------
SRCEXT  := f90
#SOURCES := $(shell ls -1 $(SRCDIR)/*.$(SRCEXT))
SOURCES := dispersion.f90 calc_susc.f90 sigma.f90 read.f90 write.f90 vardef.f90 lambda_correction.f90 Selfk_LU_parallel.f90
#OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJECTS := $(patsubst %,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

#--------------------------------------Settings for fRG flow------------------------------------------

# FFLAGS += -D SOME_DEF_VARIABLE # Set additional defines

#--------------------------------------Compiler settings------------------------------------------

FFLAGS += # General compiler flags
RUNFLAGS := -O4 # Compiler flags for quick compile and run
#DBFLAGS := -C -traceback -g -fpe0 # Compiler flags for debugging
DBFLAGS := -Og -g -fcheck=all
PROFFLAGS := -O2 -g # Compiler flags for profiling
OPTIMIZEFLAGS := -O4 # Compiler flags for optimal speed
########################################################
### be very careful here: in case of using MKL instead of LAPACK, the FFTW-libraries have
### to be loaded first, in order to avoid conflicts! (the MKL library also have FFT functions 
### with the same name than FFTW, but not all of FFTW)
########################################################
#LIB := -llapack # Specify Libraries
LIB := -L/lrz/sys/libraries/fftw/3.3.3/sse/lib -lfftw3f -lfftw3 -lfftw3l -L/lrz/sys/intel/compiler/composer_xe_2015.2.164/mkl/lib/intel64 -lmkl_rt # Specify Libraries
INC := # Additional include paths


#--------------------------------------Targets ------------------------------------------
$(TARGET): $(OBJECTS)
	@echo " Linking..."
#	@mkdir -p bin
	@echo " $(FC) $^ -o $(TARGET) $(LIB)"; $(FC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: %.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(FC) $(FFLAGS) $(INC) -c -o $@ $<"; $(FC) $(FFLAGS) $(INC) -c -o $@ $<

run:    FFLAGS += $(RUNFLAGS)
run:    $(TARGET)

debug:  FFLAGS += $(DBFLAGS)
debug:  $(TARGET)

prof:   FFLAGS += $(PROFFLAGS)
prof:   $(TARGET)

optimize: FFLAGS += $(OPTIMIZEFLAGS)
optimize: $(SOURCES) $(HEADERS)
	@mkdir -p bin
	$(FC) $(FFLAGS) $(INC) $(SOURCES) -o $(TARGET) $(LIB)
	#rm *.o

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)
	@echo " rm mod"; rm *.mod

.PHONY: clean
