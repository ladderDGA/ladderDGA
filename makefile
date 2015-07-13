
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
DBFLAGS := -C -traceback -g -fpe0 # Compiler flags for debugging
PROFFLAGS := -O2 -g # Compiler flags for profiling
OPTIMIZEFLAGS := -O4 # Compiler flags for optimal speed
LIB := -llapack # Specify Libraries
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
