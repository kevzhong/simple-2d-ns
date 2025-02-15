FC = gfortran

SRCDIR = src
OBJDIR = obj
OUTDIR = outputdir

TARGET = a.out

FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch -funroll-loops -std=f2008 -Wall -Wextra -Wpedantic -fopenmp -I$(OBJDIR) -J$(OBJDIR)
FFLAGS += -fconvert=big-endian # For dumping .vtk files, add -fconvert=big-endian 


#----FFTW library------------------- (adjust prefix if needed)--------------
# For personal Linux machine
FFTW_PREFIX := $(shell pkg-config --variable=prefix fftw3)

# For personal MAC machine
#FFTW_PREFIX := $(shell brew --prefix fftw)

#--- End prefix modification -----------------------------------------------

FFLAGS += -I$(FFTW_PREFIX)/include 
LIBS = -L$(FFTW_PREFIX)/lib  -lfftw3_omp -lfftw3 -lm

#----- END FFTW library-------------------

# Source files
SRC =	$(SRCDIR)/main.f90 \
		$(SRCDIR)/memAlloc.f90 \
		$(SRCDIR)/grid.f90 \
		$(SRCDIR)/init_FFT.f90 \
		$(SRCDIR)/initialCondition.f90 \
		$(SRCDIR)/rhsVelocity.f90 \
		$(SRCDIR)/pressurePoisson.f90 \
		$(SRCDIR)/fileIO.f90 \

MOD_FILES = $(SRCDIR)/my_fftw.f90 \
			$(SRCDIR)/params.f90 \
			$(SRCDIR)/ghost.f90

# Object files
OBJ = $(addprefix $(OBJDIR)/, $(notdir $(SRC:.f90=.o)))
MOBJ = $(addprefix $(OBJDIR)/, $(notdir $(MOD_FILES:.f90=.o)))

# Default target
all: $(OUTDIR) $(OBJDIR) $(TARGET)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OUTDIR):
	mkdir -p $(OUTDIR)

# Pattern rule to compile .f90 files into .o files
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Linking the final executable
$(TARGET): $(MOBJ) $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(MOBJ) $(OBJ) $(LIBS)

# Clean up generated files
clean:
	@rm -f $(OBJDIR)/*.o $(TARGET)
	@rm -f $(OBJDIR)/*.mod

veryclean: clean
	@rm -rf $(OUTDIR)
	@rm -rf $(OBJDIR)

.PHONY : clean veryclean