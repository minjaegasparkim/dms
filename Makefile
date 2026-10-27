#====================================================================================================
# DMS - Makefile
#
#   make              optimised program with OpenMP          -> bin/dms
#   make debug        run-time checks, no optimisation        -> bin/dms
#   make serial       optimised, without OpenMP
#   make python       Python module (f2py)                    -> python/dms/_dms*.so
#   make test         build both and run the test suite (needs pytest)
#   make clean        remove build files
#
# Other compilers:  make FC=ifx OMP=-qopenmp
#====================================================================================================

FC      ?= gfortran
ifeq ($(FC),f77)
  FC = gfortran
endif
PYTHON  ?= python3

OPT     ?= -O3 -march=native -funroll-loops
OMP     ?= -fopenmp
STD     = -std=f2008 -ffree-line-length-none
PIC     = -fPIC
DEBUG   = -O0 -g -fcheck=all -fbacktrace -Wall -Wextra -Wno-unused-dummy-argument -Wno-compare-reals \
          -ffpe-trap=zero,overflow -finit-real=snan

BUILD   = build
BIN     = bin
SRC     = src

# modules (dependency order), then the routines; dms.f90 is the program
MODS    = datatype constants parameters user_density functions mie_routines radiation_calculation \
          volume_splitting fits_io output_routines obs_tools
SUBS    = input_setup read_input_file input_questions dms_driver dms_reset dust_mixing pre_cals \
          scattering_angles_cals reading_spectrum reading_temperature_dependent_data mie_cals \
          Qabs_Qsca_S11_S12_cals temperature_distribution_cals particle_abundances belt_normalisation \
          max_temp_diff_cals computing_images observation_products computing_SED optical_depth_cals summary

LIBOBJ  = $(addprefix $(BUILD)/, $(addsuffix .o, $(MODS) $(SUBS)))
LIB     = $(BUILD)/libdms.a

FFLAGS  = $(OPT) $(OMP) $(STD) $(PIC) -J$(BUILD)

.PHONY: all debug serial python test clean

all: $(BIN)/dms

$(BIN)/dms: $(LIB) $(BUILD)/dms.o | $(BIN)
	$(FC) $(FFLAGS) -o $@ $(BUILD)/dms.o $(LIB)

$(LIB): $(LIBOBJ)
	ar rcs $@ $^

# modules must be compiled in order (each depends on the previous ones)
$(BUILD)/%.o: $(SRC)/%.f90 | $(BUILD)
	$(FC) $(FFLAGS) -c $< -o $@

$(BUILD)/constants.o:             $(BUILD)/datatype.o
$(BUILD)/parameters.o:            $(BUILD)/constants.o
$(BUILD)/user_density.o:          $(BUILD)/constants.o
$(BUILD)/functions.o:             $(BUILD)/parameters.o $(BUILD)/user_density.o
$(BUILD)/mie_routines.o:          $(BUILD)/datatype.o
$(BUILD)/radiation_calculation.o: $(BUILD)/functions.o
$(BUILD)/volume_splitting.o:      $(BUILD)/radiation_calculation.o
$(BUILD)/fits_io.o:               $(BUILD)/datatype.o
$(BUILD)/output_routines.o:       $(BUILD)/functions.o $(BUILD)/fits_io.o
$(BUILD)/obs_tools.o:             $(BUILD)/constants.o
$(addprefix $(BUILD)/, $(addsuffix .o, $(SUBS))) $(BUILD)/dms.o: \
    $(addprefix $(BUILD)/, $(addsuffix .o, $(MODS)))

$(BUILD) $(BIN):
	mkdir -p $@

debug:
	$(MAKE) clean
	$(MAKE) OPT="$(DEBUG)"

serial:
	$(MAKE) clean
	$(MAKE) OMP=""

python: $(LIB)
	cd python && FC=$(FC) $(PYTHON) build_extension.py --lib ../$(LIB) --mod ../$(BUILD) --omp="$(OMP)"

test: $(BIN)/dms python
	$(PYTHON) -m pytest tests -q

clean:
	rm -rf $(BUILD) $(BIN)/dms python/dms/_dms*.so python/build_f2py
