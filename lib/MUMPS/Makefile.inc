LPORDDIR = $(MUMPS_DIR)/PORD/lib
LORDER  = $(LPORDDIR)/libpord.a
IPOORD  = -I$(MUMPS_DIR)/PORD/include

LMETISDIR =
IMETIS    =

LMETIS    =

LIBBLAS = -lopenblas
LIBOTHERS = -lpthread

CC      = gcc
FC      = gfortran
FL      = gfortran
AR      = ar cvr
RANLIB  = ranlib
SCALAP  =

INCS = $(IPOORD)

OPTF    = -O2 -fallow-argument-mismatch
OPTL    = -O2
OPTC    = -O2

ORDERINGS = -Dpord
