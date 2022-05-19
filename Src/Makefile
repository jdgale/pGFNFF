################################################################################ #
#
#  pGFNFF
#
################################################################################

default: lib

#
#  Definitions
#
LIBS    = 
FCFLAGS	= -O

FC	= gfortran
F90	= $(FC) $(FCFLAGS) -c
ARCHIVE = ar -rv

LOBJECTS= lapack.o blas.o dlamch.o
FOBJECTS= m_pgfnff_types.o m_pgfnff_reallocate.o m_io.o \
          m_pgfnff.o m_pgfnff_nbr_lib.o m_pgfnff_mrec.o m_pgfnff_cfg.o m_pgfnff_disp.o m_pgfnff_topo.o \
          pgfnff_dcn.o pgfnff_init.o pgfnff_init_param.o pgfnff_radij.o pgfnff_hybrid.o pgfnff_setnbr.o \
          pgfnff_subs.o pgfnff_solve.o pgfnff_matsolve.o pgfnff_huckel.o pgfnff_huckel_pbc.o pgfnff_outpot.o \
          pgfnff_outpar.o pgfnff_excludebond.o pgfnff_get_max_cutoff.o pgfnff_check_atom_with_nobonds.o \
          pgfnff_pargen.o pgfnff_sethkpt.o fermi.o pgfnff_error.o pgfnff_warn.o changemaxatpgfnff.o pgfnff_query.o
MOBJECTS= main.o m_nbr.o 
OBJECTS = $(FOBJECTS) $(LOBJECTS)

exe:	$(MOBJECTS) $(OBJECTS) ; $(ARCHIVE) libpGFNFF.a $(OBJECTS); $(FC) -o pGFNFF main.o m_nbr.o libpGFNFF.a

lib:    $(FOBJECTS)
	$(ARCHIVE) libpGFNFF.a $(FOBJECTS)

clean:           ; rm -f $(OBJECTS) *.mod
libclean:        ; rm -f $(OBJECTS) *.mod ; rm -f ./libpGFNFF.a
tar  :           ; cd ../ ; tar -cvzf pGFNFF.tgz README Src/*.F90 \
                   Src/Makefile \

#
#  pGFNFF library routines
#
changemaxatpgfnff.o              : changemaxatpgfnff.F90              ; $(F90) changemaxatpgfnff.F90
m_pgfnff.o                       : m_pgfnff.F90                       ; $(F90) m_pgfnff.F90
m_pgfnff_cfg.o                   : m_pgfnff_cfg.F90                   ; $(F90) m_pgfnff_cfg.F90
m_pgfnff_disp.o                  : m_pgfnff_disp.F90                  ; $(F90) m_pgfnff_disp.F90
m_pgfnff_mrec.o                  : m_pgfnff_mrec.F90                  ; $(F90) m_pgfnff_mrec.F90
m_pgfnff_nbr_lib.o               : m_pgfnff_nbr_lib.F90               ; $(F90) m_pgfnff_nbr_lib.F90
m_pgfnff_topo.o                  : m_pgfnff_topo.F90                  ; $(F90) m_pgfnff_topo.F90
pgfnff_check_atom_with_nobonds.o : pgfnff_check_atom_with_nobonds.F90 ; $(F90) pgfnff_check_atom_with_nobonds.F90
pgfnff_dcn.o                     : pgfnff_dcn.F90                     ; $(F90) pgfnff_dcn.F90
pgfnff_excludebond.o             : pgfnff_excludebond.F90             ; $(F90) pgfnff_excludebond.F90
pgfnff_get_max_cutoff.o          : pgfnff_get_max_cutoff.F90          ; $(F90) pgfnff_get_max_cutoff.F90
pgfnff_huckel.o                  : pgfnff_huckel.F90                  ; $(F90) pgfnff_huckel.F90
pgfnff_huckel_pbc.o              : pgfnff_huckel_pbc.F90              ; $(F90) pgfnff_huckel_pbc.F90
pgfnff_hybrid.o                  : pgfnff_hybrid.F90                  ; $(F90) pgfnff_hybrid.F90
pgfnff_init.o                    : pgfnff_init.F90                    ; $(F90) pgfnff_init.F90
pgfnff_init_param.o              : pgfnff_init_param.F90              ; $(F90) pgfnff_init_param.F90
pgfnff_matsolve.o                : pgfnff_matsolve.F90                ; $(F90) pgfnff_matsolve.F90
pgfnff_outpar.o                  : pgfnff_outpar.F90                  ; $(F90) pgfnff_outpar.F90
pgfnff_outpot.o                  : pgfnff_outpot.F90                  ; $(F90) pgfnff_outpot.F90
pgfnff_pargen.o                  : pgfnff_pargen.F90                  ; $(F90) pgfnff_pargen.F90
pgfnff_query.o                   : pgfnff_query.F90                   ; $(F90) pgfnff_query.F90
pgfnff_radij.o                   : pgfnff_radij.F90                   ; $(F90) pgfnff_radij.F90
pgfnff_sethkpt.o                 : pgfnff_sethkpt.F90                 ; $(F90) pgfnff_sethkpt.F90
pgfnff_setnbr.o                  : pgfnff_setnbr.F90                  ; $(F90) pgfnff_setnbr.F90
pgfnff_setup.o                   : pgfnff_setup.F90                   ; $(F90) pgfnff_setup.F90
pgfnff_solve.o                   : pgfnff_solve.F90                   ; $(F90) pgfnff_solve.F90
pgfnff_subs.o                    : pgfnff_subs.F90                    ; $(F90) pgfnff_subs.F90
#
#  Dependencies for the above
#
m_pgfnff.o: m_pgfnff.F90 m_pgfnff_types.o m_io.o
m_pgfnff_cfg.o: m_pgfnff_cfg.F90 m_pgfnff_types.o m_pgfnff_reallocate.o m_io.o
m_pgfnff_disp.o: m_pgfnff_disp.F90 m_pgfnff_types.o m_pgfnff.o m_io.o
m_pgfnff_mrec.o: m_pgfnff_mrec.F90 m_pgfnff_types.o m_pgfnff_nbr_lib.o m_io.o
m_pgfnff_nbr_lib.o: m_pgfnff_nbr_lib.F90 m_pgfnff_types.o m_pgfnff_reallocate.o m_pgfnff_cfg.o m_pgfnff.o m_io.o
m_pgfnff_topo.o: m_pgfnff_topo.F90 m_pgfnff_types.o m_pgfnff_reallocate.o m_pgfnff_cfg.o m_io.o
changemaxatpgfnff.o: changemaxatpgfnff.F90 m_pgfnff_types.o m_pgfnff_cfg.o m_io.o
pgfnff_check_atom_with_nobonds.o: pgfnff_check_atom_with_nobonds.F90 m_pgfnff_types.o m_pgfnff.o
pgfnff_dcn.o: pgfnff_dcn.F90 m_pgfnff_types.o m_pgfnff.o
pgfnff_excludebond.o: pgfnff_excludebond.F90 m_pgfnff_types.o
pgfnff_get_max_cutoff.o: pgfnff_get_max_cutoff.F90 m_pgfnff_types.o m_io.o m_pgfnff.o
pgfnff_huckel.o: pgfnff_huckel.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_cfg.o m_io.o
pgfnff_huckel_pbc.o: pgfnff_huckel_pbc.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_cfg.o m_io.o
pgfnff_hybrid.o: pgfnff_hybrid.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_nbr_lib.o m_io.o
pgfnff_init.o: pgfnff_init.F90 m_pgfnff_types.o m_pgfnff.o 
pgfnff_init_param.o: pgfnff_init_param.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_cfg.o m_pgfnff_topo.o m_pgfnff_disp.o
pgfnff_matsolve.o: pgfnff_matsolve.F90 m_pgfnff_types.o m_io.o
pgfnff_outpar.o: pgfnff_outpar.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_disp.o m_pgfnff_topo.o
pgfnff_outpot.o: pgfnff_outpot.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_disp.o m_pgfnff_topo.o m_io.o
pgfnff_pargen.o: pgfnff_pargen.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_nbr_lib.o m_pgfnff_mrec.o m_pgfnff_disp.o m_pgfnff_topo.o m_io.o
pgfnff_query.o: pgfnff_query.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_cfg.o m_pgfnff_nbr_lib.o
pgfnff_radij.o: pgfnff_radij.F90 m_pgfnff_types.o m_pgfnff.o
pgfnff_sethkpt.o: pgfnff_sethkpt.F90 m_pgfnff_types.o m_pgfnff_cfg.o m_pgfnff.o m_io.o
pgfnff_setnbr.o: pgfnff_setnbr.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_nbr_lib.o
pgfnff_solve.o: pgfnff_solve.F90 m_pgfnff_types.o
pgfnff_subs.o: pgfnff_subs.F90 m_pgfnff_types.o m_pgfnff.o m_pgfnff_cfg.o m_pgfnff_topo.o m_io.o
#
#  Example calling routines for library
#
main.o                    : main.F90                    ; $(F90) main.F90
m_nbr.o                   : m_nbr.F90                   ; $(F90) m_nbr.F90
#
#  Dependencies for the above
#
m_nbr.o        : m_nbr.F90 m_pgfnff_types.o m_pgfnff_reallocate.o
main.o         : main.F90 m_pgfnff_types.o m_io.o m_nbr.o 
#
#  General purpose subroutines
#
m_pgfnff_types.o          : m_pgfnff_types.F90          ; $(F90) m_pgfnff_types.F90
m_pgfnff_reallocate.o     : m_pgfnff_reallocate.F90     ; $(F90) m_pgfnff_reallocate.F90
m_io.o                    : m_io.F90                    ; $(F90) m_io.F90
blas.o                    : blas.F90                    ; $(F90) blas.F90
dlamch.o                  : dlamch.F90                  ; $(F90) dlamch.F90
fermi.o                   : fermi.F90                   ; $(F90) fermi.F90
lapack.o                  : lapack.F90                  ; $(F90) lapack.F90
pgfnff_error.o            : pgfnff_error.F90            ; $(F90) pgfnff_error.F90
pgfnff_warn.o             : pgfnff_warn.F90             ; $(F90) pgfnff_warn.F90
#
#  Dependencies for the above
#
m_pgfnff_types.o: m_pgfnff_types.F90
m_io.o: m_io.F90 m_pgfnff_types.o 
blas.o: blas.F90 m_pgfnff_types.o
dlamch.o: dlamch.F90 m_pgfnff_types.o
fermi.o: fermi.F90 m_pgfnff_types.o
lapack.o: lapack.F90 m_pgfnff_types.o
pgfnff_error.o: pgfnff_error.F90 m_pgfnff_types.o m_io.o
pgfnff_warn.o: pgfnff_warn.F90 m_pgfnff_types.o m_io.o
m_pgfnff_reallocate.o: m_pgfnff_reallocate.F90 m_pgfnff_types.o 
