!
! pGFNFF library:
! ---------------
!
! Copyright (C) 2020-2022 Julian Gale
!
! pGFNFF is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! pGFNFF is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with pGFNFF.  If not, see <https://www.gnu.org/licenses/>.
!
  subroutine pgfnff_init_param_r(param,value,ierror)
!
!  Initialises parameters that are of type double precision
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  use m_pgfnff_cfg
  use m_pgfnff_disp
  use m_pgfnff_topo
  implicit none
!
!  Passed arguments
!
  character(len=*),     intent(in)  :: param
  real(dp),             intent(in)  :: value
  integer(i4),          intent(out) :: ierror
!
  ierror = 0
!
!  Search for parameter whose value has been passed in
!
  if (index(param,'tdist_thr').ne.0) then
    tdist_thr = value
  elseif (index(param,'atm_alpha1').ne.0) then
    atm_alpha1 = value
  elseif (index(param,'q_trap').ne.0) then
    gfnff_q_trap = value
  elseif (index(param,'pi_temp1').ne.0) then
    gfnff_pi_temp1 = value
  elseif (index(param,'pi_temp2').ne.0) then
    gfnff_pi_temp2 = value
  elseif (index(param,'ks').ne.0) then
    gfnff_ks = value
  elseif (index(param,'max_accuracy').ne.0) then
    max_gfnff_accuracy = value
  elseif (index(param,'max_acc_disp').ne.0) then
    max_gfnff_acc_disp = value
  elseif (index(param,'max_acc_rep').ne.0) then
    max_gfnff_acc_rep = value
  elseif (index(param,'max_acc_cn').ne.0) then
    max_gfnff_acc_cn = value
  elseif (index(param,'max_acc_hb1').ne.0) then
    max_gfnff_acc_hb1 = value
  elseif (index(param,'max_acc_hb2').ne.0) then
    max_gfnff_acc_hb2 = value
  elseif (index(param,'accuracy_overall').ne.0) then
    gfnff_accuracy = value
  elseif (index(param,'accuracy_disp').ne.0) then
    gfnff_accuracy_disp = value
  elseif (index(param,'accuracy_rep').ne.0) then
    gfnff_accuracy_rep = value
  elseif (index(param,'accuracy_cn').ne.0) then
    gfnff_accuracy_cn = value
  elseif (index(param,'accuracy_hb1').ne.0) then
    gfnff_accuracy_hb1 = value
  elseif (index(param,'accuracy_hb2').ne.0) then
    gfnff_accuracy_hb2 = value
  elseif (index(param,'cnc6tol').ne.0) then
    gfnff_cnc6tol = value
  elseif (index(param,'taper').ne.0) then
    gfnff_taper = value
  elseif (index(param,'wolf_eta').ne.0) then
    gfnff_wolf_eta = value
  elseif (index(param,'rtopo').ne.0) then
    rfgoed1 = value
  else
    ierror = -1
  endif
!
  end subroutine pgfnff_init_param_r

  subroutine pgfnff_init_param_i(param,value,ierror)
!
!  Initialises parameters that are of integer type 
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  use m_pgfnff_cfg
  use m_pgfnff_disp
  use m_pgfnff_topo
  implicit none
! 
!  Passed arguments
!
  character(len=*),     intent(in)  :: param
  integer(i4),          intent(in)  :: value
  integer(i4),          intent(out) :: ierror
!
  ierror = 0
!
!  Search for parameter whose value has been passed in
!
  if (index(param,'maxtoposhell1').ne.0) then
    maxtoposhell1 = value
  elseif (index(param,'maxtoposhell').ne.0) then
    maxtoposhell = value
  elseif (index(param,'pi_change').ne.0) then
    gfnff_pi_change = value
  else
    ierror = -1
  endif
!
  end subroutine pgfnff_init_param_i

  subroutine pgfnff_init_param_l(param,value,ierror)
!
!  Initialises parameters that are of logical type 
!
!  Julian Gale, Curtin University, December 2021
!
  use m_pgfnff_types
  use m_pgfnff
  use m_pgfnff_cfg
  use m_pgfnff_disp
  use m_pgfnff_topo
  implicit none
! 
!  Passed arguments
!
  character(len=*),     intent(in)  :: param
  logical,              intent(in)  :: value
  integer(i4),          intent(out) :: ierror
!
  ierror = 0
!
!  Search for parameter whose value has been passed in
!
  if (index(param,'newtopo').ne.0) then
    lgfnff_newtopo = value
  elseif (index(param,'xtbtopo').ne.0) then
    lgfnff_xtbtopo = value
  elseif (index(param,'highcn_trap').ne.0) then
    lgfnff_highcn_trap = value
  elseif (index(param,'fragment_bond').ne.0) then
    lgfnff_fragment_bond = value
  elseif (index(param,'topowolf').ne.0) then
    lgfnff_topowolf = value
  else
    ierror = -1
  endif
!
  end subroutine pgfnff_init_param_l
