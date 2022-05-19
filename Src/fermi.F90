  subroutine fermi(nk, wk, maxe, ne, e, kbt, qtot, wke, efermi, ehomo, elumo)
!
!  Finds the Fermi energy and the occupation weights of states.
!
!  Based on the routine fermid from SIESTA originally written by J.M.Soler.
!
!  On entry:
!  
!  nk           = number of K points
!  wk(nk)       = Weights of K points
!  maxe         = First dimension of eigenvalue and occupation arrays
!  ne           = Number of bands
!  e(maxe,nk)   = Eigenvalues for K points
!  kbt          = Temperature x Boltzmann constant in energy units
!  qtot         = Total number of electrons in bands
!
!  On exit:
!
!  wke(maxe,nk) = Occupations multiplied by k point weights
!  efermi       = Fermi energy
!  ehomo        = HOMO energy
!  elumo        = LUMO energy
!
  use m_pgfnff_types
  implicit none
!
!  Passed variables
!
  integer(i4),                     intent(in)   :: maxe
  integer(i4),                     intent(in)   :: ne
  integer(i4),                     intent(in)   :: nk
  real(dp),                        intent(in)   :: e(maxe,nk)
  real(dp),                        intent(in)   :: qtot
  real(dp),                        intent(in)   :: kbt
  real(dp),                        intent(in)   :: wk(nk)
  real(dp),                        intent(out)  :: wke(maxe,nk)
  real(dp),                        intent(out)  :: efermi
  real(dp),                        intent(out)  :: ehomo
  real(dp),                        intent(out)  :: elumo
!
!  Local variables
!
  integer(i4)                                   :: ie
  integer(i4)                                   :: ik
  integer(i4)                                   :: iter
  integer(i4)                                   :: nel
  integer(i4),                             save :: nitmax = 150
  real(dp),                                save :: tol = 1.0d-10
  real(dp)                                      :: defermi
  real(dp)                                      :: drange
  real(dp)                                      :: dstepf
  real(dp)                                      :: dsumq
  real(dp)                                      :: emin
  real(dp)                                      :: emax
  real(dp)                                      :: stepf
  real(dp)                                      :: sumq
  real(dp)                                      :: t
  real(dp)                                      :: x
!
!  Find number of electrons in each level
!
  nel = nint(qtot)
!
!  Initialise occupancies
!
  wke(1:ne,1:nk) = 0.0_dp
  wke(1:nel,1:nk) = 1.0_dp
!
!  Determine Fermi level
!
  sumq = 0.0_dp
  emin = e(1,1)
  emax = e(1,1)
  do ik = 1,nk
    do ie = 1,ne
      wke(ie,ik) = wk(ik)
      sumq = sumq + wke(ie,ik)
      emin = min(emin,e(ie,ik))
      emax = max(emax,e(ie,ik))
    enddo
  enddo
!
  if (nel.lt.ne) then
    efermi = 0.5_dp*(e(nel,1) + e(nel+1,1))
  else
    efermi = emax
  endif
!
  if (abs(sumq-qtot).lt.tol) then
    goto 99
  endif
  if (sumq.lt.qtot) then
    call pgfnff_error('insufficient states to find Fermi energy',0_i4)
    stop
  endif
!
  T = max(kbt,1.d-6)
  drange = T*sqrt(-log(tol*0.01_dp))
  emin = emin - drange
  emax = emax + drange
  do iter = 1,nitmax
    sumq = 0.0_dp
    dsumq = 0.0_dp
    do ik = 1,nk
      do ie = 1,ne
!
!  Fermi-Dirac distribution
!
        x = (e(ie,ik)-efermi)/T
        if (x.gt.100.0_dp) then
          stepf = 0.0_dp
          dstepf = 0.0_dp
        elseif (x.lt.-100.0_dp) then
          stepf = 1.0_dp
          dstepf = 0.0_dp
        else
          stepf = 1.0_dp / ( 1.0_dp + exp(x) )
          dstepf = exp(x) / &
              &       (T*(exp(x)+1.0_dp)**2)
        endif
        wke(ie,ik) = wk(ik)*stepf
        sumq = sumq + wke(ie,ik)
        dsumq = dsumq + wk(ik)*dstepf
      enddo
    enddo
!
    if (abs(dsumq).gt.1.0d-12) then
      defermi = (qtot - sumq)/dsumq
      efermi = efermi + defermi
    endif
!
!  If the Fermi level was found
!
    if (abs(sumq-qtot).lt.tol) goto 99
  enddo
!
  call pgfnff_error('calculation of Fermi energy has failed to converge',0_i4)
  stop
!
!  Find HOMO and LUMO
!
  99 continue
  ehomo = - 1.0d8
  elumo =   1.0d8
  do ik = 1,nk
    do ie = 1,ne
      if (e(ie,ik).gt.efermi) then
        elumo =  min(e(ie,ik),elumo)
      elseif (e(ie,ik).lt.efermi) then
        ehomo =  max(e(ie,ik),ehomo)
      endif
    enddo
  enddo

  end subroutine fermi
