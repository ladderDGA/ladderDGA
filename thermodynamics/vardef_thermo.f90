!in this module the variables used in the main program thermodynamics.f90 are defined
MODULE vardef_thermo
  IMPLICIT NONE

  !MPI-variables
  INTEGER :: myid, nprocs,ierror

  REAL(KIND=8), PARAMETER :: pi=dacos(-1d0)

  !parameters read in from ladderDGA_thermo.in
  REAL(KIND=8) :: uhub,mu,beta,nden,fermicut,epsmin,epsmax
  INTEGER :: Iwbox,k_range,epssteps,ap
  LOGICAL :: calcU0,calcDMFT,calcDGA


  !Auxilliary variables
  REAL(KIND=8) :: sigma_hartree,epoint
  COMPLEX(KIND=8) :: g,w

  !Anderson parameters of DMFT
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v
  REAL(KIND=8) :: vsum

  !local selfenergy and green's function (needed for DMFT only)
  COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: self_loc,gww

  !non-local selfenergy (needed for DGA only)
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: self_k_tot
  
  !variables for grids in momentum space
  INTEGER :: knumber_tot,knumber,krest,koffset,k_min,k_max
  INTEGER, DIMENSION(:), ALLOCATABLE :: recvnumber,offsetnum
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: kcount
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: klist_tot,klist,dcol

  !variables for controlling loops
  INTEGER :: i,j,ix,iy,ind,ieps

  !self-energy for the energy-subroutine: 1) 0 for U=0, 2) sigma_dmft for DMFT 3) sigma_DGA
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: self

  !for output: energies, occupation numbers
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: en,en_sum,n_eps,n_eps_sum

  !for testing the k-sums
  REAL(KIND=8) :: testsum,testsum_sum,testsum_eps,testsum_eps_sum

END MODULE vardef_thermo
