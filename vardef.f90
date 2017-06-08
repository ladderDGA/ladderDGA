!in this module the variables used in the main program are defined
MODULE vardef
  IMPLICIT NONE

  !MPI-variables
  INTEGER :: myid, nprocs,ierror,request
  INTEGER, DIMENSION(:), ALLOCATABLE :: sendstatus,recvstatus

  REAL(KIND=8), PARAMETER :: pi=DACOS(-1d0)
  
  !variables for file-names
  CHARACTER(LEN=:), ALLOCATABLE :: fname

  !parameters read in from ladderDGA.into 
  REAL(KIND=8) :: uhub,mu,beta,nden,xch_so,xsp_so
  INTEGER :: Iwbox,Iwbox_bose,shift,LQ,Nint,k_number
  LOGICAL :: sigma_only,chi_only,lambdaspin_only,sumallch,sumallsp

  !index of bosonic matsubara frequency (1 rank = 1 frequency, might be changed)
  INTEGER :: i

  !local green's functionns and vertices
  COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: gww,self,selfloc,selfloc_res
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: fupdown,gammach,gammasp
  
  !local susceptibilities and sum indices
  COMPLEX(KIND=8) :: chich_loc,chisp_loc,chich_loc_sum,chisp_loc_sum
  COMPLEX(KIND=8) :: chich_loc_all,chisp_loc_all
  REAL(KIND=8) :: ind_part_ch,ind_part_sp,ind_ch,ind_sp
  !tolerancech and tolerancesp are set 
  REAL(KIND=8) :: tolerancech=1d-7, tolerancesp=1d-7
  INTEGER :: sum_ind_ch,sum_ind_sp

  !nonlocal susceptibilities
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: chi_bubble,chi_bubble_pos,chi_bubble_neg
  COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: chich_x0,chisp_x0
  COMPLEX(KIND=8) :: chich_q_sum,chich_sum

  !variables for grids in momentum space
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: klist,dcok,dsik,dcol,dsil
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dcoq,dsiq,Qv
  REAL(KIND=8) :: qx,qy,qmax,Q0b

  !variables for controlling loops
  INTEGER :: ix,iy,j,ind
  REAL(KIND=8) :: a,b

  !variables for lambda-corrections
  REAL(KIND=8) :: lambdach,lambdasp
  !chi_inv_min is initialized
  COMPLEX(KIND=8) :: chich_inv_min=dcmplx(1d-100,0d0),chisp_inv_min=dcmplx(1d-100,0d0)

  !variables for DGA self-energy
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: selflist,selflist_res
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: selflistch,selflistch_res
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: selflistsp,selflistsp_res
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: selflistrest,selflistrest_res

END MODULE vardef
