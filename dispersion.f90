!this module contains the dispersion relation
MODULE dispersion
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: ng=1
  !points for the Gauss-Legendre integration
  REAL(KIND=8), PARAMETER, DIMENSION(ng) :: tstep=(/1.0d0/)
  !weights for the Gauss-Legendre integration
  REAL(KIND=8), PARAMETER, DIMENSION(ng) :: ws=(/2.0d0/)

  !here data for higher order Gauss-Legendre integration can be given
  !the parameter ng has to be set accordingly

  !data for 5th-order Legendre Integration
  !REAL(KIND=8), PRIVATE, PARAMETER, DIMENSION(ng) :: tstep= &
  !(/-0.9061798459386640,-0.53846931010568309, &
  !0.0000000000000000, 0.53846931010568309,&
  !0.9061798459386640/)
  !REAL(KIND=8), PRIVATE, PARAMETER, DIMENSION(ng) :: ws= &
  !(/ 0.2369268850561891, 0.47862867049936650, &
  !0.5688888888888889, 0.47862867049936650, &
  !0.2369268850561891/)
  
CONTAINS

  SUBROUTINE init_arrays(Nint,LQ,k_number,klist,qmin,qmax,q,Q0b,Qv,dcok,dsik,dcoq,dsiq,dcol,dsil)
    IMPLICIT NONE
    !calculates the cos() and sin() for all k', q - and k-values of the corresponding grids

    !input: qpoint where chi(q) is maximal, grid-defining parameters and external k-points
    REAL(KIND=8), INTENT(IN) :: qmin,qmax,q
    INTEGER, INTENT(IN) :: Nint,LQ,k_number
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: klist
    !output: q-list, cos() and sin() of all k-points
    REAL(KIND=8), DIMENSION(0:), INTENT(OUT) :: Qv
    REAL(KIND=8), INTENT(OUT) :: Q0b
    REAL(KIND=8), DIMENSION(0:,:), INTENT(OUT) :: dcok,dsik
    REAL(KIND=8), DIMENSION(0:), INTENT(OUT) :: dcoq,dsiq
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: dcol,dsil
    !subroutine internal variables
    INTEGER :: j,jg
    REAL(KIND=8) :: k,pi

    pi=dacos(-1.0d0)

    DO j=0,Nint-1
       DO jg=1,Ng
          k=(pi/dfloat(Nint))*(tstep(jg)+ &
               2.0d0*dfloat(j)+1.0d0)-pi
          dcok(j,jg)=dcos(k)
          dsik(j,jg)=dsin(k)
       ENDDO
    ENDDO

    Q0b=(qmax-qmin)/dfloat(LQ-1) 
    DO j=0,LQ-1
       Qv(j)=qmin+dfloat(j)*Q0b
       IF((Qv(j)-Q0b/2d0) .LT. q &
            .AND. (Qv(j)+Q0b/2d0) .GT. q) THEN
          Qv(j)=q
       END IF
       dcoq(j)=dcos(Qv(j))
       dsiq(j)=dsin(Qv(j))
    ENDDO

    DO j=1,k_number
       dcol(j,1)=dcos(klist(j,1))
       dcol(j,2)=dcos(klist(j,2))
       dsil(j,1)=dsin(klist(j,1))
       dsil(j,2)=dsin(klist(j,2))
    ENDDO
  END SUBROUTINE init_arrays

  REAL(KIND=8) PURE FUNCTION eps(ax,ay,bx,by,cx,cy,dx,dy)
    IMPLICIT NONE
    !calculates the dispersion in terms of sin() and cos() functions

    !input: arguments for the dispersion are the already evaluated sin and cos functions
    REAL(KIND=8), INTENT(IN) :: ax,ay,bx,by,cx,cy,dx,dy
    !subroutine internal variables
    !bandwidth: 2t/(2t*sqrt(4)) for 2D
    REAL(KIND=8), PARAMETER :: tsc=0.5d0

    !dispersion for the simple cubic lattice
    eps = -tsc*(ax*bx-cx*dx+ &
         ay*by-cy*dy)

  END FUNCTION eps

END MODULE dispersion
  

