!Contains the subroutines for calculating:
!-) non-interacting susceptibility (bare bubble) -> in the programm: chi_bubble(nu)
!-) full interacting susceptibility (in ladder approximation) -> in the program: chisp, und chich
!-) for a given omega and q!
MODULE calc_susc
  USE dispersion
  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: qmaxsteps=100
  REAL(KIND=8), PARAMETER, PRIVATE :: qmaxprec=1e-14

CONTAINS

  SUBROUTINE calc_bubble(mu,beta,Iwbox,i,self,LQ,Nint,dcok,dsik,dcoq,dsiq, &
       chi_bubble)
    !caluclates the non-interacing susceptibility as a function of (nu,omega,q)

    !input:
    INTEGER, INTENT(IN) :: LQ,Nint,i,Iwbox
    REAL(KIND=8), INTENT(IN) :: mu,beta
    REAL(KIND=8), DIMENSION(0:,:), INTENT(IN) :: dcok,dsik
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: dcoq,dsiq
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: self
    !output:
    COMPLEX(KIND=8), DIMENSION(0:,:), INTENT(OUT) :: chi_bubble
    !subroutine internal variables
    INTEGER :: j,ix,iy,iNx,iNy,igx,igy,ind
    REAl(KIND=8) :: pi,ek,ekq
    COMPLEX(KIND=8), ALLOCATABLE :: w(:)

    pi=dacos(-1.0d0)
    ALLOCATE(w(-2*Iwbox:2*Iwbox-1))
    DO j=-2*Iwbox,2*Iwbox-1
       w(j)=dcmplx(0.0d0,(pi/beta)*dfloat(2*j+1))+mu-self(j)
    ENDDO
    
    chi_bubble=dcmplx(0.0d0,0.0d0)

    DO igx=1, Ng
       DO igy=1, Ng
          DO iNx=0,Nint-1
             DO iNy=0,Nint-1
                ek=eps(dcok(iNx,igx),dcok(iNy,igy), &
                     1.d0,1.d0,0.d0,0.d0,0.d0,0.d0)
                DO ix=0,LQ-1
                   DO iy=0,ix
                      ind=ix*(ix+1)/2+iy+1
                      ekq=eps(dcok(iNx,igx),dcok(iNy,igy), &
                           dcoq(ix),dcoq(iy),dsik(iNx,igx),dsik(iNy,igy), &
                           dsiq(ix),dsiq(iy))
                      DO j=0,iwbox-1
                         chi_bubble(j,ind)=chi_bubble(j,ind)- &
                              ws(igx)*ws(igy)/((w(j)-ek)*(w(j+i)-ekq))
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    chi_bubble=beta*chi_bubble/ &
         ((2.0d0*dfloat(Nint))**2)
    
    DEALLOCATE(w)

  END SUBROUTINE calc_bubble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE calc_bubble_slice(mu,beta,Iwbox,i,self,LQ,Nint,dcok,dsik,dcoq,dsiq, &
       chi_bubble_slice)
    !caluclates the non-interacing susceptibility (q=(pi,qy)) as a function of (nu,omega,qy)

    !input:
    INTEGER, INTENT(IN) :: LQ,Nint,i,Iwbox
    REAL(KIND=8), INTENT(IN) :: mu,beta
    REAL(KIND=8), DIMENSION(0:,:), INTENT(IN) :: dcok,dsik
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: dcoq,dsiq
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: self
    !output:
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,:), INTENT(OUT) :: chi_bubble_slice
    !subroutine internal variables
    INTEGER :: j,ix,iy,iNx,iNy,igx,igy,ind
    REAl(KIND=8) :: pi,ek,ekq
    COMPLEX(KIND=8), ALLOCATABLE :: w(:)

    pi=dacos(-1.0d0)
    ALLOCATE(w(-2*Iwbox:2*Iwbox-1))
    DO j=-2*Iwbox,2*Iwbox-1
       w(j)=dcmplx(0.0d0,(pi/beta)*dfloat(2*j+1))+mu-self(j)
    ENDDO
    
    chi_bubble_slice=dcmplx(0.0d0,0.0d0)

    ix=LQ-1
    DO igx=1, Ng
       DO igy=1, Ng
          DO iNx=0,Nint-1
             DO iNy=0,Nint-1
                ek=eps(dcok(iNx,igx),dcok(iNy,igy), &
                     1.d0,1.d0,0.d0,0.d0,0.d0,0.d0)
                DO iy=0,ix
                   ekq=eps(dcok(iNx,igx),dcok(iNy,igy), &
                        -1d0,dcoq(iy),dsik(iNx,igx),dsik(iNy,igy), &
                        0.0d0,dsiq(iy))
                   DO j=-iwbox,iwbox-1
                      chi_bubble_slice(j,iy+1)=chi_bubble_slice(j,iy+1)- &
                           ws(igx)*ws(igy)/((w(j)-ek)*(w(j+i)-ekq))
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    chi_bubble_slice=beta*chi_bubble_slice/ &
         ((2.0d0*dfloat(Nint))**2)
    
    DEALLOCATE(w)

  END SUBROUTINE calc_bubble_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=8) FUNCTION find_qmax(mu,beta,Iwbox,self,LQ,Nint,gamma)
    ! input
    INTEGER, INTENT(IN) :: Nint, LQ, Iwbox
    REAL(KIND=8), INTENT(IN) :: mu, beta
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: self
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(IN) :: gamma
    ! subroutine internal variables
    INTEGER :: qsteps,ix,iy,iq
    INTEGER, PARAMETER :: k_number=1
    REAL(KIND=8) :: qmin, qmax, qmax_old, qmin_old
    REAL(KIND=8), DIMENSION(1,1:2) :: klist
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcok,dsik,dcol,dsil
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dcoq,dsiq,Qv
    REAL(KIND=8) :: qy,Q0b
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: chimax_x0
    COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: chi_bubble_slice
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: trilex

    ALLOCATE(chi_bubble_slice(-Iwbox:Iwbox-1,1:LQ))   !for the qmax determination
    ALLOCATE(trilex(-Iwbox:Iwbox-1))
    ALLOCATE(chimax_x0(1:LQ))
    !allocate arrays containing cos and sin functions evaluated for all momenta in the k'-, q- and k-grid
    ALLOCATE(Qv(0:LQ-1))
    !k'-grid
    ALLOCATE(dcok(0:Nint-1,Ng))
    ALLOCATE(dsik(0:Nint-1,Ng))
    !q-grid
    ALLOCATE(dcoq(0:LQ-1))
    ALLOCATE(dsiq(0:LQ-1))
    !k-grid
    ALLOCATE(dcol(k_number,2))
    ALLOCATE(dsil(k_number,2))
    klist=0d0

    qmin=0d0
    qmax=dacos(-1.0d0)     
    qmax_old=qmax
    qmin_old=qmin
    ix=LQ-1
    DO qsteps=1,qmaxsteps
       !initialize the cos()- and sin()-arrays for the three momenta
       CALL init_arrays(Nint,LQ,1,klist,qmin,qmax,qmax,Q0b,Qv,dcok,dsik,dcoq,dsiq,dcol,dsil)
       !calculate bare susceptibility (bubble term)
       CALL calc_bubble_slice(mu,beta,Iwbox,0,self,LQ,Nint,dcok,dsik,dcoq,dsiq,chi_bubble_slice)
       DO iy=0,ix
          qy = Qv(iy)
          
          !calculate chi (without lambda correction)
          chimax_x0(iy+1)=calc_chi(Iwbox,beta,gamma, &
               chi_bubble_slice(:,iy+1))
          
          IF(REAL(chimax_x0(iy+1)) .EQ. 0d0) THEN
             WRITE(6,*)'Warning: ChiS is zero!'
          ENDIF
          
       ENDDO
       iq=MINLOC(1d0/REAL(chimax_x0),1)-1
       IF((iq.EQ.LQ-1 .OR. iq.EQ.0) .AND. qsteps.EQ.1) THEN
          qmax=Qv(LQ-1)
          qmin=Qv(0)
          EXIT
       ELSE IF(iq.EQ.LQ-1) THEN
          qmax=qmax_old
          qmin=Qv(iq-1)
       ELSE IF(iq.EQ.0) THEN
          qmin=qmin_old
          qmax=Qv(iq+1)
       ELSE
          qmin_old=qmin
          qmin=Qv(iq-1)
          qmax_old=qmax
          qmax=Qv(iq+1)
       ENDIF
       WRITE(6,*) qsteps, qmin, qmax
       IF(ABS(qmax_old-qmax).LT.qmaxprec) THEN
          EXIT
       ENDIF
    ENDDO
    find_qmax=qmax
  END FUNCTION find_qmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !computes subsceptibilites (without lambda correction!)
  !input: indices ix, iy for the q-point (Qv(ix),Qv(iy)), lambda_spin (x), lambda_charge(x_ch)
  !output: susceptibilities in chis, chich
  COMPLEX(KIND=8) FUNCTION calc_chi(Iwbox,beta,gamma,chi_bubble)
    IMPLICIT NONE
    !calculates full susceptibility as a function of q and omega
    
    !input:
    INTEGER, INTENT(IN) :: Iwbox
    REAL(KIND=8), INTENT(IN) :: beta
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(IN) :: gamma
    COMPLEX(KIND=8), DIMENSION(-Iwbox:), INTENT(IN) :: chi_bubble
    !function internal variables
    INTEGER :: j,k
    COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: selftemp
    !variable for lapack inversion subroutines
    INTEGER :: infoinv,Mmat
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: workinv
    
    ALLOCATE(selftemp(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
    ALLOCATE(ipiv(1:2*Iwbox))
    ALLOCATE(workinv(1:20*Iwbox))
    Mmat=2*Iwbox
    
    DO k=-Iwbox,Iwbox-1
       DO j=-Iwbox,Iwbox-1
          selftemp(j,k)=-gamma(j,k)
          IF (j.EQ.k) THEN
             selftemp(j,j)=selftemp(j,j)+1.d0/chi_bubble(j)
          ENDIF
       ENDDO
    ENDDO
    
    ! compute inversion chi_q=(chi_bubble_q**-1-Gamma_loc)**-1
    CALL ZGETRF(Mmat,Mmat,selftemp,Mmat,ipiv,infoinv)
    CALL ZGETRI(Mmat,selftemp,Mmat,ipiv,workinv,10*Mmat,infoinv)
    
    ! compute Chi
    calc_chi=(0.d0,0.d0)
    DO k=-iwbox,iwbox-1 
       DO j=-iwbox,iwbox-1
          calc_chi=calc_chi+selftemp(j,k)
       ENDDO
    ENDDO
    
    calc_chi=calc_chi/beta**2

    DEALLOCATE(selftemp)
    DEALLOCATE(ipiv)
    DEALLOCATE(workinv)

  END FUNCTION calc_chi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !computes subsceptibilites (without lambda correction!) and the trilex vertex (for obtaining the selfenergy)
  !input: indices ix, iy for the q-point (Qv(ix),Qv(iy)), lambda_spin (x), lambda_charge(x_ch)
  !output: susceptibilities in chis, chich and trilex vertex gammas, gammach 
  SUBROUTINE calc_chi_trilex(Iwbox,U,beta,gamma,chi_bubble,chi,trilex)
    IMPLICIT NONE
    !calculates full susceptibility for a given q and omega
    !and the trilex vertex for a given q and omega as a function of \nu
    
    !input:
    INTEGER, INTENT(IN) :: Iwbox
    REAL(KIND=8), INTENT(IN) :: U,beta
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(IN) :: gamma
    COMPLEX(KIND=8), DIMENSION(-Iwbox:), INTENT(IN) :: chi_bubble
    !output:
    COMPLEX(KIND=8), INTENT(OUT) :: chi
    COMPLEX(KIND=8), DIMENSION(-Iwbox:), INTENT(OUT) :: trilex
    !function internal variables
    INTEGER :: j,k
    COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: selftemp
    !variable for lapack inversion subroutines
    INTEGER :: infoinv,Mmat
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: workinv
    
    ALLOCATE(selftemp(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
    ALLOCATE(ipiv(1:2*Iwbox))
    ALLOCATE(workinv(1:20*Iwbox))
    Mmat=2*Iwbox
    
    DO k=-Iwbox,Iwbox-1
       DO j=-Iwbox,Iwbox-1
          selftemp(j,k)=-gamma(j,k)
          IF (j.EQ.k) THEN
             selftemp(j,j)=selftemp(j,j)+1.d0/chi_bubble(j)
          ENDIF
       ENDDO
    ENDDO
    
    ! compute inversion chi_q=(chi_bubble_q**-1-Gamma_loc)**-1
    CALL ZGETRF(Mmat,Mmat,selftemp,Mmat,ipiv,infoinv)
    CALL ZGETRI(Mmat,selftemp,Mmat,ipiv,workinv,10*Mmat,infoinv)
    
    ! compute Chi and trilex vertex
    chi=(0.d0,0.d0)
    DO k=-iwbox,iwbox-1
       trilex(k)=(0.d0,0.d0)
       DO j=-iwbox,iwbox-1
          chi=chi+selftemp(j,k)
          trilex(k)=trilex(k)+selftemp(k,j)
       ENDDO
       trilex(k)=trilex(k)/chi_bubble(k)
    ENDDO
    
    chi=chi/beta**2

    trilex=trilex/(1.0d0-U*chi)

    DEALLOCATE(selftemp)
    DEALLOCATE(ipiv)
    DEALLOCATE(workinv)

  END SUBROUTINE calc_chi_trilex
  
END MODULE calc_susc
