!Contains the subroutines for calculating:
!-) non-interacting susceptibility (bare bubble) -> in the programm: chi_bubble(nu)
!-) full interacting susceptibility (in ladder approximation) -> in the program: chisp, und chich
!-) for a given omega and q!
MODULE calc_susc
  USE dispersion
  IMPLICIT NONE

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
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,:), INTENT(OUT) :: chi_bubble
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
                      DO j=-iwbox,iwbox-1
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
  
END MODULE calc_susc
