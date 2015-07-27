!     This module calculates local and nonlocal self-energy with the DGA-equation

MODULE sigma
   USE dispersion
   IMPLICIT NONE

CONTAINS
  
  SUBROUTINE calc_self_loc(Iwbox,i,uhub,beta,gww,fupdown,gammach,gammasp, &
       self_loc)
    !this subroutine calculates the local self-energy with the DGA equation for one bosonic frequency
    !should be extended to more bosonic frequencies, if there is not on core per bosonic frequency available!!!
    
    !input: # fermionix frequencies, bosonic frequency, interaction strength, inverse temperature, 
    !DMFT-self-energy,DMFT-green's function,DMFT-irreducible vertex
    INTEGER, INTENT(IN) :: Iwbox,i
    REAL(KIND=8), INTENT(IN) :: uhub,beta
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: gww
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(IN) :: fupdown
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(IN) :: gammach,gammasp
    !output: chi_ch_loc(omega), chi_sp_loc(omega), local self-energy calculated with DGA equation
    COMPLEX(KIND=8), DIMENSION(0:), INTENT(OUT) :: self_loc
!    COMPLEX(KIND=8), INTENT(OUT) :: chich_loc,chisp_loc
    !subroutine internal variables
    INTEGER :: j,k
    COMPLEX(KIND=8) :: chich_loc,chisp_loc
    COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: selftempch,selftempsp
    COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: selfphich,selfphisp
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: vrgch,vrgsp
    !variable for lapack inversion subroutines
    INTEGER :: infoinv,Mmat
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: workinv

    ALLOCATE(selftempch(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
    ALLOCATE(selftempsp(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
    ALLOCATE(selfphich(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
    ALLOCATE(selfphisp(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
    ALLOCATE(vrgch(0:Iwbox-1))
    ALLOCATE(vrgsp(0:Iwbox-1))
    ALLOCATE(ipiv(1:2*Iwbox))
    ALLOCATE(workinv(1:20*Iwbox))
    Mmat=2*Iwbox

    DO j=-Iwbox,Iwbox-1
       DO k=-Iwbox,Iwbox-1
          selftempch(j,k)=-gammach(j,k)
          selftempsp(j,k)=-gammasp(j,k)
          selfphich(j,k)=-gammach(j,k)-uhub/beta**2
          selfphisp(j,k)=-gammasp(j,k)+uhub/beta**2
          IF(j .eq. k) THEN
             selftempch(j,j)= selftempch(j,j)+ & 
                  1.d0/(-beta*gww(j)*gww(j+i))
             selftempsp(j,j)= selftempsp(j,j)+ &
                  1.d0/(-beta*gww(j)*gww(j+i))
             selfphich(j,j)= selfphich(j,j)+ &
                  1.d0/(-beta*gww(j)*gww(j+i))
             selfphisp(j,j)= selfphisp(j,j)+ &
                  1.d0/(-beta*gww(j)*gww(j+i))
          ENDIF
       ENDDO
    ENDDO
    
    ! Compute inversion chi_q=(chi_0_loc**-1-Gamma_loc)**-1
    CALL ZGETRF(Mmat,Mmat,selftempch,Mmat,ipiv,infoinv)
    CALL ZGETRI(Mmat,selftempch,Mmat,ipiv,workinv,10*Mmat,infoinv)
    
    CALL ZGETRF(Mmat,Mmat,selftempsp,Mmat,ipiv,infoinv) 
    CALL ZGETRI(Mmat,selftempsp,Mmat,ipiv,workinv,10*Mmat,infoinv)
    
    CALL ZGETRF(Mmat,Mmat,selfphich,Mmat,ipiv,infoinv)
    CALL ZGETRI(Mmat,selfphich,Mmat,ipiv,workinv,10*Mmat,infoinv)
    
    CALL ZGETRF(Mmat,Mmat,selfphisp,Mmat,ipiv,infoinv)
    CALL ZGETRI(Mmat,selfphisp,Mmat,ipiv,workinv,10*Mmat,infoinv)
    
    DO j=0,Iwbox-1 
       vrgch(j)=dcmplx(0.d0,0.d0)
       vrgsp(j)=dcmplx(0.d0,0.d0)
       DO k=-Iwbox,Iwbox-1
          vrgch(j)=vrgch(j)+selfphich(k,j)
          vrgsp(j)=vrgsp(j)+selfphisp(k,j) 
       ENDDO
       vrgch(j)=vrgch(j)/(-beta*gww(j)*gww(j+i))
       vrgsp(j)=vrgsp(j)/(-beta*gww(j)*gww(j+i))
    ENDDO

    chich_loc=dcmplx(0.d0,0.d0)
    chisp_loc=dcmplx(0.d0,0.d0) 
    DO j=-Iwbox,Iwbox-1 
       DO k=-Iwbox,Iwbox-1
          chich_loc=chich_loc+selftempch(k,j)
          chisp_loc=chisp_loc+selftempsp(k,j)
       ENDDO
    ENDDO
    
    chich_loc=chich_loc/beta**2
    chisp_loc=chisp_loc/beta**2
    
    DO j=0,Iwbox-1
       self_loc(j)= &
            (1.5d0*vrgsp(j)*(1.d0+uhub*chisp_loc)- &
            0.5d0*vrgch(j)*(1.d0-uhub*chich_loc)- &
            1.0d0)*gww(i+j)
    ENDDO
    
    DO j=0,Iwbox-1
       DO k=-Iwbox,Iwbox-1
          self_loc(j)=self_loc(j)-gww(k)*gww(k+i)* &
               fupdown(j,k)*gww(i+j)*beta
       ENDDO
    ENDDO

    DEALLOCATE(selftempch)
    DEALLOCATE(selftempsp)
    DEALLOCATE(selfphich)
    DEALLOCATE(selfphisp)
    DEALLOCATE(vrgch)
    DEALLOCATE(vrgsp)
    DEALLOCATE(ipiv)
    DEALLOCATE(workinv)
    
  END SUBROUTINE calc_self_loc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE calc_self(Iwbox,myid,i,LQ,k_number,uhub,mu,beta,lambdach,lambdasp, &
       self,fupdown,gammach,gammasp, chi_bubble,chich_x0,chisp_x0, &
       dcoq,dsiq,dcol,dsil, &
       selflist,selflistch,selflistsp,selflistrest)
    !this subroutine calculates the DGA self-energy  for one bosonic frequency
    !should be extended to more bosonic frequencies, if there is not on core per bosonic frequency available!!!
    
    !input: # fermionix frequencies, bosonic frequency, interaction strength, inverse temperature, 
    !DMFT-self-energy,DMFT-green's function,DMFT-irreducible vertex
    INTEGER, INTENT(IN) :: Iwbox,myid,i,k_number,LQ
    REAL(KIND=8), INTENT(IN) :: uhub,mu,beta,lambdach,lambdasp
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: dcoq,dsiq
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: dcol,dsil
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: self
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(IN) :: fupdown
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(IN) :: gammach,gammasp
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,:), INTENT(IN) :: chi_bubble
    COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: chich_x0,chisp_x0
    !output: chi_ch_loc(omega), chi_sp_loc(omega), local self-energy calculated with DGA equation
    COMPLEX(KIND=8), DIMENSION(:,0:), INTENT(OUT) :: selflist,selflistch
    COMPLEX(KIND=8), DIMENSION(:,0:), INTENT(OUT) :: selflistsp,selflistrest
    !subroutine internal variables
    INTEGER :: j,k,ix,iy,s,t,i1,perm,ind
    REAL(KIND=8) :: a,b,chich,chisp,pi,Q0b
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: eklist
    COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: selfphich,selfphisp
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: vrgch,vrgsp
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: w,chi0,chi0ch,chi0sp,chi0rest
    !variable for lapack inversion subroutines
    INTEGER :: infoinv,Mmat
    INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: workinv

    ALLOCATE(w(-2*Iwbox:2*Iwbox-1))
    ALLOCATE(selfphich(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
    ALLOCATE(selfphisp(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
    ALLOCATE(vrgch(0:Iwbox-1))
    ALLOCATE(vrgsp(0:Iwbox-1))
    ALLOCATE(ipiv(1:2*Iwbox))
    ALLOCATE(workinv(1:20*Iwbox))
    Mmat=2*Iwbox
    ALLOCATE(chi0(0:Iwbox-1))
    ALLOCATE(chi0ch(0:Iwbox-1))
    ALLOCATE(chi0sp(0:Iwbox-1))
    ALLOCATE(chi0rest(0:Iwbox-1))
    ALLOCATE(eklist(1:k_number,4,2))

    pi=dacos(-1.0d0)
    Q0b=pi/dfloat(LQ-1)

    !Initialize w
    DO j=-2*Iwbox,2*Iwbox-1
       w(j)=dcmplx(0.0d0,(pi/beta)*dfloat(2*j+1))+mu-self(j)
    ENDDO

    !Initalize selflists
    DO j=0,Iwbox-1
       DO i1=1,k_number
          selflist(i1,j)=dcmplx(0.0d0,0.0d0)
          selflistch(i1,j)=dcmplx(0.0d0,0.0d0)
          selflistsp(i1,j)=dcmplx(0.0d0,0.0d0)
          selflistrest(i1,j)=dcmplx(0.0d0,0.0d0)
       ENDDO
    ENDDO

    DO ix=0,LQ-1
       a=1.d0
       IF(ix.EQ.0.OR.ix.EQ.LQ-1) THEN
          a=a*0.5d0
       ENDIF
       DO iy=0,ix
          b=a
          IF(iy.EQ.0.OR.iy.EQ.LQ-1) THEN
             b=b*0.5d0
          ENDIF
            
          b=b*dfloat(2/((1+iy)/(1+ix)+1))/2.0d0
          
          ind=ix*(ix+1)/2+iy+1

          IF (myid.EQ.0) THEN
             WRITE(6,*) 'i,j=',ix,iy
          ENDIF
          
          !compute Lambda-corrected Chi_lambda
          chich=1.d0/(1.d0/chich_x0(ind)+lambdach)
          chisp=1.d0/(1.d0/chisp_x0(ind)+lambdasp)
          
          !symmetry conditions for the fully irreducible BZ
          DO s=0,1
             DO t=0,1
                DO i1=1,k_number
                   eklist(i1,2*s+t+1,1)=eps( &
                        dcol(i1,1),dcol(i1,2), &
                        dcoq(ix),dcoq(iy), &
                        dsil(i1,1),dsil(i1,2), &
                        (-1)**s*dsiq(ix),(-1)**t*dsiq(iy))
                   eklist(i1,2*s+t+1,2)=eps( &
                        dcol(i1,1),dcol(i1,2), &
                        dcoq(iy),dcoq(ix), &
                        dsil(i1,1),dsil(i1,2), &
                        (-1)**t*dsiq(iy),(-1)**s*dsiq(ix))
                ENDDO
             ENDDO
          ENDDO
          
          DO k=-Iwbox,Iwbox-1
             DO j=-Iwbox,Iwbox-1
                selfphich(j,k)=-gammach(j,k)-uhub/beta**2
                selfphisp(j,k)=-gammasp(j,k)+uhub/beta**2
                IF(j .EQ. k) THEN
                   selfphich(j,j)=selfphich(j,j)+ &
                        1.d0/chi_bubble(j,ind)
                   selfphisp(j,j)=selfphisp(j,j)+ &
                        1.d0/chi_bubble(j,ind)
                ENDIF
             ENDDO
          ENDDO
          
          !Compute inversion chi_q=(chi_0_loc**-1-Gamma_loc)**-1
          CALL ZGETRF(Mmat,Mmat,selfphich,Mmat,ipiv,infoinv)
          CALL ZGETRI(Mmat,selfphich,Mmat,ipiv,workinv,10*Mmat,infoinv)
          
          CALL ZGETRF(Mmat,Mmat,selfphisp,Mmat,ipiv,infoinv)
          CALL ZGETRI(Mmat,selfphisp,Mmat,ipiv,workinv,10*Mmat,infoinv)
          
          DO j=0,Iwbox-1 
             vrgch(j)=dcmplx(0.d0,0.d0)
             vrgsp(j)=dcmplx(0.d0,0.d0)
             DO k=-Iwbox,Iwbox-1
                vrgch(j)=vrgch(j)+selfphich(k,j)
                vrgsp(j)=vrgsp(j)+selfphisp(k,j) 
             ENDDO
             vrgch(j)=vrgch(j)/chi_bubble(j,ind)
             vrgsp(j)=vrgsp(j)/chi_bubble(j,ind)
          ENDDO
          
          !compute Sigma_q
          DO j=0,Iwbox-1
             chi0(j)=dcmplx(0.0d0,0.0d0)
             chi0ch(j)=dcmplx(0.0d0,0.0d0)
             chi0sp(j)=dcmplx(0.0d0,0.0d0)
             chi0rest(j)=dcmplx(0.0d0,0.0d0)
             DO k=-Iwbox,Iwbox-1
                chi0(j)=chi0(j)+chi_bubble(k,ind)*fupdown(j,k)
                chi0ch(j)=chi0ch(j)+0.5d0*chi_bubble(k,ind)*gammach(j,k)
                chi0sp(j)=chi0sp(j)-1.5d0*chi_bubble(k,ind)*gammasp(j,k)
                chi0rest(j)=chi0rest(j)+chi_bubble(k,ind)* &
                     (1.5d0*gammasp(j,k)-0.5d0*gammach(j,k)+fupdown(j,k))
             ENDDO
          ENDDO
          DO j=0,Iwbox-1
             chi0(j)=(chi0(j)+1.5d0*vrgsp(j)* &
                  (1.d0+uhub*chisp)-0.5d0*vrgch(j)* &
                  (1.d0-uhub*chich)-1.5d0+0.5d0)*uhub/beta
             chi0ch(j)=(chi0ch(j)-0.5d0*vrgch(j)* &
                  (1.d0-uhub*chich)+0.5d0)*uhub/beta
             chi0sp(j)=(chi0sp(j)+1.5d0*vrgsp(j)* &
                  (1.d0+uhub*chisp)-1.5d0)*uhub/beta
             chi0rest(j)=(chi0rest(j))*uhub/beta
          ENDDO
          
          DO j=0,Iwbox-1
             DO s=1,4
                DO i1=1,k_number
                   DO perm=1,2 ! permutations of fully irreducible BZ
                      selflist(i1,j)=selflist(i1,j)+chi0(j)/ &
                           (w(i+j)-Eklist(i1,s,perm))* &
                           Q0b**2/(4.0d0*pi**2)*b
                      selflistch(i1,j)=selflistch(i1,j)+chi0ch(j)/ &
                           (w(i+j)-Eklist(i1,s,perm))* &
                           Q0b**2/(4.0d0*pi**2)*b
                      selflistsp(i1,j)=selflistsp(i1,j)+chi0sp(j)/ &
                           (w(i+j)-Eklist(i1,s,perm))* &
                           Q0b**2/(4.0d0*pi**2)*b
                      selflistrest(i1,j)=selflistrest(i1,j)+chi0rest(j)/ &
                           (w(i+j)-Eklist(i1,s,perm))* &
                           Q0b**2/(4.0d0*pi**2)*b
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(w)
    DEALLOCATE(selfphich)
    DEALLOCATE(selfphisp)
    DEALLOCATE(vrgch)
    DEALLOCATE(vrgsp)
    DEALLOCATE(ipiv)
    DEALLOCATE(workinv)
    DEALLOCATE(chi0)
    DEALLOCATE(chi0ch)
    DEALLOCATE(chi0sp)
    DEALLOCATE(chi0rest)
    DEALLOCATE(eklist)
    
  END SUBROUTINE calc_self
  
END MODULE sigma
