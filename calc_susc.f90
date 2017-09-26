!Contains the subroutines for calculating:
!-) non-interacting susceptibility (bare bubble) -> in the programm: chi1(nu)
!-) full interacting susceptibility (in ladder approximation) -> in the program: chisp, und chich
!-) for a given omega and q!
MODULE calc_susc
  USE dispersion
  IMPLICIT NONE

  include '/lrz/sys/libraries/fftw/3.3.3/include/fftw3.f'

  INTEGER, PARAMETER :: qmaxsteps=100
  REAL(KIND=8), PARAMETER :: qmaxprec=1E-14

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
    INTEGER :: j,ind,ix,iy,iz,iNx,iNy,iNz,igx,igy,igz
    REAl(KIND=8) :: pi,ek,ekq
    COMPLEX(KIND=8), ALLOCATABLE :: w(:)

    pi=dacos(-1.0d0)
    ALLOCATE(w(-2*Iwbox:2*Iwbox-1))
    DO j=-2*Iwbox,2*Iwbox-1
       w(j)=dcmplx(0.0d0,(pi/beta)*dfloat(2*j+1))+mu-self(j)
    ENDDO
    
    chi_bubble=dcmplx(0.0d0,0.0d0)

    if(i.eq.0)write(6,*) "----------------"
    if(i.eq.0)write(6,*) "bubble plain"
    DO igx=1, Ng
       DO igy=1, Ng
          DO igz=1,Ng
             DO iNx=0,Nint-1
                if(i.eq.0)write(*,*) "   ", inx, "/", Nint-1
                DO iNy=0,Nint-1
                   DO iNz=0,Nint-1
                      ek=eps(dcok(iNx,igx),dcok(iNy,igy),dcok(iNz,igz), &
                           1.d0,1.d0,1.0d0,0.0d0,0.0d0,0.d0,0.d0,0.d0,0.d0)
                      DO ix=0,LQ-1
                         DO iy=0,ix
                            DO iz=0,iy
                               ind=ix*(ix+1)*(ix+2)/6+iy*(iy+1)/2+iz+1
                               ekq=eps(dcok(iNx,igx),dcok(iNy,igy),dcok(iNz,igz), &
                                    dcoq(ix),dcoq(iy),dcoq(iz), &
                                    dsik(iNx,igx),dsik(iNy,igy),dsik(iNz,igz), &
                                    dsiq(ix),dsiq(iy),dsiq(iz))
                               DO j=0,iwbox-1
                                  chi_bubble(j,ind)=chi_bubble(j,ind)- &
                                       ws(igx)*ws(igy)*ws(igz)/((w(j)-ek)*(w(j+i)-ekq))
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(w)
    chi_bubble=beta*chi_bubble/ &
         ((2.0d0*dfloat(Nint))**3)
  END SUBROUTINE calc_bubble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE calc_bubble_fft(mu,beta,Iwbox,i,self,LQ,Nint,dcok,dsik,dcoq,dsiq, &
       chi_bubble_pos)
    !caluclates the non-interacing susceptibility as a function of (nu,omega,q)
  
    !input:
    INTEGER, INTENT(IN) :: LQ,Nint,i,Iwbox
    REAL(KIND=8), INTENT(IN) :: mu,beta
    REAL(KIND=8), DIMENSION(0:,:), INTENT(IN) :: dcok,dsik
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: dcoq,dsiq
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: self
    !output:
    COMPLEX(KIND=8), DIMENSION(0:,:), INTENT(OUT) :: chi_bubble_pos
    !subroutine internal variables
    INTEGER :: j,ix,iy,iz,iNx,iNy,iNz,igx,igy,igz,ind
    REAl(KIND=8) :: pi,ek,ekq
    COMPLEX(KIND=8), ALLOCATABLE :: w(:)
    COMPLEX(KIND=8), ALLOCATABLE :: Gk1(:,:,:), Gr1(:,:,:)
    COMPLEX(KIND=8), ALLOCATABLE :: Gk2(:,:,:), Gr2(:,:,:)
  
    COMPLEX(KIND=8), ALLOCATABLE :: phasematrix(:,:,:)
    COMPLEX(KIND=8), ALLOCATABLE :: phasevector(:)
  
    integer*8 :: plan
  
    !!! some checks that the parameters are compatible with fft algorithm
    write(*,*) "Ng", Ng
    write(*,*) "Nint", Nint
    write(*,*) "Lq", Lq
    if(Ng.ne.1)then
       write(6,*) "Gauss-Legendre not compatible with FFT bubble!"
       stop
    endif
    if((lq)*2-2.ne.nint)then
       write(6,*) "(LQ)*2-2.ne.Nint: FFT not possible!"
       !write(*,*) "lq", lq
       !write(*,*) "nint", nint
       stop
    endif
  
    pi=dacos(-1.0d0)
    ALLOCATE(w(-2*Iwbox:2*Iwbox-1))
    DO j=-2*Iwbox,2*Iwbox-1
       w(j)=dcmplx(0.0d0,(pi/beta)*dfloat(2*j+1))+mu-self(j)
    ENDDO
  
    !!! initialize Greens functions in real and reciprocal space
    ALLOCATE(Gk1(0:Nint-1,0:Nint-1,0:Nint-1))
    ALLOCATE(Gk2(0:Nint-1,0:Nint-1,0:Nint-1))
    ALLOCATE(Gr1(0:Nint-1,0:Nint-1,0:Nint-1))
    ALLOCATE(Gr2(0:Nint-1,0:Nint-1,0:Nint-1))
    Gk1=0d0
    Gk2=0d0
    Gr1=0d0
    Gr2=0d0
  
    !!! phases necessary for use of fftw package
    ALLOCATE(phasematrix(0:Nint-1,0:Nint-1,0:Nint-1))
    ALLOCATE(phasevector(0:Nint-1))
  
    do inx=0,nint-1
       phasevector(inx)=exp(-2*pi*cmplx(0d0,1d0)*float(inx)/float(nint)) * (-1)**inx
    enddo
  
    do inx=0,nint-1
    do iny=0,nint-1
    do inz=0,nint-1
       phasematrix(inx,iny,inz)=phasevector(inx)*phasevector(iny)*phasevector(inz)
    enddo
    enddo
    enddo
    
    chi_bubble_pos=dcmplx(0.0d0,0.0d0)
  
    !!! bubble with fft
    if(i.eq.0)write(6,*) "----------------"
    if(i.eq.0)write(6,*) "bubble with fftw"
    DO j=0,iwbox-1
       if(i.eq.0)write(6,*) "iw / iwbox", j, "/", iwbox
       DO igx=1, Ng
          DO igy=1, Ng
             DO igz=1, Ng

                DO iNx=0,Nint-1
                   DO iNy=0,Nint-1
                      DO iNz=0,Nint-1
                         ek=eps(dcok(iNx,igx),dcok(iNy,igy),dcok(iNz,igz), &
                            1.d0,1.d0,1.0d0,0.0d0,0.0d0,0.d0,0.d0,0.d0,0.d0)
                         Gk1(inx,iny,inz)=ws(igx)/(w(j)-ek)
                         Gk2(inx,iny,inz)=ws(igy)/(w(j+i)-ek)
                      ENDDO
                   ENDDO
                ENDDO
  
                if(i.eq.0)then
                DO iNx=0,Nint-1
               !DO iNy=0,Nint-1
               !DO iNz=0,Nint-1
               iny=1
               inz=1
                write(100,*) gk1(inx,iny,inz)
                enddo
               !enddo
               !enddo
                endif
  
                !!! forward Fourier transform
                call dfftw_plan_dft_3d(plan,Nint,Nint,Nint,Gk1,Gr1,FFTW_FORWARD,FFTW_ESTIMATE)
                call dfftw_execute(plan)
                call dfftw_destroy_plan(plan)
                !!! forward Fourier transform
                call dfftw_plan_dft_3d(plan,Nint,Nint,Nint,Gk2,Gr2,FFTW_FORWARD,FFTW_ESTIMATE)
                call dfftw_execute(plan)
                call dfftw_destroy_plan(plan)

                 if(i.eq.0)then
                 DO iNx=0,Nint-1
                !DO iNy=0,Nint-1
                !DO iNz=0,Nint-1
                iny=1
                inz=1
                 write(101,*) gr1(inx,iny,inz)
                 enddo
                !enddo
                !enddo
                 endif
  
                !!! normalize
                !Gr1=Gr1/sqrt(float(nint**3))
                !Gr2=Gr2/sqrt(float(nint**3))
  
                !!! apply phase factor (TODO documentation of this!)
                DO iNx=0,Nint-1
                DO iNy=0,Nint-1
                DO iNz=0,Nint-1
                   Gr1(inx,iny,inz)=Gr1(inx,iny,inz)*phasematrix(inx,iny,inz)
                   Gr2(inx,iny,inz)=Gr2(inx,iny,inz)*phasematrix(inx,iny,inz)
                ENDDO
                ENDDO
                ENDDO
  
                !!! convolute in real space
                DO iNx=0,Nint-1
                DO iNy=0,Nint-1
                DO iNz=0,Nint-1
                   Gr1(inx,iny,inz)=Gr1(inx,iny,inz)*Gr2(inx,iny,inz)
                ENDDO
                ENDDO
                ENDDO
  
                !!! backward Fourier transform
                Gk1=0d0
                call dfftw_plan_dft_3d(plan,Nint,Nint,Nint,Gr1,Gk1,FFTW_BACKWARD,FFTW_ESTIMATE)
                call dfftw_execute(plan)
                call dfftw_destroy_plan(plan)

                !!! normalize
                Gk1=Gk1/(float(nint**(9/2)))

                if(i.eq.0)then
                DO iNx=0,Nint-1
               !DO iNy=0,Nint-1
               !DO iNz=0,Nint-1
               iny=1
               inz=1
                write(102,*) gk1(inx,iny,inz)
                enddo
               !enddo
               !enddo
                endif
  
                !!! write values of full-BZ bubble in reduced-BZ bubble
                do ix=0,LQ-1
                do iy=0,ix
                do iz=0,iy
                   !ind=ix*(ix+1)/2+iy+1   !!! index for fully irred. BZ
                   ind=ix*(ix+1)*(ix+2)/6+iy*(iy+1)/2+iz+1
                   chi_bubble_pos(j,ind)=chi_bubble_pos(j,ind)-Gk1(ix,iy,iz)
                enddo
                enddo
                enddo
  
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  
    chi_bubble_pos=beta*chi_bubble_pos/ &
         ((2.0d0**2*dfloat(Nint**2)))
  
    DEALLOCATE(w,Gk1,Gr1,Gk2,Gr2)
  
  END SUBROUTINE calc_bubble_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE calc_bubble_slice(mu,beta,Iwbox,i,self,LQ,Nint,dcok,dsik,dcoq,dsiq, &
       chi_bubble_slice)
    !caluclates the non-interacing susceptibility as a function of (nu,omega,q=(pi,pi,qz))

    !input:
    INTEGER, INTENT(IN) :: LQ,Nint,i,Iwbox
    REAL(KIND=8), INTENT(IN) :: mu,beta
    REAL(KIND=8), DIMENSION(0:,:), INTENT(IN) :: dcok,dsik
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: dcoq,dsiq
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: self
    !output:
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,:), INTENT(OUT) :: chi_bubble_slice
    !subroutine internal variables
    INTEGER :: j,ind,ix,iy,iz,iNx,iNy,iNz,igx,igy,igz
    REAl(KIND=8) :: pi,ek,ekq
    COMPLEX(KIND=8), ALLOCATABLE :: w(:)

    pi=dacos(-1.0d0)
    ALLOCATE(w(-2*Iwbox:2*Iwbox-1))
    DO j=-2*Iwbox,2*Iwbox-1
       w(j)=dcmplx(0.0d0,(pi/beta)*dfloat(2*j+1))+mu-self(j)
    ENDDO

    chi_bubble_slice=dcmplx(0.0d0,0.0d0)
    ix=LQ-1
    iy=LQ-1
    DO igx=1, Ng
       DO igy=1, Ng
          DO igz=1,Ng
             DO iNx=0,Nint-1
                DO iNy=0,Nint-1
                   DO iNz=0,Nint-1
                      ek=eps(dcok(iNx,igx),dcok(iNy,igy),dcok(iNz,igz), &
                           1.d0,1.d0,1.0d0,0.0d0,0.0d0,0.d0,0.d0,0.d0,0.d0)
                      DO iz=0,iy
                         ekq=eps(dcok(iNx,igx),dcok(iNy,igy),dcok(iNz,igz), &
                              -1.0d0,-1.0d0,dcoq(iz), &
                              dsik(iNx,igx),dsik(iNy,igy),dsik(iNz,igz), &
                              0.0d0,0.0d0,dsiq(iz))
                         DO j=-iwbox,iwbox-1
                            chi_bubble_slice(j,iz+1)=chi_bubble_slice(j,iz+1)- &
                                 ws(igx)*ws(igy)*ws(igz)/((w(j)-ek)*(w(j+i)-ekq))
                            
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(w)

    chi_bubble_slice=beta*chi_bubble_slice/ &
         ((2.0d0*dfloat(Nint))**3)

  END SUBROUTINE calc_bubble_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(KIND=8) FUNCTION find_qmax(mu,beta,Iwbox,self,LQ,Nint,gamma)
    ! input
    INTEGER, INTENT(IN) :: Nint, LQ, Iwbox
    REAL(KIND=8), INTENT(IN) :: mu, beta
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: self
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(IN) :: gamma
    ! subroutine internal variables
    INTEGER :: qsteps,ix,iy,iz,iq
    INTEGER, PARAMETER :: k_number=1
    REAL(KIND=8) :: qmin, qmax, qmax_old, qmin_old
    REAL(KIND=8), DIMENSION(1,1:3) :: klist
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dcok,dsik,dcol,dsil
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dcoq,dsiq,Qv
    REAL(KIND=8) :: qz,Q0b
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
    ALLOCATE(dcol(k_number,3))
    ALLOCATE(dsil(k_number,3))
    klist=0d0

    qmin=0d0
    qmax=dacos(-1.0d0)     
    qmax_old=qmax
    qmin_old=qmin
    ix=LQ-1
    iy=LQ-1
    DO qsteps=1,qmaxsteps
       !initialize the cos()- and sin()-arrays for the three momenta
       CALL init_arrays(Nint,LQ,1,klist,qmin,qmax,qmax,Q0b,Qv,dcok,dsik,dcoq,dsiq,dcol,dsil)
       !calculate bare susceptibility (bubble term)
       CALL calc_bubble_slice(mu,beta,Iwbox,0,self,LQ,Nint,dcok,dsik,dcoq,dsiq,chi_bubble_slice)
       DO iz=0,iy
          qz = Qv(iz)
          
          !calculate chi (without lambda correction)
          chimax_x0(iz+1)=calc_chi(Iwbox,beta,gamma, &
               chi_bubble_slice(:,iz+1))
          
          IF(REAL(chimax_x0(iz+1)) .EQ. 0d0) THEN
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
       WRITE(6,*)qsteps,qmin,qmax
       IF(ABS(qmax_old-qmax).LT.qmaxprec) THEN
          EXIT
       ENDIF
    ENDDO
    find_qmax=qmax
  END FUNCTION find_qmax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
