!Contains the subroutines for calculating:
!-) non-interacting susceptibility (bare bubble) -> in the programm: chi_bubble(nu)
!-) full interacting susceptibility (in ladder approximation) -> in the program: chisp, und chich
!-) for a given omega and q!
MODULE calc_susc
  USE dispersion
  IMPLICIT NONE

  include '/lrz/sys/libraries/fftw/3.3.3/include/fftw3.f'

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

   if(i.eq.0)write(6,*) "----------------"
   if(i.eq.0)write(6,*) "bubble plain"
    DO igx=1, Ng
       DO igy=1, Ng
          DO iNx=0,Nint-1
             if(i.eq.0)write(6,*) "inx", inx, "/", Nint
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
    INTEGER :: j,ix,iy,iNx,iNy,igx,igy,ind
    REAl(KIND=8) :: pi,ek,ekq
    COMPLEX(KIND=8), ALLOCATABLE :: w(:)
    COMPLEX(KIND=8), ALLOCATABLE :: Gk1(:,:), Gr1(:,:)
    COMPLEX(KIND=8), ALLOCATABLE :: Gk2(:,:), Gr2(:,:)
  
    COMPLEX(KIND=8), ALLOCATABLE :: phasematrix(:,:)
    COMPLEX(KIND=8), ALLOCATABLE :: phasevector(:)
  
    integer*8 :: plan
  
    !!! some checks that the parameters are compatible with fft algorithm
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
    ALLOCATE(Gk1(0:Nint-1,0:Nint-1))
    ALLOCATE(Gk2(0:Nint-1,0:Nint-1))
    ALLOCATE(Gr1(0:Nint-1,0:Nint-1))
    ALLOCATE(Gr2(0:Nint-1,0:Nint-1))
    Gk1=0d0
    Gk2=0d0
    Gr1=0d0
    Gr2=0d0
  
    !!! phases necessary for use of fftw package
    ALLOCATE(phasematrix(0:Nint-1,0:Nint-1))
    ALLOCATE(phasevector(0:Nint-1))
  
    do inx=0,nint-1
       phasevector(inx)=exp(-2*pi*cmplx(0d0,1d0)*float(inx)/float(nint)) * (-1)**inx
    enddo
  
    do inx=0,nint-1
    do iny=0,nint-1
       phasematrix(inx,iny)=phasevector(inx)*phasevector(iny)
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
          DO iNx=0,Nint-1
                DO iNy=0,Nint-1
                   ek=eps(dcok(iNx,igx),dcok(iNy,igy), &
                        1.d0,1.d0,0.d0,0.d0,0.d0,0.d0)
                   Gk1(inx,iny)=ws(igx)/(w(j)-ek)
                   Gk2(inx,iny)=ws(igy)/(w(j+i)-ek)
                ENDDO
             ENDDO
  
            !if(i.eq.0)then
            !DO iNx=0,Nint-1
            !write(*,*) "gk1(inx,1)", gk1(inx,1)
            !enddo
            !endif
  
            !!! forward Fourier transform
            call dfftw_plan_dft_2d(plan,Nint,Nint,Gk1,Gr1,FFTW_FORWARD,FFTW_ESTIMATE)
            call dfftw_execute(plan)
            call dfftw_destroy_plan(plan)
            !!! forward Fourier transform
            call dfftw_plan_dft_2d(plan,Nint,Nint,Gk2,Gr2,FFTW_FORWARD,FFTW_ESTIMATE)
            call dfftw_execute(plan)
            call dfftw_destroy_plan(plan)
  
            !!! normalize
            Gr1=Gr1/sqrt(float(nint**2))
            Gr2=Gr2/sqrt(float(nint**2))
  
            !!! apply phase factor (TODO documentation of this!)
            DO iNx=0,Nint-1
            DO iNy=0,Nint-1
               Gr1(inx,iny)=Gr1(inx,iny)*phasematrix(inx,iny)
               Gr2(inx,iny)=Gr2(inx,iny)*phasematrix(inx,iny)
            ENDDO
            ENDDO
  
            !!! convolute in real space
            DO iNx=0,Nint-1
            DO iNy=0,Nint-1
               Gr1(inx,iny)=Gr1(inx,iny)*Gr2(inx,iny)
            ENDDO
            ENDDO
  
            !!! backward Fourier transform
            Gk1=0d0
            call dfftw_plan_dft_2d(plan,Nint,Nint,Gr1,Gk1,FFTW_BACKWARD,FFTW_ESTIMATE)
            call dfftw_execute(plan)
            call dfftw_destroy_plan(plan)
  
            !!! normalize
            Gk1=Gk1/sqrt(float(nint**2))
  
            !!! write values of full-BZ bubble in reduced-BZ bubble
            do ix=0,LQ-1
            do iy=0,ix
               ind=ix*(ix+1)/2+iy+1   !!! index for fully irred. BZ
               chi_bubble_pos(j,ind)=chi_bubble_pos(j,ind)-Gk1(ix,iy)
            enddo
            enddo
  
          ENDDO
       ENDDO
    ENDDO
  
    chi_bubble_pos=beta*chi_bubble_pos/ &
         ((2.0d0**2*dfloat(Nint)))
  
    DEALLOCATE(w,Gk1,Gr1,Gk2,Gr2)
  
  END SUBROUTINE calc_bubble_fft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE calc_bubble_fft_real(mu,beta,Iwbox,i,self,LQ,Nint,dcok,dsik,dcoq,dsiq, &
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
    COMPLEX(KIND=8), ALLOCATABLE :: Gk1(:,:), Gr1(:,:)
    COMPLEX(KIND=8), ALLOCATABLE :: Gk2(:,:), Gr2(:,:)
    COMPLEX(KIND=8), ALLOCATABLE :: Gk_tmp(:,:)
    REAL(KIND=8), ALLOCATABLE :: Gr1_real(:,:), Gk1_real(:,:)
    REAL(KIND=8), ALLOCATABLE :: Gr1_imag(:,:), Gk1_imag(:,:)
    REAL(KIND=8), ALLOCATABLE :: Gr2_real(:,:), Gk2_real(:,:)
    REAL(KIND=8), ALLOCATABLE :: Gr2_imag(:,:), Gk2_imag(:,:)
    integer :: shiftarray(Nint)

    integer*8 :: plan

    !!! some checks that the parameters are compatible with fft algorithm
    if(Ng.ne.1)then
       write(6,*) "Gauss-Legendre not compatible with FFT bubble!"
       stop
    endif
    if((lq)*2-2.ne.nint)then
       write(6,*) "(LQ)*2-2.ne.Nint: FFT not possible!"
       write(6,*) "lq=", lq
       write(6,*) "nint=", nint
       stop
    endif

    pi=dacos(-1.0d0)
    ALLOCATE(w(-2*Iwbox:2*Iwbox-1))
    DO j=-2*Iwbox,2*Iwbox-1
       w(j)=dcmplx(0.0d0,(pi/beta)*dfloat(2*j+1))+mu-self(j)
    ENDDO

    !!! allocate Greens functions in real and reciprocal space
    ALLOCATE(Gk1_real(0:Nint/2,0:Nint/2))
    ALLOCATE(Gk1_imag(0:Nint/2,0:Nint/2))
    ALLOCATE(Gk2_real(0:Nint/2,0:Nint/2))
    ALLOCATE(Gk2_imag(0:Nint/2,0:Nint/2))

    ALLOCATE(Gr1_real(0:Nint/2,0:Nint/2))
    ALLOCATE(Gr1_imag(0:Nint/2,0:Nint/2))
    ALLOCATE(Gr2_real(0:Nint/2,0:Nint/2))
    ALLOCATE(Gr2_imag(0:Nint/2,0:Nint/2))

    ALLOCATE(Gk1(0:Nint-1,0:Nint-1))
    ALLOCATE(Gk2(0:Nint-1,0:Nint-1))
    ALLOCATE(Gr1(0:Nint-1,0:Nint-1))
    ALLOCATE(Gr2(0:Nint-1,0:Nint-1))

    ALLOCATE(Gk_tmp(0:Nint-1,0:Nint-1))

    Gk1=0d0
    Gk2=0d0
    Gr1=0d0
    Gr2=0d0

    Gk1_real=0d0
    Gk2_real=0d0
    Gk1_imag=0d0
    Gk2_imag=0d0

    Gr1_real=0d0
    Gr2_real=0d0
    Gr1_imag=0d0
    Gr2_imag=0d0

    Gk_tmp=0d0

    chi_bubble=dcmplx(0.0d0,0.0d0)

    !!! bubble with fft
    if(i.eq.0)write(6,*) "----------------"
    if(i.eq.0)write(6,*) "bubble with fftw"
    DO j=0,iwbox-1
       if(i.eq.0)write(6,*) "iw / iwbox", j, "/", iwbox
       DO igx=1, Ng
          DO igy=1, Ng
          
            !!! generate Greens functions
            DO iNx=0,Nint-1
               DO iNy=0,Nint-1
                  ek=eps(dcok(iNx,igx),dcok(iNy,igy), &
                       1.d0,1.d0,0.d0,0.d0,0.d0,0.d0)
                  Gk1(inx,iny)=ws(igx)/(w(j)-ek)
                  Gk2(inx,iny)=ws(igy)/(w(j+i)-ek)
               ENDDO
            ENDDO

            !!! shift value without parter to beginning (for fftw)
            shiftarray=-1
            Gk_tmp=cshift(Gk1,shift=shiftarray,dim=1)
            Gk_tmp=cshift(Gk_tmp,shift=shiftarray,dim=2)
            Gk1=Gk_tmp
            Gk_tmp=cshift(Gk2,shift=shiftarray,dim=1)
            Gk_tmp=cshift(Gk_tmp,shift=shiftarray,dim=2)
            Gk2=Gk_tmp

            !!! pick one quadrant
            Gk1_real=real(Gk1(0:nint/2,0:nint/2))
            Gk1_imag=imag(Gk1(0:nint/2,0:nint/2))

            Gk2_real=real(Gk2(0:nint/2,0:nint/2))
            Gk2_imag=imag(Gk2(0:nint/2,0:nint/2))

            !!! Fourier transform
            call dfftw_plan_r2r_2d(plan,Nint/2+1,Nint/2+1,Gk1_real,Gr1_real,FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE)
            call dfftw_execute(plan)
            call dfftw_destroy_plan(plan)

            call dfftw_plan_r2r_2d(plan,Nint/2+1,Nint/2+1,Gk1_imag,Gk1_imag,FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE)
            call dfftw_execute(plan)
            call dfftw_destroy_plan(plan)

            !!! 
            call dfftw_plan_r2r_2d(plan,Nint/2+1,Nint/2+1,Gk2_real,Gk2_real,FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE)
            call dfftw_execute(plan)
            call dfftw_destroy_plan(plan)

            call dfftw_plan_r2r_2d(plan,Nint/2+1,Nint/2+1,Gk2_imag,Gk2_imag,FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE)
            call dfftw_execute(plan)
            call dfftw_destroy_plan(plan)

            !!! normalize
            Gr1_real=Gr1_real/float(nint)
            Gk1_imag=Gk1_imag/float(nint)
            Gk2_real=Gk2_real/float(nint)
            Gk2_imag=Gk2_imag/float(nint)

            !!! apply phase factor
            DO iNx=0,Nint/2
            DO iNy=0,Nint/2
               Gr1_real(inx,iny)=Gr1_real(inx,iny)*(-1)**(inx+iny)
               Gk1_imag(inx,iny)=Gk1_imag(inx,iny)*(-1)**(inx+iny)
               Gk2_real(inx,iny)=Gk2_real(inx,iny)*(-1)**(inx+iny)
               Gk2_imag(inx,iny)=Gk2_imag(inx,iny)*(-1)**(inx+iny)
            ENDDO
            ENDDO

            !!! convolute in real space
            DO iNx=0,Nint/2
            DO iNy=0,Nint/2
               Gr1_real(inx,iny)=Gr1_real(inx,iny)*Gk2_real(inx,iny)-Gk1_imag(inx,iny)*Gk2_imag(inx,iny)
            ENDDO
            ENDDO

            !!! Fourier transform back
            call dfftw_plan_r2r_2d(plan,Nint/2+1,Nint/2+1,Gr1_real,Gk1_real,FFTW_REDFT00,FFTW_REDFT00,FFTW_ESTIMATE)
            call dfftw_execute(plan)
            call dfftw_destroy_plan(plan)

            !!! normalize
            Gk1_real=Gk1_real/float(nint)

            !!! write values of full-BZ bubble in reduced-BZ bubble
            do ix=0,LQ-1
            do iy=0,ix
               ind=ix*(ix+1)/2+iy+1   !!! index for fully irred. BZ
               chi_bubble(j,ind)=chi_bubble(j,ind)-Gk1_real(ix,iy)
            enddo
            enddo

          ENDDO
       ENDDO
    ENDDO

    chi_bubble=beta*chi_bubble/ &
         ((2.0d0**2*dfloat(Nint)))

    DEALLOCATE(w,Gk1,Gr1,Gk2,Gr2)
    DEALLOCATE(Gk_tmp,Gk1_real,Gk2_real,Gk1_imag,Gk2_imag,Gr1_real,Gr2_real,Gr1_imag,Gr2_imag)

  END SUBROUTINE calc_bubble_fft_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_fftw()
      integer :: n, i
      integer*8 :: plan
      double complex, allocatable :: in(:), out(:)

      write(*,*) 'Input data:'
      n = 10
      allocate(in(n))
      allocate(out(n))
      do i=1,n
        in(i) = cmplx(i,0.0)
        write(*,*) in(i)
      enddo

      ! Forward Fourier transform
      call dfftw_plan_dft_1d(plan,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)

      ! normalize
      out = out/sqrt(float(n))

      write(*,*) 'Fourier transform of the input data:'
      do i=1,n
        write(*,*) out(i)
      enddo

      ! Inverse Fourier transform
      call dfftw_plan_dft_1d(plan,n,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)

      ! normalize
      in = in*sqrt(float(n))

      write(*,*) 'Recovered data from inverse Fourier transform:'
      do i=1,n
        write(*,*) real(in(i)/n)
      enddo

  end subroutine test_fftw

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
