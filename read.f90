!This module contains all routines for reading in parameters and datafiles
MODULE read
  IMPLICIT NONE

CONTAINS
  
  SUBROUTINE read_parameters(fname,uhub,xmu,beta,nden, &
       Iwbox,Iwbox_bose,shift,LQ,Nint,k_number,sigma_only,chi_only,lambdaspin_only,sumallch,sumallsp,xch_so,xsp_so)
    
    !input: filename of parameter file
    CHARACTER(LEN=*), INTENT(IN) :: fname
    !input: Parameters
    REAL(KIND=8), INTENT(OUT) :: uhub,xmu,beta,nden
    INTEGER, INTENT(OUT) :: Iwbox,Iwbox_bose,shift,LQ,Nint,k_number
    LOGICAL, INTENT(OUT) :: sigma_only,chi_only,lambdaspin_only,sumallch,sumallsp
    REAL(KIND=8), INTENT(OUT) :: xch_so,xsp_so
    
    !reading AIM/DMFT and lambda-correction parameters from file "fname"
    OPEN(30,file=fname,form='formatted',status='old')
    READ(30,*) 
    READ(30,*) uhub, xmu, beta, nden
    READ(30,*) 
    READ(30,*) Iwbox,Iwbox_bose,shift
    READ(30,*) 
    READ(30,*) LQ,Nint,k_number
    READ(30,*) 
    READ(30,*) sigma_only, xch_so, xsp_so
    READ(30,*) 
    READ(30,*) chi_only
    READ(30,*) 
    READ(30,*) lambdaspin_only
    READ(30,*)
    READ(30,*) sumallch,sumallsp
    CLOSE(30)
    
  END SUBROUTINE read_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_parameters_thermo(fname,calcU0,calcDMFT,calcDGA, &
       uhub,mu,beta,nden,ap,Iwbox,k_range,epssteps,epsmin,epsmax,fermicut)
    
    !input: filename of parameter file
    CHARACTER(LEN=*), INTENT(IN) :: fname
    !input: Parameters
    LOGICAL, INTENT(OUT) :: calcU0,calcDMFT,calcDGA
    REAL(KIND=8), INTENT(OUT) :: uhub,mu,beta,nden
    INTEGER, INTENT(OUT) :: Iwbox,k_range,epssteps,ap
    REAL(KIND=8), INTENT(OUT) :: epsmin,epsmax,fermicut
    
    !reading AIM/DMFT and lambda-correction parameters from file "fname"
    OPEN(30,file=fname,form='formatted',status='old')
    READ(30,*) 
    READ(30,*) uhub, mu, beta, nden, ap  !(ap...Number of Anderson Parameters in DMFT)
    READ(30,*) 
    READ(30,*) Iwbox
    READ(30,*) 
    READ(30,*) k_range, epssteps
    READ(30,*) 
    READ(30,*) epsmin, epsmax, fermicut
    READ(30,*) 
    READ(30,*) calcU0,calcDMFT,calcDGA
    CLOSE(30)
    
  END SUBROUTINE read_parameters_thermo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_anderson_parameters(fname,ap,v)
    
    !input:
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: ap
    !output:
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: v
    !subroutine internal variables
    INTEGER :: i
    REAL(KIND=8) :: indat

    OPEN(30,file=fname,form='formatted',status='old')
    DO i=1,10+ap
       READ(30,*)
    ENDDO
    DO i=1,ap
       READ(30,*)indat
       v(i)=indat
    ENDDO
    CLOSE(30)

  END SUBROUTINE read_anderson_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_local_sigma(Iwbox,gww,self)
    !reads local interacting and non-interacting green's functions and calculates the local self-energy
    !this routine might be replaced be reading directly a local or nonlocal self-enregy!!!
    
    !input: number of fermionic frequencies
    INTEGER, INTENT(IN) :: Iwbox
    !output: (local) self-energy
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(OUT) :: self
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(OUT) :: gww
    !subroutine internal variables
    INTEGER :: i
    REAL(KIND=8) :: freq,re,im
    COMPLEX(KIND=8) :: g0wwand
    
    OPEN(30,file='gm_wim',form='formatted',status='old')
    OPEN(31,file='g0mand',form='formatted',status='old')
    
    DO i=0,2*Iwbox-1
       READ(30,*)freq,re,im
       gww(i)=dcmplx(re,im)
       gww(-i-1)=dconjg(gww(i))
       READ(31,*)freq,re,im
       g0wwand=dcmplx(re,im)
       
       self(i)=g0wwand-dcmplx(1.0d0,0.0d0)/gww(i)
       self(-i-1)=dconjg(self(i))
    ENDDO
    
    CLOSE(30)
    CLOSE(31)
    
  END SUBROUTINE read_local_sigma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_fupdown(fname,Iwbox,i,beta,gww,chich_loc,chisp_loc,fupdown)
    !reads the local generalized susceptibility chi_updown for one bosonic frequency i and
    !extracts the vertex fupdown
    !->will be changed to the number of bosonic frequencies treated by the actual rank!!!
    
    !input: fname, number of fermionic frequencies,bosonic frequence,inverse temperature, local green's function
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: Iwbox,i
    REAL(KIND=8), INTENT(IN) :: beta
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: gww
    !output
    COMPLEX(KIND=8), INTENT(OUT) :: chich_loc,chisp_loc
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(OUT) :: fupdown
    !subroutine internal variables
    INTEGER :: j,k,ierror
    REAL(KIND=8) :: freq,reup,imup,redown,imdown
    

    OPEN(30,file=fname,form='formatted',status='old')
    chich_loc=dcmplx(0.0d0,0.0d0)
    chisp_loc=dcmplx(0.0d0,0.0d0)
    DO j=-Iwbox,Iwbox-1
       DO k=-Iwbox,Iwbox-1
          READ(30,*)freq,freq,freq,reup,imup,redown,imdown
          fupdown(j,k)=dcmplx(redown,imdown)/(beta**2* &
               gww(j)*gww(i+j)*gww(k)*gww(i+k))
          chich_loc=chich_loc+dcmplx(reup+redown,imup+imdown)
          chisp_loc=chisp_loc+dcmplx(reup-redown,imup-imdown)
       ENDDO
    ENDDO
    chich_loc=chich_loc/beta**2
    chisp_loc=chisp_loc/beta**2

    CLOSE(30)
    
  END SUBROUTINE read_fupdown
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_gamma(fname,Iwbox,gammach,gammasp)
    !reads local irreducible vertex in the spin- and charge-channel for one bosonic frequency
    !->will be changed to the number of bosonic frequencies treated by the actual rank!!!
    
    !input: fname,number of fermionic frequencies
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: Iwbox
    !output
    COMPLEX(KIND=8), DIMENSION(-Iwbox:,-Iwbox:), INTENT(OUT) :: gammach,gammasp
    !subroutine internal variables
    INTEGER :: j,k
    REAL(KIND=8) :: freq,rech,imch,resp,imsp

    OPEN(30,file=fname,form="formatted",status="old")
    DO j=-Iwbox,Iwbox-1
       DO k=-Iwbox,Iwbox-1
          READ(30,*)freq,freq,freq, &
               rech,imch,resp,imsp
          gammach(j,k)=dcmplx(-rech,-imch)  !asymptotics of gammach -> -U here
          gammasp(j,k)=dcmplx(-resp,-imsp)  !asymptotics of gammasp -> +U here
       ENDDO
    ENDDO
    
    CLOSE(30)
    
  END SUBROUTINE read_gamma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_klist(fname,k_number,klist)
    !reads the external k-points, for which sigma should be calculated,  from a list
    
    !input: fname, number of external k-points
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: k_number
    !output: list of all k-points
    REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: klist
    !subroutine internal variables
    INTEGER :: i
    REAL(KIND=8) :: kx,ky,kz
    
    OPEN(30,file=fname,form='formatted',status='old')
    DO i=1,k_number
       READ(30,*)kx,ky,kz
       klist(i,1)=kx
       klist(i,2)=ky
       klist(i,3)=kz
    ENDDO
    CLOSE(30)
    
  END SUBROUTINE read_klist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_chix0(fname,LQ,chi_x0)
    !writes the physical susceptibilities

    !input:
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: LQ
    !output
    COMPLEX(KIND=8), DIMENSION(:), INTENT(OUT) :: chi_x0
    !subroutine internal variables
    INTEGER :: ix,iy,iz
    REAL(KIND=8) :: mom,re,im,rel0,iml0

    OPEN(30,file=fname,form='formatted',status='old')

    DO ix=0,LQ-1
       DO iy=0,ix
          DO iz=0,iy
             READ (30,*)mom,mom,mom,re,im,rel0,iml0
             chi_x0(ix*(ix+1)*(ix+2)/6+iy*(iy+1)/2+iz+1)=dcmplx(rel0,iml0)
          ENDDO
       ENDDO
    ENDDO
    CLOSE(30)

  END SUBROUTINE read_chix0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE read_sigma(Iwbox,k_number,self_k)
    !reads the nonlocal self-energy from files ./klist/SELF_Q_#.dat
    
    !input: number of frequencies and k-points
    INTEGER, INTENT(IN) :: Iwbox,k_number
    !output: self-energy for momentum k
    COMPLEX(KIND=8), DIMENSION(:,0:), INTENT(OUT) :: self_k
    !subroutine internal variables
    CHARACTER(LEN=50) :: fname
    INTEGER :: i1,j
    REAL(KIND=8) :: freq,re,im
    
    DO i1=1,k_number
       WRITE(fname,'(A13,I6.6,A4)'),'klist/SELF_Q_',i1,'.dat'
       OPEN(400+i1,file=fname,form='formatted',status='old')
       READ(400+i1,*)
       DO j=0,Iwbox-1
          READ(400+i1,*)freq,re,im
          self_k(i1,j)=dcmplx(re,im)
       ENDDO
       CLOSE(400+i1)
    ENDDO
   
  END SUBROUTINE read_sigma

END MODULE read
