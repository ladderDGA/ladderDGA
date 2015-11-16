!this module includes subroutines for calculating energies and occupation numbers
MODULE thermo
  USE dispersion
  IMPLICIT NONE

  REAL(KIND=8), PARAMETER, PRIVATE :: deltino=0.010d0

CONTAINS

  SUBROUTINE calc_energies(Iwbox,k_min,k_max,k_range,kcount,epssteps, &
       uhub,mu,beta,nden,fermicut,epsmin,epsmax,dcol,self,en,n_eps,testsum,testsum_eps)
    !calculates kinetic (en(:,1:2)) and potential (en(:,3:4)) energies as function of cut-off frequency

    !input
    INTEGER, INTENT(IN) :: Iwbox,k_min,k_max,k_range,epssteps
    INTEGER, DIMENSION(k_min:,:), INTENT(IN) :: kcount
    REAL(KIND=8), INTENT(IN) :: uhub,mu,beta,nden,fermicut,epsmin,epsmax
    REAL(KIND=8), DIMENSION(k_min:,:), INTENT(IN) :: dcol
    COMPLEX(KIND=8), DIMENSION(k_min:,0:), INTENT(IN) :: self
    !output
    REAL(KIND=8), INTENT(OUT) :: testsum,testsum_eps
    COMPLEX(KIND=8), DIMENSION(-1:,:), INTENT(OUT) :: en
    COMPLEX(KIND=8), DIMENSION(0:,:), INTENT(OUT) :: n_eps
    !subroutine internal variables
    INTEGER :: j,ix,iy,iz,ind,ieps
    REAL(KIND=8) :: pi,sigma_hartree,energy,mult
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: weight
    COMPLEX(KIND=8) :: ekin_correct1,ekin_correct2,epot_correct1,epot_correct2,gk
    COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: w,egrid,n_k
    COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: enfreq

    pi=dacos(-1.0d0)
    sigma_hartree=uhub*nden/2

    ALLOCATE(w(0:Iwbox-1))
    DO j=0,Iwbox-1
       w(j)=dcmplx(0.0d0,dfloat(2*j+1))*pi/beta
    ENDDO
    ALLOCATE(weight(0:k_range))
    weight=1.0d0
    weight(0)=0.50d0
    weight(k_range)=0.50d0
    
    ALLOCATE(enfreq(-1:Iwbox-1,4))
    en=dcmplx(0.0d0,0.0d0)
    enfreq=dcmplx(0.0d0,0.0d0)
    ekin_correct1=dcmplx(0.0d0,0.0d0)
    ekin_correct2=dcmplx(0.0d0,0.0d0)
    epot_correct1=dcmplx(0.0d0,0.0d0)
    epot_correct2=dcmplx(0.0d0,0.0d0)
    n_eps=dcmplx(0.0d0,0.0d0)
    testsum=0.0d0
    testsum_eps=0.0d0
    ALLOCATE(n_k(2))
    ALLOCATE(egrid(0:epssteps))
    DO ieps=0,epssteps
       egrid(ieps)=epsmin+dfloat(ieps)*(epsmax-epsmin)/dfloat(epssteps)
    ENDDO

    DO ind=k_min,k_max
       energy=eps(dcol(ind,1),dcol(ind,2),dcol(ind,3),1d0,1d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0)
       n_k=dcmplx(0.0d0,0.0d0)
       IF (DABS(energy).LT.fermicut) THEN
          ix=kcount(ind,1)
          iy=kcount(ind,2)
          iz=kcount(ind,3)
          mult=weight(ix)*weight(iy)*weight(iz)*dfloat(6/ &
               ((1+iy)/(1+ix)+(1+iz)/(1+iy)+ &
               3*((1+iz)/(1+ix))+1))
          DO j=0,Iwbox-1
             !define the Greens-funktion
             gk=1.0d0/(w(j)+mu-energy-self(ind,j))
             !plain sum for kinetic energy (1/nu-term does not contribute in finite symmetric sum!)
             enfreq(j,1)=enfreq(j-1,1)+2.0d0*mult*energy*dreal(gk)
             !1/nu^2-term of the summand subtracted and analytically evaluated
             enfreq(j,2)=enfreq(j-1,2)+2.0d0*mult*energy*dreal(gk-(energy+sigma_hartree-mu)/w(j)**2)
             !plain sum for potential energy (1/nu-term does not contribute in finite symmetric sum!)
             enfreq(j,3)=enfreq(j-1,3)+2.0d0*mult*dreal(self(ind,j)*gk)
             !1/nu^2-term subtracted and analytically evaluated
             enfreq(j,4)=enfreq(j-1,4)+2.0d0*mult*dreal(self(ind,j)*gk- &
                  (uhub**2*0.50d0*nden*(1.0d0-0.50d0*nden)+ &
                  sigma_hartree*(energy+sigma_hartree-mu))/w(j)**2)
             
             !calculation of occupation number
             !full occupation number (1/nu-term does not contribute in finite symmetric sum!)
             n_k(1)=n_k(1)+2.0d0*mult*dreal(gk)
             !occupation number without U=0 contribution
             n_k(2)=n_k(2)+2.0d0*mult*dreal(gk*(self(ind,j)-mu)/(w(j)-energy)- &
                  (sigma_hartree-mu)/w(j)**2)
          ENDDO
          en=en+enfreq
          !Asumptotic corrections of energies
          ekin_correct1=ekin_correct1+mult*energy
          ekin_correct2=ekin_correct2+mult*energy*(energy+sigma_hartree-mu)
          epot_correct1=epot_correct1+mult*sigma_hartree
          epot_correct2=epot_correct2+mult*(uhub**2*0.50d0*nden*(1.0d0-0.50d0*nden)- &
               sigma_hartree*(mu-energy-sigma_hartree))
          !Asymptotic correction of n_k
          n_k(1)=n_k(1)+0.50d0*beta*mult
          n_k(2)=n_k(2)-0.250d0*(sigma_hartree-mu)*mult*beta**2
          !Loop over energies for occupation numbers
          DO ieps=0,epssteps
             n_eps(ieps,:)=n_eps(ieps,:)+n_k*deltino/((egrid(ieps)-energy)**2+deltino**2)
          ENDDO
          testsum=testsum+mult
          testsum_eps=testsum_eps+mult*energy**2
       ENDIF
    ENDDO

    en(:,1)=2.0d0*(en(:,1)/beta+ekin_correct1*0.50d0)
    en(:,2)=2.0d0*(en(:,2)/beta+0.50d0*ekin_correct1- &
         0.250d0*beta*ekin_correct2)
    en(:,3)=(en(:,3)/beta+0.5*epot_correct1)
    en(:,4)=(en(:,4)/beta+0.5*epot_correct1- &
         0.250d0*beta*epot_correct2)
    !Normalization of n_eps; factor (1/pi) comes from Lorentzian!
    n_eps=2.0d0*n_eps/(pi*beta)

  END SUBROUTINE calc_energies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE thermo
