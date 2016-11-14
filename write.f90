!this module contains subroutines for writing data
MODULE write
  IMPLICIT NONE

CONTAINS

  SUBROUTINE write_self_loc(fname,Iwbox,beta,self,selfloc_res)
    !writes the DMFT self-energy and the local self-energy from the DGA-equation

    !input:
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: Iwbox
    REAL(KIND=8), INTENT(IN) :: beta
    COMPLEX(KIND=8), DIMENSION(0:), INTENT(IN) :: self,selfloc_res
    !subroutine internal variables
    INTEGER :: j
    REAL(KIND=8) :: pi
    
    pi=dacos(-1.0d0)
    OPEN(30,file=fname,form='formatted',status='replace')
    DO j=0,Iwbox-1
       WRITE(30,'(5f17.10)') dfloat(2*j+1)*pi/beta, &
            dreal(self(j)),dimag(self(j)), &
            dreal(selfloc_res(j)),dimag(selfloc_res(j))
    ENDDO
    CLOSE(30)

  END SUBROUTINE write_self_loc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE write_chi_bubble(fname,Iwbox,chi_bubble)
    !writes the bubble term for a given q

    !input:
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: Iwbox
    COMPLEX(KIND=8), DIMENSION(-Iwbox:), INTENT(IN) :: chi_bubble
    !subroutine internal variables
    INTEGER j

    OPEN(30,file=fname,form='formatted',status='replace')

    DO j=-Iwbox,Iwbox-1
       WRITE (30,'(3f17.10)')real(j),dreal(chi_bubble(j)),dimag(chi_bubble(j))
    ENDDO
    CLOSE(30)

  END SUBROUTINE write_chi_bubble

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE write_chi(fname,LQ,Qv,lambda,chi)
    !writes the physical susceptibilities

    !input:
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: LQ
    REAL(KIND=8), INTENT(IN) :: lambda
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: Qv
    COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: chi
    !subroutine internal variables
    INTEGER :: ix,iy
    COMPLEX(KIND=8) :: chival,chival_lambda

    OPEN(30,file=fname,form='formatted',status='replace')

    DO ix=0,LQ-1
       DO iy=0,ix
          chival=chi(ix*(ix+1)/2+iy+1)
          chival_lambda=1.0d0/(1.0d0/chival+lambda)
          WRITE (30,'(6f30.20)')Qv(ix),Qv(iy), &
               dreal(chival_lambda),dimag(chival_lambda), &
               dreal(chival),dimag(chival)
       ENDDO
    ENDDO
    CLOSE(30)

  END SUBROUTINE write_chi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE write_sigma(Iwbox,k_number,beta,selflist_res,selflistch_res, &
       selflistsp_res,selflistrest_res,self,selfloc_res)
    !writes the DGA self-energye (total, charge-contribution, spin-contribution, rest)

    !input
    INTEGER, INTENT(IN) :: Iwbox,k_number
    REAL(KIND=8), INTENT(IN) :: beta
    COMPLEX(KIND=8), DIMENSION(:,0:), INTENT(IN) :: selflist_res,selflistch_res
    COMPLEX(KIND=8), DIMENSION(:,0:), INTENT(IN) :: selflistsp_res,selflistrest_res
    COMPLEX(KIND=8), DIMENSION(-2*Iwbox:), INTENT(IN) :: self
    COMPLEX(KIND=8), DIMENSION(0:), INTENT(IN) :: selfloc_res
    !subroutine internal variables
    CHARACTER(LEN=50) :: fname
    INTEGER :: i1,j
    REAL(KIND=8) :: pi
    COMPLEX(KIND=8) :: self_out,selfrest_out

    pi=dacos(-1.0d0)
    DO i1=1,k_number
       WRITE(fname,'(A13,I6.6,A4)'),'klist/SELF_Q_',i1,'.dat'
       OPEN(400+i1,file=fname,form='formatted',status='replace')
       WRITE(400+i1,'(A25)')'# iv_n S S_ch S_sp S_rest'
       DO j=0,iwbox-1
          self_out=selflist_res(i1,j)-selfloc_res(j)+self(j)
          selfrest_out=selflistrest_res(i1,j)-selfloc_res(j)+self(j)
          WRITE(400+i1,'(9f17.10)')dfloat(2*j+1)*pi/beta, &
               dreal(self_out),dimag(self_out), &
               dreal(selflistch_res(i1,j)),dimag(selflistch_res(i1,j)), &
               dreal(selflistsp_res(i1,j)),dimag(selflistsp_res(i1,j)), &
               dreal(selfrest_out),dimag(selfrest_out)
       ENDDO
       CLOSE(400+i1)
    ENDDO

    END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE write_energies(fname,Iwbox,beta,en)
      !this subroutine writes kinetic and potential energies
      
      !input:
      CHARACTER(LEN=*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: Iwbox
      REAL(KIND=8), INTENT(IN) :: beta
      COMPLEX(KIND=8), DIMENSION(-1:,:) :: en
      !subroutine internal variables
      INTEGER :: j
      REAL(KIND=8) :: pi

      pi=dacos(-1.0d0)
      OPEN (30,file=fname,form='formatted',status='replace')
      DO j=0,Iwbox-1
         WRITE(30,'(9f17.10)')dfloat(2*j+1)*pi/beta, &
              dreal(en(j,1)),dimag(en(j,1)),dreal(en(j,2)),dimag(en(j,2)), &
              dreal(en(j,3)),dimag(en(j,3)),dreal(en(j,4)),dimag(en(j,4))
      ENDDO
      CLOSE(30)
      
    END SUBROUTINE write_energies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE write_occupation(fname,epssteps,epsmin,epsmax,n_eps)
      !this subroutine writes the occupation with/without U=0 contribution
      
      !input:
      CHARACTER(LEN=*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: epssteps
      REAL(KIND=8), INTENT(IN) :: epsmin,epsmax
      COMPLEX(KIND=8), DIMENSION(0:,:) :: n_eps
      !subroutine internal variables
      INTEGER :: j

      OPEN (30,file=fname,form='formatted',status='replace')
      DO j=0,epssteps
         WRITE(30,'(5f17.10)')epsmin+dfloat(j)*(epsmax-epsmin)/dfloat(epssteps), &
              dreal(n_eps(j,1)),dimag(n_eps(j,1)),dreal(n_eps(j,2)),dimag(n_eps(j,2))
      ENDDO
      CLOSE(30)

    END SUBROUTINE write_occupation
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE write
