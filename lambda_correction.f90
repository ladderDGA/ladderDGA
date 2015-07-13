!the value of lambda for the lambda-correction is calculated
MODULE lambda_correction
  IMPLICIT NONE

  !parameters for charge and spin correction
  
  ! number of steps for newton iteration
  INTEGER, PARAMETER, DIMENSION(2), PRIVATE :: xsteps=(/1000,1000/)  
  !precision goal for newton iteration
  REAL(KIND=8), PARAMETER, DIMENSION(2), PRIVATE :: prec=(/1d-14,1d-14/)
  !starting value for lambda: -min(chi^-1)+deltino
  REAL(KIND=8), PARAMETER, DIMENSION(2), PRIVATE :: deltino=(/0.10d0,0.10d0/)

CONTAINS

  REAL(KIND=8) FUNCTION lambda(fname,chsp,myid,i,sum_ind,beta, &
       LQ,Chi_loc_sum,Chi_inv_min,chi_x0)
    INCLUDE 'mpif.h'
    !calculates the lambda-correction value in the corresponding channel chsp (chsp=1->charge, chsp=2->sp)
    !results for lambda are writeen in file fname
    
    !input: chsp=1->calculation for charge channel, chsp=2->calculation for spin-channel
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: chsp,myid,i,sum_ind,LQ
    REAL(KIND=8), INTENT(IN) :: beta
    COMPLEX(KIND=8), INTENT(IN) :: Chi_loc_sum,Chi_inv_min
    COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: chi_x0
    !function internal variables
    INTEGER :: xstep,ix,iy,ind
    REAL(KIND=8) :: x0,x,norm,a,b,deltino_used
    COMPLEX(KIND=8) :: Chi_DGA,Chi_DGA_der,Chi_DGA_der2
    COMPLEX(KIND=8) :: Chi_DGA_sum,Chi_DGA_der_sum,Chi_DGA_der2_sum
    !MPI-variable
    INTEGER ::  ierror
    
    OPEN(30,FILE=fname,form='formatted',status='replace')
    deltino_used=deltino(chsp)
    x0=-1.0d0/dreal(Chi_inv_min)+deltino_used
    x=x0
    DO xstep=0,xsteps(chsp)
       Chi_DGA=(0d0,0d0)
       Chi_DGA_der=(0d0,0d0)
       Chi_DGA_der2=(0d0,0d0)
       norm=0d0
       DO ix=0,LQ-1
          a=1.d0
          IF(ix.EQ.0.OR.ix.EQ.LQ-1) THEN
             a=a*0.5d0
          ENDIF
          DO iy=0,ix
             b=a
             IF(iy.eq.0.or.iy.eq.LQ-1) THEN
                b=b*0.5d0
             ENDIF
             
             b=b*dfloat(2/((1+iy)/(1+ix)+1))

             ind=ix*(ix+1)/2+iy+1

             IF((i.GE. -sum_ind).AND.(i.LE.sum_ind)) THEN
                Chi_DGA=Chi_DGA+4*b/ &
                     (1.0d0/Chi_x0(ind)+x)
                Chi_DGA_der=Chi_DGA_der-4*b/ &
                     (1.0d0/Chi_x0(ind)+x)**2
                Chi_DGA_der2=Chi_DGA_der2+8*b/ &
                     (1d0/Chi_x0(ind)+x)**3
             ENDIF
             norm=norm+4*b
          ENDDO
       ENDDO
       Chi_DGA=Chi_DGA/(norm*beta)
       Chi_DGA_der=Chi_DGA_der/(norm*beta)
       Chi_DGA_der2=Chi_DGA_der2/(norm*beta)
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
       CALL MPI_ALLREDUCE(Chi_DGA,Chi_DGA_sum, &
            1,MPI_COMPLEX16,MPI_SUM,MPI_COMM_WORLD,ierror)
       CALL MPI_ALLREDUCE(Chi_DGA_der,Chi_DGA_der_sum, &
            1,MPI_COMPLEX16,MPI_SUM,MPI_COMM_WORLD,ierror)
       CALL MPI_ALLREDUCE(Chi_DGA_der2,Chi_DGA_der2_sum, &
            1,MPI_COMPLEX16,MPI_SUM,MPI_COMM_WORLD,ierror)
       IF(myid.EQ.0) THEN
          WRITE(6,*),'xstep=',xstep+1,' of ',xsteps
          WRITE(6,*),', Chi_DGA_sum=',Chi_DGA_sum
          WRITE(6,*),', Chi_AIM_sum=',Chi_loc_sum
          WRITE(30,'(I10,5f25.20,4f15.10)')xstep+1,x, &
               Chi_DGA_sum,Chi_loc_sum, &
               Chi_DGA_der_sum,Chi_DGA_der2_sum
       ENDIF
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
       
       IF(ABS(Chi_DGA_sum-Chi_loc_sum).LT.prec(chsp)) THEN
          EXIT
       ENDIF
         
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
       
       x=x0-dreal((Chi_DGA_sum-Chi_loc_sum)/Chi_DGA_der_sum)
       
       IF (x.LT.x0) THEN
          deltino_used=deltino_used/2.0d0
          x0=-1.0d0/dreal(Chi_inv_min)+deltino_used
          x=x0
       ELSE
          x0=x
       ENDIF
       
       Chi_DGA_der2_sum=(0d0,0d0)
       Chi_DGA_der_sum=(0d0,0d0)
       Chi_DGA_sum=(0d0,0d0)
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
    ENDDO
    lambda=x
    CLOSE(30)
  END FUNCTION lambda

END MODULE lambda_correction
