PROGRAM make_klist
  USE dispersion
  IMPLICIT NONE
  INTEGER, PARAMETER :: k_range=50
  INTEGER :: ix,iy
  REAL(KIND=8) :: pi,kx,ky,energy

  OPEN(30,file='klist.dat',form='formatted',status='replace')

  pi=dacos(-1.0d0)

  DO ix=0,k_range
     kx=dfloat(ix)*pi/dfloat(k_range)
     DO iy=0,ix
        ky=dfloat(iy)*pi/dfloat(k_range)
        energy=eps(dcos(kx),dcos(ky),1d0,1d0,0d0,0d0,0d0,0d0)
        WRITE(30,'(3f20.15)')kx,ky,energy
     ENDDO
  ENDDO

END PROGRAM make_klist
