PROGRAM make_klist
  USE dispersion
  IMPLICIT NONE
  INTEGER, PARAMETER :: k_range=50
  INTEGER :: ix,iy,iz
  REAL(KIND=8) :: pi,kx,ky,kz,energy

  OPEN(30,file='klist.dat',form='formatted',status='replace')

  pi=dacos(-1.0d0)

  DO ix=0,k_range
     kx=dfloat(ix)*pi/dfloat(k_range)
     DO iy=0,ix
        ky=dfloat(iy)*pi/dfloat(k_range)
        DO iz=0,iy
           kz=dfloat(iz)*pi/dfloat(k_range)
           energy=eps(dcos(kx),dcos(ky),dcos(kz),1d0,1d0,1d0,0d0,0d0,0d0,0d0,0d0,0d0)
           WRITE(30,'(4f20.15)')kx,ky,kz,energy
        ENDDO
     ENDDO
  ENDDO

END PROGRAM make_klist
