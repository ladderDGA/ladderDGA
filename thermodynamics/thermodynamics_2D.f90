!calculation of kinetic and potential energy and occuptaion-number for a given energy (for U=0, DMFT and DGA)
PROGRAM thermodynamics
  USE read
  USE write
  USE dispersion
  USE vardef_thermo
  USE thermo

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  !MPI initialization
  CALL MPI_INIT(ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)
  
  !allocate arrays for distributing the k-points among the processes
  ALLOCATE(recvnumber(0:nprocs-1))
  ALLOCATE(offsetnum(0:nprocs-1))

  !only rank 0 reads the parameters
  IF (myid.EQ.0) THEN
     CALL read_parameters_thermo('ladderDGA_thermo.in',calcU0,calcDMFT,calcDGA, &
          uhub,mu,beta,nden,ap,Iwbox,k_range,epssteps,epsmin,epsmax,fermicut)
  ENDIF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(uhub,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(mu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(beta,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(nden,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(Iwbox,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(k_range,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(epssteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(epsmin,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(epsmax,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(fermicut,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(calcU0,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(calcDMFT,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(calcDGA,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

  knumber_tot=(k_range+1)*(k_range+2)/2
  knumber=(knumber_tot)/nprocs
  krest=MOD(knumber_tot,nprocs)
  DO i=0,nprocs-1
     IF (i.LT.krest) THEN
        recvnumber(i)=knumber+1
        offsetnum(i)=i*(knumber+1)
     ELSE
        recvnumber(i)=knumber
        offsetnum(i)=krest*(knumber+1)+(i-krest)*knumber
     ENDIF
  ENDDO
  IF (myid.LT.krest) THEN
     knumber=knumber+1
  ENDIF
  IF (myid.LT.krest) THEN
     koffset=myid*knumber
  ELSE
     koffset=krest*(knumber+1)+(myid-krest)*knumber
  ENDIF
  k_min=koffset+1
  k_max=koffset+knumber
  !mapping: total index for a given k-point -> (ix,iy)-coordinate in k-grid
  ALLOCATE(kcount(k_min:k_max,2))
  DO ix=0,k_range
     DO iy=0,ix
        ind=ix*(ix+1)/2+iy+1
        IF ((ind.GE.k_min).AND.(ind.LE.k_max)) THEN
           kcount(ind,1)=ix
           kcount(ind,2)=iy
        ENDIF
     ENDDO
  ENDDO
  ALLOCATE(klist_tot(knumber_tot,2))
  ALLOCATE(klist(k_min:k_max,2))
  IF (myid.EQ.0) THEN
     CALL read_klist('klist.dat',knumber_tot,klist_tot)
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  CALL MPI_SCATTERV(klist_tot(:,1),recvnumber,offsetnum,MPI_REAL8, &
       klist(:,1),knumber,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_SCATTERV(klist_tot(:,2),recvnumber,offsetnum,MPI_REAL8, &
       klist(:,2),knumber,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  DEALLOCATE(klist_tot)

  !cos(k) for all k-values of each process
  ALLOCATE(dcol(k_min:k_max,2))
  DO j=k_min,k_max
     dcol(j,1)=dcos(klist(j,1))
     dcol(j,2)=dcos(klist(j,2))
  ENDDO

  !self-energy (self) for which kinetic and potential energy (en(:,1:4)) are calculated
  ALLOCATE(self(k_min:k_max,0:Iwbox-1))
  ALLOCATE(en(-1:Iwbox-1,4))
  ALLOCATE(en_sum(-1:Iwbox-1,4))
  ALLOCATE(n_eps(0:epssteps,2))
  ALLOCATE(n_eps_sum(0:epssteps,2))

  !Self-energy for U=0
  IF (calcU0) THEN
     self=dcmplx(0.0d0,0.0d0)
     CALL calc_energies(Iwbox,k_min,k_max,k_range,kcount,epssteps, &
          0.0d0,0.0d0,beta,nden,fermicut,epsmin,epsmax,dcol,self,en,n_eps,testsum,testsum_eps)
     CALL MPI_REDUCE(en,en_sum,4*(Iwbox+1),MPI_COMPLEX16, &
          MPI_SUM,0,MPI_COMM_WORLD,ierror)
     en_sum=en_sum/dfloat(k_range)**2
     CALL MPI_REDUCE(n_eps,n_eps_sum,2*(epssteps+1),MPI_COMPLEX16, &
          MPI_SUM,0,MPI_COMM_WORLD,ierror)
     n_eps_sum=n_eps_sum/dfloat(k_range)**2
     IF (myid.EQ.0) THEN
        CALL write_energies('energiesU0.dat',Iwbox,beta,en_sum)
        CALL write_occupation('occupationU0.dat',epssteps,epsmin,epsmax,n_eps_sum)
     ENDIF
  ENDIF

  !Local self-energy of DMFT
  IF (calcDMFT) THEN
     !allocate local arrays for input data
     ALLOCATE(self_loc(-2*Iwbox:2*Iwbox-1))
     !read local green's function and local energy of DMFT (only rank 0)
     IF (myid.eq.0) THEN
        ALLOCATE(gww(-2*Iwbox:2*Iwbox-1))
        CALL read_local_sigma(Iwbox,gww,self_loc)
     ENDIF
     !broadcast local green's function and local self-energy to all ranks
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
     CALL MPI_BCAST(self_loc,4*Iwbox,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierror)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
     DO j=k_min,k_max
        self(j,:)=self_loc(0:Iwbox-1)
     ENDDO
     CALL calc_energies(Iwbox,k_min,k_max,k_range,kcount,epssteps, &
          uhub,mu,beta,nden,fermicut,epsmin,epsmax,dcol,self,en,n_eps,testsum,testsum_eps)
     CALL MPI_REDUCE(en,en_sum,4*(Iwbox+1),MPI_COMPLEX16, &
          MPI_SUM,0,MPI_COMM_WORLD,ierror)
     en_sum=en_sum/dfloat(k_range)**2
     CALL MPI_REDUCE(n_eps,n_eps_sum,2*(epssteps+1),MPI_COMPLEX16, &
          MPI_SUM,0,MPI_COMM_WORLD,ierror)
     n_eps_sum=n_eps_sum/dfloat(k_range)**2
     IF (myid.EQ.0) THEN
        CALL write_energies('energiesDMFT.dat',Iwbox,beta,en_sum)
        CALL write_occupation('occupationDMFT.dat',epssteps,epsmin,epsmax,n_eps_sum)
     ENDIF
     IF (myid.EQ.0) THEN
        ALLOCATE(v(1:ap))
        CALL read_anderson_parameters('hubb.andpar',ap,v)
        vsum=0.0d0
        DO j=1,ap
           vsum=vsum+(v(j))**2
        ENDDO
        en=dcmplx(0.0d0,0.0d0)
        sigma_hartree=uhub*nden/2.0d0
        !DMFT self-energy without k-sum (calculated with g_loc=gww)
        DO j=0,Iwbox-1
           w=dcmplx(0.0d0,dfloat(2*j+1))*pi/beta
           !plain sum for kinetic energy
           en(j,1)=en(j-1,1)+2.0d0*dreal((w+mu-self_loc(j))*gww(j)-1.0d0)
           !1/nu^2-term subtracted and analytically evaluated
           en(j,2)=en(j-1,2)+2.0d0*dreal((w+mu-self_loc(j))*gww(j)-1.0d0-vsum/w**2)
           !plain sum for potential energy (1/nu-term does not contribute in finite symmetric sum!)
           en(j,3)=en(j-1,3)+2.0d0*dreal(self_loc(j)*gww(j))
           !1/nu^2-term subtracted and analytically evaluated
           en(j,4)=en(j-1,4)+2.0d0*dreal(self_loc(j)*gww(j)- &
                (uhub**2*0.50d0*nden*(1.0d0-0.50d0*nden)+ &
                sigma_hartree*(sigma_hartree-mu))/w**2)
        ENDDO
        en(:,1)=en(:,1)*2.0d0/beta
        en(:,2)=(en(:,2)-vsum*beta**2/4)*2.0d0/beta
        en(:,3)=en(:,3)/beta+0.50d0*sigma_hartree
        en(:,4)=en(:,4)/beta+0.50d0*sigma_hartree- &
             0.250d0*beta*(uhub**2*0.50d0*nden*(1.0d0-0.50d0*nden)- &
             sigma_hartree*(mu-sigma_hartree))
        CALL write_energies('energies_dmft_nusum.dat',Iwbox,beta,en)
        !DMFT occupation without k-sum -> has to be multiplied with D(eps)
        DO ieps=0,epssteps
           epoint=epsmin+dfloat(ieps)*(epsmax-epsmin)/dfloat(epssteps)
           n_eps(ieps,:)=dcmplx(0.0d0,0.0d0)
           DO j=0,Iwbox-1
              w=dcmplx(0.0d0,dfloat(2*j+1))*pi/beta
              g=1.0d0/(w+mu-epoint-self(1,j))
              n_eps(ieps,1)=n_eps(ieps,1)+dreal(g*(self(1,j)-mu)/(w-epoint)- &
                   (sigma_hartree-mu)/w**2)
           ENDDO
           n_eps(ieps,1)=2.0d0*(n_eps(ieps,1)*2.0d0/beta-0.250d0*(sigma_hartree-mu)*beta)
        ENDDO
        CALL write_occupation('occup_dmft_nusum.dat',epssteps,epsmin,epsmax,n_eps)
        DEALLOCATE(self_loc)
        DEALLOCATE(gww)
        DEALLOCATE(v)
     ENDIF
  ENDIF

  !Self-energy of DGA
  IF (calcDGA) THEN
     ALLOCATE(self_k_tot(knumber_tot,0:Iwbox-1))
     IF (myid.EQ.0) THEN
        CALL read_sigma(Iwbox,knumber_tot,self_k_tot)
     ENDIF
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
     DO i=0,Iwbox-1
        CALL MPI_SCATTERV(self_k_tot(:,i),recvnumber,offsetnum,MPI_COMPLEX16, &
             self(:,i),knumber,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierror)
     ENDDO
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
     DEALLOCATE(self_k_tot)
     !!!!Add Hartree-term to self_k -> should be removed if already present in the input data!!!!
     self=self+uhub*nden/2.0d0   !Hartree term explicitely added!!
     CALL calc_energies(Iwbox,k_min,k_max,k_range,kcount,epssteps, &
          uhub,mu,beta,nden,fermicut,epsmin,epsmax,dcol,self,en,n_eps,testsum,testsum_eps)
     CALL MPI_REDUCE(en,en_sum,4*(Iwbox+1),MPI_COMPLEX16, &
          MPI_SUM,0,MPI_COMM_WORLD,ierror)
     en_sum=en_sum/dfloat(k_range)**2
     CALL MPI_REDUCE(n_eps,n_eps_sum,2*(epssteps+1),MPI_COMPLEX16, &
          MPI_SUM,0,MPI_COMM_WORLD,ierror)
     n_eps_sum=n_eps_sum/dfloat(k_range)**2
     IF (myid.EQ.0) THEN
        CALL write_energies('energiesDGA.dat',Iwbox,beta,en_sum)
        CALL write_occupation('occupationDGA.dat',epssteps,epsmin,epsmax,n_eps_sum)
     ENDIF
  ENDIF

  CALL MPI_REDUCE(testsum,testsum_sum,1,MPI_REAL8, &
       MPI_SUM,0,MPI_COMM_WORLD,ierror)
  testsum_sum=testsum_sum/dfloat(k_range)**2
  CALL MPI_REDUCE(testsum_eps,testsum_eps_sum,1,MPI_REAL8, &
       MPI_SUM,0,MPI_COMM_WORLD,ierror)
  testsum_eps_sum=testsum_eps_sum/dfloat(k_range)**2
  IF (myid.EQ.0) THEN
     WRITE(6,*)testsum_sum
     WRITE(6,*)testsum_eps_sum
  ENDIF

  CALL MPI_FINALIZE(ierror)

  STOP

END PROGRAM thermodynamics



