! ladder-dga with automatic lambda-correction
! ver 1.0
! for the implementation of the Dyson-Schwinger equation of motion
! cf. PRB 80, 075104 (2009), especially eq. 8
PROGRAM self_k
  USE vardef
  USE read
  USE write
  USE sigma
  USE dispersion
  USE calc_susc
  USE lambda_correction

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'

  
  !MPI initialization
  CALL MPI_INIT(ierror)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)

  !reading parameter from file ladderDGA.in:
  !-) uhub...Hubbard interaction
  !-) mu...chemical potential
  !-) beta...inverse temperature
  !-) nden...average density per lattice site
  !-) Iwbox...number of fermionic Matsubara frequencies
  !-) shift...offset for bosonic frequencies
  !-) LQ...interval (0,pi) for internal q-summation is diveded into (LQ-1) parts (LQ=0 -> 0, LQ->pi)
  !-) Nint...interval (-pi,pi) for internal k' summation is split into Nint subintervals
  !-) k_number...number of external k-points for which sigma_dga is calculated
  !-) sigma_only...only sigma_dga is calculated -> chi(q,omega), lambda_(ch,sp) have to be read from files
  !-) schi_so,xsp_so...lambda-corrections which should be used in the sigma_only-case
  !-) chi_only...only chi(q,omega) is calculated

  !only rank 0 reads the parameters
  IF (myid.EQ.0) THEN
     CALL read_parameters('ladderDGA.in',uhub,mu,beta,nden, &
          Iwbox,shift,LQ,Nint,k_number,sigma_only,chi_only,xch_so,xsp_so)
     !Check parameters
     WRITE(6,*) 'U= ', uhub
     WRITE(6,*) 'MU=',mu
     WRITE(6,*) 'BETA=',beta
     WRITE(6,*) '< n >=', nden 
     WRITE(6,*) 'fermionic frequency box=',Iwbox
     WRITE(6,*) 'shift of bosonic frequency box=',shift
     WRITE(6,*) "number of q-points=",LQ
     WRITE(6,*) "number of k'-intervals=",Nint
     WRITE(6,*) "number of external k-points=",k_number
     WRITE(6,*) 'sigma only=',sigma_only
     WRITE(6,*)'lambda correcting for sigma-only, charge=',xch_so
     WRITE(6,*)'lambda correcting for sigma-only, spin=',xsp_so
     WRITE(6,*) 'chi only=',chi_only
  ENDIF
  !Broadcast parameters to all ranks
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  CALL MPI_BCAST(uhub,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(mu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(beta,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(nden,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(Iwbox,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(shift,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(LQ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(Nint,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(k_number,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(xch_so,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(xsp_so,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(sigma_only,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(chi_only,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

  !determine bosonic frequency index i
  !-> should be generalized to an array of bosonic frequenies if #ranks<#bosonic frequencies!!!!!
  i=myid-Iwbox+1+shift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Local part of the program. Calculates: selfloc, chich_loc, chisp_loc,
  !                                                            chich_loc_sum, chisp_loc_sum,
  !                                                            sum_ind_ch, sum_ind_sp          

  !allocate local arrays for input data
  ALLOCATE(gww(-2*Iwbox:2*Iwbox-1))
  ALLOCATE(self(-2*Iwbox:2*Iwbox-1))
  ALLOCATE(fupdown(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
  ALLOCATE(gammach(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))
  ALLOCATE(gammasp(-Iwbox:Iwbox-1,-Iwbox:Iwbox-1))

  !read local green's function and local energy of DMFT (only rank 0)
  IF (myid.eq.0) THEN
     CALL read_local_sigma(Iwbox,gww,self)
  ENDIF
  !broadcast local green's function and local self-energy to all ranks
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  CALL MPI_BCAST(gww,4*Iwbox,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(self,4*Iwbox,MPI_COMPLEX16,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  !create filename for reading the local susceptibilities (chi_up_up, chi_up_down)
  ALLOCATE(CHARACTER(LEN=14)::fname)
  IF ((myid+shift).LT.10) THEN
     WRITE(fname,'(A8,5Hchi00,I1)')'chi_dir/',myid+shift
  ELSEIF ((myid+shift).LT.100) THEN
     WRITE(fname,'(A8,4Hchi0,I2)')'chi_dir/',myid+shift
  ELSE
     WRITE(fname,'(A8,3Hchi,I3)')'chi_dir/',myid+shift
  ENDIF
  !read local vertex fupdown (one bosonic frequency/rank)
  CALL read_fupdown(fname,Iwbox,i,beta,gww,chich_loc,chisp_loc,fupdown)
  DEALLOCATE(fname)
  !create filename for reading the local irreducible vertices (gammach and gammasp)
  ALLOCATE(CHARACTER(LEN=18)::fname)
  IF ((myid+shift).LT.10) THEN
     WRITE(fname,'(A10,7Hgamma00,I1)')'gamma_dir/',myid+shift
  ELSEIF ((myid+shift).LT.100) THEN
     WRITE(fname,'(A10,6Hgamma0,I2)')'gamma_dir/',myid+shift
  ELSE
     WRITE(fname,'(A10,5Hgamma,I3)')'gamma_dir/',myid+shift
  ENDIF
  !read local irreducible vertices gammach and gammasp (one bosonic frequency/rank)
  CALL read_gamma(fname,Iwbox,gammach,gammasp)
  DEALLOCATE(fname)

  !calculate local self-energy with DGA equation
  !output: selfloc(nu,i), chich_loc(i), chisp_loc(i) for bosonic index i
  ALLOCATE(selfloc(-Iwbox:Iwbox-1))
  CALL calc_self_loc(Iwbox,i,uhub,beta,gww,fupdown,gammach,gammasp, &
       selfloc)
  !sum self(nu,i) over bosonic index i
  ALLOCATE(selfloc_res(-Iwbox:Iwbox-1))
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  CALL MPI_REDUCE(selfloc,selfloc_res,2*Iwbox, MPI_COMPLEX16, &
       MPI_SUM,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  DO j=-Iwbox,Iwbox-1
     selfloc_res(j)=selfloc_res(j)*uhub/beta+uhub/2.d0*nden
  ENDDO
  ! output of local self-energy
  IF (myid.EQ.0) THEN
     CALL write_self_loc('klist/SELF_LOC_parallel',Iwbox,beta,self(-Iwbox:Iwbox-1),selfloc_res)
  ENDIF

  !determine bosonices indices for which chich_loc_res, chisp_loc_res < 0 (sum_ind_ch, sum_ind_sp)
  IF (i.GE.0) THEN
     ind_part_ch=1.0d0/(DSIGN(dfloat(i),dreal(chich_loc)+tolerancech)+0.5d0)
     ind_part_sp=1.0d0/(DSIGN(dfloat(i),dreal(chisp_loc)+tolerancesp)+0.5d0)
  ELSE
     ind_part_ch=1.0d0
     ind_part_sp=1.0d0
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  CALL MPI_ALLREDUCE(ind_part_ch,ind_ch,1,MPI_REAL8, &
       MPI_MIN,MPI_COMM_WORLD,ierror)
  CALL MPI_ALLREDUCE(ind_part_sp,ind_sp,1,MPI_REAL8, &
       MPI_MIN,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  sum_ind_ch=INT(1.0d0+DABS(1.0d0/ind_ch))-1
  sum_ind_sp=INT(1.0d0+DABS(1.0d0/ind_sp))-1

  IF ((sum_ind_ch.EQ.0).OR.(sum_ind_sp.EQ.0)) THEN
     WRITE(6,*)'bad chi_loc'
     WRITE(6,*)sum_ind_ch,sum_ind_sp
     CALL MPI_FINALIZE(ierror)
     STOP
  ENDIF
     
  !set chi to 0 for i > sum_ind and save old values in chich_loc_all and chisp_loc_all
  chich_loc_all=chich_loc
  chisp_loc_all=chisp_loc
  IF (ABS(i).GT.sum_ind_ch) THEN
     chich_loc=dcmplx(0.0d0,0.0d0)
  ENDIF
  IF (ABS(i).GT.sum_ind_sp) THEN
     chisp_loc=dcmplx(0.0d0,0.0d0)
  ENDIF

  !perform sum of chich_loc and chisp_loc (via MPI_REDUCE)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  CALL MPI_ALLREDUCE(chich_loc,chich_loc_sum,1,MPI_REAL8, &
       MPI_SUM,MPI_COMM_WORLD,ierror)
  CALL MPI_ALLREDUCE(chisp_loc,chisp_loc_sum,1,MPI_REAL8, &
       MPI_SUM,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  chich_loc_sum=chich_loc_sum/beta
  chisp_loc_sum=chisp_loc_sum/beta

  !write bosonic index, chi and chi_all to standard output
  WRITE(6,*)sum_ind_ch,chich_loc,chich_loc_all
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  WRITE(6,*)sum_ind_sp,chisp_loc,chisp_loc_all
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

  !end of local part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !nonlocal part

  !allocate nonlocal susceptibilities (bare bubble, q-dependent DMFT susceptibility)
  ALLOCATE(chi_bubble(-Iwbox:Iwbox-1,1:(LQ+1)*LQ/2))   !only qx >= qy
  ALLOCATE(chich_x0(1:(LQ+1)*LQ/2))   !only qx >= qy
  ALLOCATE(chisp_x0(1:(LQ+1)*LQ/2))   !only qx >= qy
  
  qmax=-dacos(1.0d0)   !max_q[chi(q,omega=0)], not yet implemented -> set to pi!!!!

  !allocate klist (contains list of external k-points)
  ALLOCATE(klist(k_number,2))
  !read external k-points (only rank 0)
  IF (myid.EQ.0) THEN
     CALL read_klist('klist.dat',k_number,klist)
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  CALL MPI_BCAST(klist,2*k_number,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)  
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
  !initialize the cos()- and sin()-arrays for the three momenta
  CALL init_arrays(Nint,LQ,k_number,klist,qmax,Q0b,Qv,dcok,dsik,dcoq,dsiq,dcol,dsil)

  !calculate bare susceptibility (bubble term)
  CALL calc_bubble(mu,beta,Iwbox,i,self,LQ,Nint,dcok,dsik,dcoq,dsiq,chi_bubble)

  !here: determination of qmax is possible!!!
  !CALL init_arrays -> with new qmax!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IF(.NOT.sigma_only) THEN
     !calculation of the (not lambda-corrected) chi(x=0)
     !in spin and charge channel
     DO ix=0,LQ-1
        !multiplicity (borders of BZ)
        a=1.d0
        IF (ix.EQ.0.OR.ix.EQ.LQ-1) THEN
           a=a*0.5d0
        ENDIF
        qx = Qv(ix)
        DO iy=0,ix
           b=a
           IF (iy.eq.0.or.iy.eq.LQ-1) then
              b=b*0.5d0
           ENDIF
           qy = Qv(iy)
           
           !multiplicity (fully irreducible BZ)
           b=b*dfloat(2/((1+iy)/(1+ix)+1))
           
           ind=ix*(ix+1)/2+iy+1

           !write status to standard output
           IF (myid.EQ.0) THEN
              WRITE(6,*) 'i,j=',ix,iy,'Qx =',qx,'Qy =',qy,'b',b 
           ENDIF
           
           !calculate chi (without lambda correction)
           chich_x0(ind)=calc_chi(Iwbox,beta,gammach, &
                chi_bubble(:,ind))
           chisp_x0(ind)=calc_chi(Iwbox,beta,gammasp, &
                chi_bubble(:,ind))
           !determine starting value for lambda-correction
           !in order to get positive chi(w=0,q)
           IF (i.EQ.0) THEN
              IF ((1.0d0/dreal(chich_x0(ind))).LT. &
                   (1.0d0/dreal(chich_inv_min))) THEN
                 chich_inv_min=chich_x0(ind)
              ENDIF
              IF ((1.0d0/dreal(chisp_x0(ind))).LT. &
                   (1.0d0/dreal(chisp_inv_min))) THEN
                 chisp_inv_min=chisp_x0(ind)
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
     CALL MPI_BCAST(chich_inv_min,1,MPI_COMPLEX16,Iwbox-1-shift, &
          MPI_COMM_WORLD,ierror)
     CALL MPI_BCAST(chisp_inv_min,1,MPI_COMPLEX16,Iwbox-1-shift, &
          MPI_COMM_WORLD,ierror)
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !calculate lambda corrections
     
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
     lambdach=lambda('lambda_correction_ch.dat',1,myid,i,sum_ind_ch,beta,LQ,chich_loc_sum, &
          chich_inv_min,chich_x0)
     lambdasp=lambda('lambda_correction_sp.dat',2,myid,i,sum_ind_sp,beta,LQ,chisp_loc_sum, &
          chisp_inv_min,chisp_x0)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
     
     !end lambda correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     !create filename for charge susceptibility
     ALLOCATE(CHARACTER(LEN=18)::fname)
     IF ((myid+shift).LT.10) THEN
        WRITE(fname,'(A12,5Hchi00,I1)')'chich_omega/',myid+shift
     ELSEIF ((myid+shift).LT.100) THEN
        WRITE(fname,'(A12,4Hchi0,I2)')'chich_omega/',myid+shift
     ELSE
        WRITE(fname,'(A12,3Hchi,I3)')'chich_omega/',myid+shift
     ENDIF
     !write chich(q,omega) for calculated lambdach and for lambdach
     CALL write_chi(fname,LQ,Qv,lambdach,chich_x0)
     DEALLOCATE(fname)
     !create filename for spin susceptibility
     ALLOCATE(CHARACTER(LEN=18)::fname)
     IF ((myid+shift).LT.10) THEN
        WRITE(fname,'(A12,5Hchi00,I1)')'chisp_omega/',myid+shift
     ELSEIF ((myid+shift).LT.100) THEN
        WRITE(fname,'(A12,4Hchi0,I2)')'chisp_omega/',myid+shift
     ELSE
        WRITE(fname,'(A12,3Hchi,I3)')'chisp_omega/',myid+shift
     ENDIF
     !write chich(q,omega) for calculated lambdach and for lambdach
     CALL write_chi(fname,LQ,Qv,lambdasp,chisp_x0)
     DEALLOCATE(fname)

  ELSE !else from "if (.not.sigma_only)"
  
     !create filename for charge susceptibility
     ALLOCATE(CHARACTER(LEN=18)::fname)
     IF ((myid+shift).LT.10) THEN
        WRITE(fname,'(A12,5Hchi00,I1)')'chich_omega/',myid+shift
     ELSEIF ((myid+shift).LT.100) THEN
        WRITE(fname,'(A12,4Hchi0,I2)')'chich_omega/',myid+shift
     ELSE
        WRITE(fname,'(A12,3Hchi,I3)')'chich_omega/',myid+shift
     ENDIF
     !read chich(q,omega) for lambdach=0
     CALL read_chix0(fname,LQ,chich_x0)
     DEALLOCATE(fname)
     
     !create filename for spin susceptibility
     ALLOCATE(CHARACTER(LEN=18)::fname)
     IF ((myid+shift).LT.10) THEN
        WRITE(fname,'(A12,5Hchi00,I1)')'chisp_omega/',myid+shift
     ELSEIF ((myid+shift).LT.100) THEN
        WRITE(fname,'(A12,4Hchi0,I2)')'chisp_omega/',myid+shift
     ELSE
        WRITE(fname,'(A12,3Hchi,I3)')'chisp_omega/',myid+shift
     ENDIF
     !read chich(q,omega) for lambdach=0
     CALL read_chix0(fname,LQ,chisp_x0)
     DEALLOCATE(fname)

     !set lambda values to the ones read from input file
     lambdach=xch_so
     lambdasp=xsp_so

  ENDIF   !endif from "if (.not.sigma_only)"
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Calculation of the lambda-corrected Self-energy with the lambda-value read or obtained above

  IF(.NOT.chi_only) THEN

     ALLOCATE(selflist(k_number,-Iwbox:Iwbox-1))
     ALLOCATE(selflist_res(k_number,-Iwbox:Iwbox-1))
     ALLOCATE(selflistch(k_number,-Iwbox:Iwbox-1))
     ALLOCATE(selflistch_res(k_number,-Iwbox:Iwbox-1))
     ALLOCATE(selflistsp(k_number,-Iwbox:Iwbox-1))
     ALLOCATE(selflistsp_res(k_number,-Iwbox:Iwbox-1))
     ALLOCATE(selflistrest(k_number,-Iwbox:Iwbox-1))
     ALLOCATE(selflistrest_res(k_number,-Iwbox:Iwbox-1))

     CALL calc_self(Iwbox,myid,i,LQ,k_number,uhub,mu,beta,lambdach,lambdasp, &
       self,fupdown,gammach,gammasp, chi_bubble,chich_x0,chisp_x0, &
       dcoq,dsiq,dcol,dsil, &
       selflist,selflistch,selflistsp,selflistrest)

     !sum over omega
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
     CALL MPI_REDUCE(selflist,selflist_res,2*Iwbox*k_number, &
          MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     CALL MPI_REDUCE(selflistch,selflistch_res,2*Iwbox*k_number, &
          MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     CALL MPI_REDUCE(selflistsp,selflistsp_res,2*Iwbox*k_number, &
          MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     CALL MPI_REDUCE(selflistrest,selflistrest_res,2*Iwbox*k_number, &
          MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
      
     !output of the lambda-corrected self-energy
     IF (myid.EQ.0) THEN
        CALL write_sigma(Iwbox,k_number,beta,selflist_res,selflistch_res, &
             selflistsp_res,selflistrest_res,self,selfloc_res)
     ENDIF

  ENDIF   !endif (.not.chi_only)

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  IF(i .EQ. 0) THEN
     WRITE(6,*),'Finished.'
  ENDIF

  CALL MPI_FINALIZE(ierror)
  
  STOP
  
END PROGRAM self_k