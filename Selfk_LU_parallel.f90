! ladder-dga with automatic lambda-correction
! ver 3.0
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
  ALLOCATE(sendstatus(MPI_STATUS_SIZE))
  ALLOCATE(recvstatus(MPI_STATUS_SIZE))

  !reading parameter from file ladderDGA.in:
  !-) uhub...Hubbard interaction
  !-) mu...chemical potential
  !-) beta...inverse temperature
  !-) nden...average density per lattice site
  !-) Iwbox...number of fermionic Matsubara frequencies
  !-) Iwbox_bose...number of bosonic Matsubara frequencies
  !-) shift...offset for bosonic frequencies
  !-) LQ...interval (0,pi) for internal q-summation is diveded into (LQ-1) parts (LQ=0 -> 0, LQ->pi)
  !-) Nint...interval (-pi,pi) for internal k' summation is split into Nint subintervals
  !-) k_number...number of external k-points for which sigma_dga is calculated
  !-) chi_only...only chi(q,omega) is calculated
  !-) lambdaspin_only..lambda_correction is performed only in the magnetic channel

  !only rank 0 reads the parameters
  IF (myid.EQ.0) THEN
     CALL read_parameters('ladderDGA.in',uhub,mu,beta,nden, &
          Iwbox,Iwbox_bose,shift,LQ,Nint,k_number,chi_only,lambdaspin_only,sumallch,sumallsp)
     !Check parameters
     WRITE(6,*) 'U= ', uhub
     WRITE(6,*) 'MU=',mu
     WRITE(6,*) 'BETA=',beta
     WRITE(6,*) '< n >=', nden 
     WRITE(6,*) 'fermionic frequency box=',Iwbox
     WRITE(6,*) 'bosonic frequency box=',Iwbox_bose
     WRITE(6,*) 'shift of bosonic frequency box=',shift
     WRITE(6,*) "number of q-points=",LQ
     WRITE(6,*) "number of k'-intervals=",Nint
     WRITE(6,*) "number of external k-points=",k_number
     WRITE(6,*) 'chi only=',chi_only
     WRITE(6,*) 'Lambda correction only for the magnetic channel (lambda_charge=0)?',lambdaspin_only
     WRITE(6,*) 'Sum for lambda correction in the density channel over all bosonic frequencies?',sumallch
     WRITE(6,*) 'Sum for lambda correction in the magnetic channel over all bosonic frequencies?',sumallsp
  ENDIF

  !Broadcast parameters to all ranks
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  CALL MPI_BCAST(uhub,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(mu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(beta,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(nden,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(Iwbox,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(Iwbox_bose,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(shift,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(LQ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(Nint,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(k_number,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(chi_only,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(lambdaspin_only,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(sumallch,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(sumallsp,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

  !determine bosonic frequency index i
  !-> should be generalized to an array of bosonic frequenies if #ranks<#bosonic frequencies!!!!!
  i=myid-Iwbox_bose+shift
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
  ALLOCATE(selfloc(0:Iwbox-1))
  CALL calc_self_loc(Iwbox,i,uhub,beta,gww,fupdown,gammach,gammasp, &
       selfloc)
  !sum self(nu,i) over bosonic index i
  ALLOCATE(selfloc_res(0:Iwbox-1))
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  CALL MPI_REDUCE(selfloc,selfloc_res,Iwbox,MPI_COMPLEX16, &
       MPI_SUM,0,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

  IF (myid.EQ.0) THEN
     DO j=0,Iwbox-1
        selfloc_res(j)=selfloc_res(j)*uhub/beta+uhub/2.d0*nden
     ENDDO
  ENDIF

  ! output of local self-energy
  IF (myid.EQ.0) THEN
     CALL write_self_loc('klist/SELF_LOC_parallel',Iwbox,beta,self(0:Iwbox-1),selfloc_res)
  ENDIF

  !determine bosonices indices for which chich_loc_res, chisp_loc_res < 0 (sum_ind_ch, sum_ind_sp)
  IF (i.GE.0) THEN
     ind_part_ch=1.0d0/DSIGN(dfloat(i)+0.5d0,dreal(chich_loc)+tolerancech)
     ind_part_sp=1.0d0/DSIGN(dfloat(i)+0.5d0,dreal(chisp_loc)+tolerancesp)
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
  sum_ind_ch=INT(DABS(1.0d0/ind_ch))+(INT(DSIGN(1.0d0,ind_ch))-1)/2
  sum_ind_sp=INT(DABS(1.0d0/ind_sp))+(INT(DSIGN(1.0d0,ind_sp))-1)/2

  !If lambda-correction should be calculated with all bosonic frequencies sumallch/sumallsp = .TRUE.
  IF (sumallch) THEN
     sum_ind_ch=Iwbox_bose
  ENDIF
  IF (sumallsp) THEN
     sum_ind_sp=Iwbox_bose
  ENDIF
  
  IF (myid.EQ.0) THEN
     WRITE(6,*)'Maximal index for bosonic frequency summation in the charge-channel:',sum_ind_ch
     WRITE(6,*)'Maximal index for bosonic frequency summation in the spin-channel:',sum_ind_sp
  ENDIF

  IF (myid.EQ.0) THEN
     IF (sum_ind_ch.LT.0) THEN
        WRITE(6,*)'No bosonic frequencies for charge-sum'
     ENDIF
     IF (sum_ind_sp.LT.0) THEN
        WRITE(6,*)'No bosonic frequencies for spin-sum'
     ENDIF
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
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

  !end of local part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !nonlocal part

  !allocate nonlocal susceptibilities (bare bubble, q-dependent DMFT susceptibility)
  ALLOCATE(chi_bubble(-Iwbox:Iwbox-1,1:(LQ+2)*(LQ+1)*LQ/6))   !only qx >= qy
  ALLOCATE(chi_bubble_pos(0:Iwbox-1,1:(LQ+2)*(LQ+1)*LQ/6))
  ALLOCATE(chi_bubble_neg(0:Iwbox-1,1:(LQ+2)*(LQ+1)*LQ/6))
  ALLOCATE(chich_x0(1:(LQ+2)*(LQ+1)*LQ/6))   !only qx >= qy
  ALLOCATE(chisp_x0(1:(LQ+2)*(LQ+1)*LQ/6))   !only qx >= qy
  
  !if the self-energy should be also calculated, allocate trilex vertices
  IF (.NOT.chi_only) THEN
     ALLOCATE(trilexch(-Iwbox:Iwbox-1,1:(LQ+2)*(LQ+1)*LQ/6))
     ALLOCATE(trilexsp(-Iwbox:Iwbox-1,1:(LQ+2)*(LQ+1)*LQ/6))
  ENDIF

  !allocate klist (contains list of external k-points)
  ALLOCATE(klist(k_number,3))
  !read external k-points (only rank 0)
  IF (myid.EQ.0) THEN
     CALL read_klist('klist.dat',k_number,klist)
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  CALL MPI_BCAST(klist,3*k_number,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
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
  ALLOCATE(dcol(k_number,3))
  ALLOCATE(dsil(k_number,3))

  !determination of qmax for spin channel (w=0) (maximizes Re(ChiS(w=0,q=(pi,pi,qz))))
  !assuming an ordering vector of (pi,pi,qz)
  IF(i.EQ.0) THEN
     WRITE(6,*)'determination of qmax'
     qmax=find_qmax(mu,beta,Iwbox,self,LQ,Nint,gammasp)
     WRITE(6,*)'qmax=',qmax
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  CALL MPI_BCAST(qmax,1,MPI_REAL8,Iwbox_bose-shift,MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

  CALL init_arrays(Nint,LQ,k_number,klist,0d0,pi,qmax,Q0b,Qv,dcok,dsik,dcoq,dsiq,dcol,dsil)

  !calculate bare susceptibility (bubble term)
  CALL calc_bubble(mu,beta,Iwbox,i,self,LQ,Nint,dcok,dsik,dcoq,dsiq,chi_bubble_pos)
  !copy chi_bubble for nu>0 to nu<0 by using the symmetry chi_bubble(-nu,omega,q)=dcmplx(chi_bubble(nu,-omega,q))
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  request=MPI_REQUEST_NULL
  CALL MPI_ISSEND(chi_bubble_pos,Iwbox*LQ*(LQ+1)*(LQ+2)/6,MPI_COMPLEX16,-myid+2*Iwbox_bose,98,MPI_COMM_WORLD,request,ierror)
  CALL MPI_RECV(chi_bubble_neg,Iwbox*LQ*(LQ+1)*(LQ+2)/6,MPI_COMPLEX16,-myid+2*Iwbox_bose,98,MPI_COMM_WORLD,recvstatus,ierror)
  CALL MPI_WAIT(request,sendstatus,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  DO ind=1,LQ*(LQ+1)*(LQ+2)/6
     DO j=-Iwbox,-1
        chi_bubble(j,ind)=CONJG(chi_bubble_neg(-j-1,ind))
     ENDDO
     DO j=0,Iwbox-1
        chi_bubble(j,ind)=chi_bubble_pos(j,ind)
     ENDDO
  ENDDO
  DEALLOCATE(chi_bubble_pos)
  DEALLOCATE(chi_bubble_neg)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !calculation of the (not lambda-corrected) chi(x=0)
  !in spin and charge channel
  chich_q_sum=dcmplx(0.0d0,0.0d0)
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
        DO iz=0,iy
           c=b
           IF (iz.eq.0.or.iz.eq.LQ-1) then
              c=c*0.5d0
           ENDIF
           qz = Qv(iz)
           
           !multiplicity (fully irreducible BZ)
           c=c*dfloat(6/ &
                ((1+iy)/(1+ix)+(1+iz)/(1+iy)+ &
                3*((1+iz)/(1+ix))+1))

           ind=ix*(ix+1)*(ix+2)/6+iy*(iy+1)/2+iz+1
           
           !write status to standard output
           IF (myid.EQ.0) THEN
              WRITE(6,*) 'i,j,k=',ix,iy,iz,'Qx =',qx,'Qy =',qy,'Qz =',qz,'c',c 
           ENDIF
           
           !if chi_only=.true., calculate only chi (without lambda correction)
           IF (chi_only) THEN
              chich_x0(ind)=calc_chi(Iwbox,beta,gammach, &
                   chi_bubble(:,ind))
              chisp_x0(ind)=calc_chi(Iwbox,beta,gammasp, &
                   chi_bubble(:,ind))
              !else calculate chi (without lambda correction) and trilex
           ELSE
              CALL calc_chi_trilex(Iwbox,uhub,beta,gammach,chi_bubble(:,ind), &
                   chich_x0(ind),trilexch(:,ind))
              CALL calc_chi_trilex(Iwbox,-uhub,beta,gammasp,chi_bubble(:,ind), &
                   chisp_x0(ind),trilexsp(:,ind))
           ENDIF
           

           !Calculate reference value for lambda_correction in the spin-channel only
           IF ((i.GE.-sum_ind_ch).AND.(i.LE.sum_ind_ch)) THEN
              chich_q_sum=chich_q_sum+c*chich_x0(ind)
           ENDIF

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
  ENDDO
  
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  CALL MPI_ALLREDUCE(chich_q_sum,chich_sum,1,MPI_COMPLEX16,MPI_SUM, &
       MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(chich_inv_min,1,MPI_COMPLEX16,Iwbox_bose-shift, &
       MPI_COMM_WORLD,ierror)
  CALL MPI_BCAST(chisp_inv_min,1,MPI_COMPLEX16,Iwbox_bose-shift, &
       MPI_COMM_WORLD,ierror)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
  chich_sum=chich_sum/(beta*dfloat(LQ-1)**3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calculate lambda corrections
  
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
  IF (.NOT.lambdaspin_only) THEN
     lambdach=lambda('lambda_correction_ch.dat',1,myid,i,sum_ind_ch,beta,LQ,chich_loc_sum, &
          chich_inv_min,chich_x0)
     lambdasp=lambda('lambda_correction_sp.dat',2,myid,i,sum_ind_sp,beta,LQ,chisp_loc_sum, &
          chisp_inv_min,chisp_x0)
  ELSE
     lambdach=0.0d0
     lambdasp=lambda('lambda_correction_sp.dat',2,myid,i,sum_ind_sp,beta,LQ,chich_loc_sum+chisp_loc_sum-chich_sum, &
          chisp_inv_min,chisp_x0)
  ENDIF
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
  !write chisp(q,omega) for calculated lambdasp and for lambdasp
  CALL write_chi(fname,LQ,Qv,lambdasp,chisp_x0)
  DEALLOCATE(fname)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  Calculation of the lambda-corrected Self-energy with the lambda-value read or obtained above

  IF(.NOT.chi_only) THEN

     ALLOCATE(selflist(k_number,0:Iwbox-1))
     ALLOCATE(selflist_res(k_number,0:Iwbox-1))
     ALLOCATE(selflistch(k_number,0:Iwbox-1))
     ALLOCATE(selflistch_res(k_number,0:Iwbox-1))
     ALLOCATE(selflistsp(k_number,0:Iwbox-1))
     ALLOCATE(selflistsp_res(k_number,0:Iwbox-1))
     ALLOCATE(selflistrest(k_number,0:Iwbox-1))
     ALLOCATE(selflistrest_res(k_number,0:Iwbox-1))

     CALL calc_self(Iwbox,myid,i,LQ,k_number,uhub,mu,beta,lambdach,lambdasp, &
       self,fupdown,gammach,gammasp, chi_bubble,chich_x0,chisp_x0,trilexch,trilexsp, &
       dcoq,dsiq,dcol,dsil, &
       selflist,selflistch,selflistsp,selflistrest)

     !sum over omega
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
     CALL MPI_REDUCE(selflist,selflist_res,Iwbox*k_number, &
          MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     CALL MPI_REDUCE(selflistch,selflistch_res,Iwbox*k_number, &
          MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     CALL MPI_REDUCE(selflistsp,selflistsp_res,Iwbox*k_number, &
          MPI_COMPLEX16,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     CALL MPI_REDUCE(selflistrest,selflistrest_res,Iwbox*k_number, &
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
