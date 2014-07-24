      PROGRAM VecBosProd 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Input and output strings.
      CHARACTER FRAME*12,BEAM*12,TARGET*12,PARAM*100
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      INTEGER iboson,iferm1,iferm2,iferm3,iferm4,iferm5,ngamma
      INTEGER igamma(20)
      INTEGER idboson,idferm1,idferm2,idferm3,idferm4,idferm5 ! third fermion for tau->enunu
      INTEGER total
      INTEGER i,ii
      INTEGER idebug
      real qsq,shat,that,uhat,x1,x2

      idebug=0

C...Read parameters for PYINIT call.
      READ(*,*) FRAME,BEAM,TARGET,ENERGY

C...Read number of events to generate, and to print.
      READ(*,*) NEV,NPRT

C...Loop over reading and setting parameters/switches.
  100 READ(*,'(A)',END=200) PARAM
      CALL PYGIVE(PARAM)
      GOTO 100

 200  continue

C  open output file
      OPEN(UNIT=18,FILE='bosons.txt',FORM='FORMATTED')
      total = 0

C...Initialize PYTHIA.
      CALL PYINIT(FRAME,BEAM,TARGET,ENERGY)

C...Warmup
      DO 400 IEV=1,1000
        CALL PYEVNT
  400 CONTINUE

C...Event generation loop
      DO 300 IEV=1,NEV
        CALL PYEVNT
C        IF (IEV.LE.NPRT) CALL PYLIST(1)

C       Analyse event
        iboson=-1
        ibosonint=-1
        iferm1=-1
        iferm2=-1
        iferm3=-1
        iferm4=-1
        iferm5=-1
        ngamma=0

        if(mod(IEV, 10000).eq.0) then
           print *, 'Processing event: ', IEV
        endif

        qsq  = PARI(22)
        that = PARI(15)
        uhat = PARI(16)
        x1 = PARI(33)
        x2 = PARI(34)
        flav1 = MSTI(15)
        flav2 = MSTI(16)

        do i=1,n

           if (K(i,2).eq.23) then     ! we have found a Z
              if (idebug.gt.0) print *,'Z: ',i,K(i,3),K(i,1)
              if (K(i,1).eq.11) then  ! this is the final decayed Z
                 iboson=i
                 idboson=90
                 if (K(i,3).eq.0) then
                    print *,'Orphan Z !'
                    goto 300
                 end if
                 if (K(K(i,3),2).ne.23) then
                    print *,'Cannot find intermediate Z entry !'
                    goto 300
                 else
                    ibosonint=K(i,3)
                 end if
                 do ii=i+1,n
                    if (K(ii,3).eq.iboson) then  ! have direct daughter from Z
                       print *,'Z with unexpected daughters !'
                       goto 300
                    end if
                    if (K(ii,3).ne.0) then
                       if (K(K(ii,3),3).eq.ibosonint) then  ! have daughter
                                                            ! from intermediate Z
                          if (idebug.gt.0) 
     *                         print *,'dau: ',K(ii,1),K(ii,2),K(ii,3)
                          if ((abs(K(ii,2)).eq.11) .or. (abs(K(ii,2)).eq.13).or.(abs(K(ii,2)).eq.12).or.
     * (abs(K(ii,2)).eq.14).or.(abs(K(ii,2)).eq.16)) then      ! electron, muon, tau or neutrinos
                             if (iferm1.lt.0) then
                                iferm1=ii
                                if (K(ii,2).eq.-11) then    ! e+
                                   idferm1=-12
                                else if (K(ii,2).eq. 11) then !e-
                                   idferm1=12
                                else if (K(ii,2).eq.-13) then !mu+
                                   idferm1=-14
                                else if (K(ii,2).eq. 13) then !mu-
                                   idferm1=14
                                else if (K(ii,2).eq. -12) then !nu_e_bar
                                   idferm1=-11
                                else if (K(ii,2).eq. 12) then !nu_e
                                   idferm1=11
                                else if (K(ii,2).eq. -14) then !nu_mu_bar
                                   idferm1=-13
                                else if (K(ii,2).eq. 14) then !nu_mu
                                   idferm1=13
                                else if (K(ii,2).eq. -16) then !nu_tau_bar
                                   idferm1=-15
                                else if (K(ii,2).eq. 16) then !nu_tau
                                   idferm1=15
                                end if
                             else if (iferm2.lt.0) then
                                iferm2=ii
                                if (K(ii,2).eq.-11) then    ! e+
                                   idferm2=-12
                                else if (K(ii,2).eq. 11) then !e-
                                   idferm2=12
                                else if (K(ii,2).eq.-13) then !mu+
                                   idferm2=-14
                                else if (K(ii,2).eq. 13) then !mu-
                                   idferm2=14
                                else if (K(ii,2).eq. -12) then !nu_e_bar
                                   idferm2=-11
                                else if (K(ii,2).eq. 12) then !nu_e
                                   idferm2=11
                                else if (K(ii,2).eq. -14) then !nu_mu_bar
                                   idferm2=-13
                                else if (K(ii,2).eq. 14) then !nu_mu
                                   idferm2=13
                                else if (K(ii,2).eq. -16) then !nu_tau_bar
                                   idferm2=-15
                                else if (K(ii,2).eq. 16) then !nu_tau
                                   idferm2=15
                                end if
                             else
                                print *,'More than two fermions !'
                                goto 300
                             end if
                          else if (abs(K(ii,2)).eq.22) then  ! photon
                             ngamma=ngamma+1
                             igamma(ngamma)=ii
                          else if (idebug.gt.0) then
                             print *,'Unexpected daughter: ',abs(K(ii,2))
                          end if
                       end if
                    end if
                 end do
              end if
           end if

           if (abs(K(i,2)).eq.24) then ! we have found a W
              if (idebug.gt.0) print *,'W: ',i,K(i,3),K(i,1)
              if (K(i,1).eq.11) then  ! this is the final decayed W
                 iboson=i
                 if (K(i,2).eq.24) then ! W+
                    idboson=80
                 else
                    idboson=-80
                 end if
                 if (K(i,3).eq.0) then
                    print *,'Orphan W !'
                    goto 300
                 end if
                 if (abs(K(K(i,3),2)).ne.24) then
                    print *,'Cannot find intermediate W entry !'
                    goto 300
                 else
                    ibosonint=K(i,3)
                 end if
                 do ii=i+1,n
                    if (K(ii,3).eq.iboson) then  ! have direct daughter from W
                       print *,'W with unexpected daughters !'
                       goto 300
                    end if
                    if (K(ii,3).ne.0) then
                       if (K(K(ii,3),3).eq.ibosonint) then  ! have daughter
                                                            ! from intermediate W
                          if (idebug.gt.0) 
     *                         print *,'dau: ',K(ii,1),K(ii,2),K(ii,3)
C DB : this is the case for W->enu
                          if (abs(K(ii,2)).eq.11) then      ! electron
                             if (iferm1.lt.0) then
                                iferm1=ii
                                if (K(ii,2).eq.-11) then    ! e+
                                   idferm1=-12
                                else
                                   idferm1=12
                                end if
                             else if (iferm2.lt.0) then
                                iferm2=ii
                                if (K(ii,2).eq.-11) then    ! e+
                                   idferm2=-12
                                else
                                   idferm2=12
                                end if
                             else
                                print *,'More than two fermions !'
                                goto 300
                             end if
                          else if (abs(K(ii,2)).eq.12) then ! neutrino
                             if (iferm1.lt.0) then
                                iferm1=ii
                                if (K(ii,2).eq.-12) then    ! nu_e_bar
                                   idferm1=-11
                                else
                                   idferm1=11
                                end if
                             else if (iferm2.lt.0) then
                                iferm2=ii
                                if (K(ii,2).eq.-12) then    ! nu_e_bar
                                   idferm2=-11
                                else
                                   idferm2=11
                                end if
                             else
                                print *,'More than two fermions !'
                                goto 300
                             end if
C DB : this is the case for w->taunu, put nu_tau first
                          else if ((abs(K(ii,2)).eq.16).and.(abs(K(K(K(ii,3),3),2)).eq.24)) then
                             if (idebug.gt.0) then
                                print *,'nu_tau debug',K(ii,2),K(ii,3),K(K(ii,3),2),K(K(K(ii,3),3),2)
                                print *,'nu_tau daught pt',P(ii,1),P(ii,2),P(ii,3),P(ii,4)
                             end if
                             if(iferm1.lt.0) then
                                iferm1=ii
                                if(K(ii,2).eq.-16) then
                                   idferm1=-15
                                else if(K(ii,2).eq.16) then
                                   idferm1=15
                                end if
                             end if
                          else if (abs(K(ii,2)).eq.15) then  ! tau
			     if (idebug.gt.0) print *,'tau debug',K(ii,3),K(K(ii,3),2)
			     do jj=ii+1,n
				if (abs(K(K(jj,3),2)).eq.15) then
				 if (idebug.gt.0) then
				    print *,'tau daught debug',K(jj,2),K(jj,3),K(K(jj,3),2)
				    print *,'tau daught pt',P(jj,1),P(jj,2),P(jj,3),P(jj,4)
				 end if
C DB : put nu_tau first so we can handle both tau->e nu_en u_tau AND tau->pi nu_tau
                                 if (abs(K(jj,2)).eq.16) then ! neutrino
                                    if(iferm2.lt.0) then
                                       iferm2=jj
                                       if(K(jj,2).eq.-16) then
                                          idferm2=-15
                                       else if(K(jj,2).eq.16) then
                                          idferm2=15
                                       end if
                                    end if
                                 else if (abs(K(jj,2)).eq.211) then ! pion
                                    if(iferm3.lt.0) then
                                       iferm3=jj
                                       if(K(jj,2).eq.-211) then
                                          idferm3=-120
                                       else if(K(jj,2).eq.211) then
                                          idferm3=120
                                       end if
                                    end if
				 else if (abs(K(jj,2)).eq.11) then      ! electron
				    if (iferm3.lt.0) then
				       iferm3=jj
				       if (K(jj,2).eq.-11) then    ! e+
					  idferm3=-12
				       else
					  idferm3=12
				       end if
				    else if (iferm4.lt.0) then
				       iferm4=jj
				       if (K(jj,2).eq.-11) then    ! e+
					  idferm4=-12
				       else
					  idferm4=12
				       end if
				    else
				       print *,'More than two fermions tau !'
				       goto 300
				    end if
				 else if (abs(K(jj,2)).eq.12) then ! neutrino
				    if (iferm3.lt.0) then
				       iferm3=jj
				       if (K(jj,2).eq.-12) then    ! nu_e_bar
					  idferm3=-11
				       else
					  idferm3=11
				       end if
				    else if (iferm4.lt.0) then
				       iferm4=jj
				       if (K(jj,2).eq.-12) then    ! nu_e_bar
					  idferm4=-11
				       else
					  idferm4=11
				       end if
				    else
				       print *,'More than two fermions nue !'
				       goto 300
				    end if
				 end if 
				end if
			     end do
                             if (abs(K(ii,2)).eq.15) then ! tau
                               if(iferm5.lt.0) then
                                   iferm5=ii
                                   if(K(ii,2).eq.-15) then
                                     idferm5=-16
                                   else if(K(ii,2).eq.15) then
                                     idferm5=16
                                   end if
                               end if
                             end if		     
                          else if (abs(K(ii,2)).eq.22) then  ! photon
                             ngamma=ngamma+1
                             igamma(ngamma)=ii
                          else if (idebug.gt.0) then
                             print *,'Unexpected daughter: ',abs(K(ii,2))
                          end if
                       end if
                    end if
                 end do
              end if
           end if

        end do

c       write out event information
	if (idebug.gt.0) print *,'debug output',iboson,iferm1,iferm2
        if ((iboson.ge.0).and.(iferm1.ge.0).and.(iferm2.ge.0)) then
           total=total+1
c          event number and weight
c          write(18, *) total, 1.
           write(18, *) total,1.,qsq,that,uhat,x1,x2,flav1,flav2 
c          vertex
           write(18, *) 0., 0., 0.
c          boson
           write(18, *) 
     *          idboson, 
     *          P(iboson,1), P(iboson,2), P(iboson,3), P(iboson,4), 0, 0
c          fermions
           write(18, *) 
     *          idferm1, P(iferm1,1), P(iferm1,2),
     *          P(iferm1,3), P(iferm1,4), 1, 0
           write(18, *) 
     *          idferm2, P(iferm2,1), P(iferm2,2),
     *          P(iferm2,3), P(iferm2,4), 1, 0
	   if (iferm3.gt.0) then
	        write(18, *) 
     *           idferm3, P(iferm3,1), P(iferm3,2),
     *           P(iferm3,3), P(iferm3,4), 1, 0
	   end if
	   if (iferm4.gt.0) then
	        write(18, *) 
     *           idferm4, P(iferm4,1), P(iferm4,2),
     *           P(iferm4,3), P(iferm4,4), 1, 0
	   end if
           if (iferm5.gt.0) then
                write(18, *) 
     *           idferm5, P(iferm5,1), P(iferm5,2),
     *           P(iferm5,3), P(iferm5,4), 0, 0
           end if
c          gammas
           do ii=1,ngamma
              write(18, *) 
     *             10, P(igamma(ii),1), P(igamma(ii),2),
     *             P(igamma(ii),3), P(igamma(ii),4), 1, 0
           end do
c          done with particles
           write(18,*) 0
        end if

 300  CONTINUE

C   write two 0s to indicate the end of the file
      WRITE(18, *) 0, 0

C...Print cross sections.
      CALL PYSTAT(1)

      END

