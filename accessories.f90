module accessories
    
    		
   IMPLICIT NONE

	CONTAINS

    
   

   
    
    subroutine preread_input(file_name,gr,dn,en,n,r)
        
        implicit none
        
        !Externals
        character(30), intent(in):: file_name       !Name of the input file
        integer, intent(out):: dn                   !Nr of subpopulations (demes)
        integer, intent(out):: gr                   !Nr of sampling groups
        integer, intent(out):: en                   !Nr of demographic events
        integer, intent(out):: n                    !Overall sample size (over all samples) 
        integer, intent(out):: r                    !Nr of priors    
        !Internals
        integer:: ierror 

        !---------------------------------------------------------------------------------
        OPEN(UNIT=1, FILE=file_name, STATUS='OLD', ACTION='READ', IOSTAT=ierror) ! Opening the file
	      openif: IF (ierror==0) THEN                                              
	         READ(1,*,IOSTAT=ierror) gr, dn, en, n, r
          ELSE
             CALL ErrorMessage(1005)
          ENDIF openif               
        CLOSE(UNIT=1)
     
        return
     
	endsubroutine preread_input 
	
	
	
	
	
	subroutine read_info(file_name,dn,gr,en,n,Ne,growth,bid,samplinfo,events,gentime,sexratio,marker,labels,r,PrInfo,headings,empCoal,chainLength,int1,seed1)
        
        
        implicit none
        
        !Externals
        character(30), intent(in):: file_name                           !Name of the input file
        integer, intent(in):: dn                                        != # of subpopulations (demes)
        integer, intent(in):: gr                                        != # of sampling groups
        integer, intent(in):: en                                        != # of demographic events
        integer, intent(in):: n                                         != total samlpe size (over all samples) 
        real(8), intent(out),dimension(dn):: Ne                         != vector of population sizes (one for each deme)  
        real(8), intent(out), dimension(dn):: growth                    != growth rate of the corresponding population (deme)
        real(8), intent(out), dimension(dn):: bid                       != identifier of population number
        real(8), intent(out), dimension(3,gr):: samplinfo               !An array with: 1st col:= subsample sizes; 2nd col:= deme; and 3rd col:= sampling time (if ancient, 0 if not)  
        real(8), intent(out), dimension(en,8):: events                  !Array with the info of ancient demographic events: 
                                                                        !(1) time;
                                                                        !(2) Ne; 
                                                                        !(3) Growth rate; 
                                                                        !(4) 1st deme involved; 
                                                                        !(5) 2nd deme involved; 
                                                                        !(6) % of migrating lineages from 1st deme; 
                                                                        !(7) % of migrating lineages from 2nd deme; 
                                                                        !(8) Block id    
        real, intent(out) :: gentime                                    !Generation time(number of years per generation) 
        real, intent(out):: sexratio                                    !Sexratio
        integer, intent (out):: marker                                  !Type of marker: 1)Autosomal haplo; 2)Autosomal diplo; 3)Haplo-diploid/X-linked; 4)Mitochondrial/Y-linked        
        character(30), dimension(n), intent(out):: labels               !Labels of the samples, i.e. the names of the real sequences they represent in the simulations 
        integer, intent(in):: r                                         !# of priors
        real(8), dimension(r,9), intent(out):: PrInfo                   !Array with the info of the priors  
        character(18), dimension(r+2), intent(out):: headings           !Headings of the results file
        real(8), dimension(n-1), intent(out):: empCoal                  !Array with the "observed" coalescent times
		integer, intent(out):: chainLength                              !Length of the MCMC chain
        integer, intent(out):: int1                                     !Interval for writing out the results
		integer, intent(out):: seed1                                    !Seed of pseudo-random number generator  
        
        !Internals        
	    integer:: dn2, gr2, en2, n2, r2                                 !Same as dn, gr, en, n, r, but obtained from input file for verification
        character(12), dimension(3):: charA                             !For capturing the info of the sampling and demes section
        character(12), dimension(8):: charB                             !For capturing the info of the events section
        character(12), dimension(2):: char2                             !For capturing pairs of values
        character(12), dimension(4):: char4                             !For capturing four values
        character(12), dimension(6):: char6                             !For capturing six values
        character(12):: char1                                           !For capturing a value
        integer:: ierror 
	    integer:: i, j, k, l, x, z                                      !Counters                                                 
        character(12), dimension(8):: charPrInfo                        !For reading the priors info  
        character(6):: xstring
        integer:: dummy
        !-----------------------------------------------------------------------------------------------------------------------
	    !  S  T  A  R  T
        !-----------------------------------------------------------------------------------------------------------------------
        headings(1)='Simulation'
        headings(r+2)='L'
        z=0
        !-------------------------------------------------------------------------------------------------------------
        OPEN(UNIT=1, FILE=file_name, STATUS='OLD', ACTION='READ', IOSTAT=ierror) ! Opening the file
	    !-------------------------------------------------------------------------------------------------------------
        openif: IF (ierror==0) THEN                                              ! ierror=0 means Open was succesfull
	  	    WRITE(*,*)'Succesfull!'                                              
		    !---------------------------------------------------------------------------------------------------------
            !(1) Reading the number of demes, sampling groups and events
            READ(1,*,IOSTAT=ierror) gr2, dn2, en2, n2, r2 
            WRITE(*,1000) gr2
		    1000 FORMAT(I5,' sampling groups')
            WRITE(*,1001) dn2
		    1001 FORMAT(I5,' demes sampled')
            WRITE(*,1002) en2
            1002 FORMAT(I5,' events')  
            WRITE(*,1003) n2
            1003 FORMAT(I7,' individuals') 
            WRITE(*,1004) r2
            1004 FORMAT(I5,' priors')  
            
            !---------------------------------------------------------------------------------------------------------
            !(2) Reading sampling values (SAMPLINFO)   
            DO i=1,gr,1
               READ(1,*,IOSTAT=ierror)  (charA(l), l=1,3) 
               READ(charA(1),*) samplinfo(1,i)
               READ(charA(2),*) samplinfo(2,i)
               IF ((TRIM(charA(3))=='Prior').OR.(TRIM(charA(3))=='prior').OR.(TRIM(charA(3))=='PRIOR')) THEN
                  samplinfo(3,i)=-6666666
                  !--
                  z=z+1
                  WRITE( char1, '(I12)' )  i
                  headings(z)='Age_'//TRIM(ADJUSTL(char1))
               ELSE  
                  READ(charA(3),*) samplinfo(3,i) 
               ENDIF 
               IF (ierror/=0) STOP
            ENDDO
            WRITE(*,1020) samplinfo 
            1020 FORMAT(' ',F4.0,' ind from deme ',F4.0,' w/age ', F7.1)
            IF (n2 /= SUM(samplinfo(1,:))) THEN	
               CALL ErrorMessage(1007) 
            ENDIF    
            IF (MAXVAL(samplinfo(2,:))/=dn) THEN
               CALL ErrorMessage(1016) 
            ENDIF
            
            !---------------------------------------------------------------------------------------------------------           
            !(3) Reading population sizes and growth rates
            DO i=1,dn,1
               READ(1,*,IOSTAT=ierror)  (charA(l), l=1,3) 
               IF ((TRIM(charA(1))=='Prior').OR.(TRIM(charA(1))=='prior').OR.(TRIM(charA(1))=='PRIOR')) THEN
                  Ne(i)=-6666666
                  !--
                  z=z+1
                  WRITE( char1, '(I12)' )  i
                  headings(z)='Ne_'//TRIM(ADJUSTL(char1))
               ELSE  
                  READ(charA(1),*) Ne(i) 
               ENDIF 
               IF ((TRIM(charA(2))=='Prior').OR.(TRIM(charA(2))=='prior').OR.(TRIM(charA(2))=='PRIOR')) THEN
                  growth(i)=-6666666
                  !--
                  z=z+1
                  WRITE( char1, '(I12)' )  i
                  headings(z)='Growth_'//TRIM(char1)
               ELSEIF ((TRIM(charA(2))=='Match').OR.(TRIM(charA(2))=='match').OR.(TRIM(charA(2))=='MATCH')) THEN
                  growth(i)=-7777777
               ELSE    
                  READ(charA(2),*) growth(i) 
               ENDIF  
               READ(charA(3),*) bid(i) 
               IF (ierror/=0) STOP
               !--
               WRITE(*,1010) Ne(i), growth(i), bid(i) 
               1010 FORMAT('Ne is',I10,' with growth rate of ',F11.6, ' in pop # ', F3.0) 
            ENDDO
            IF (MAXVAL(bid)/=dn) THEN
               CALL ErrorMessage(1017)
            ENDIF 
            DO i=1,dn,1
               IF ((Ne(i)<1).AND.(Ne(i)/=-6666666)) THEN
                  CALL ErrorMessage(1018)
               ENDIF
            ENDDO 
            
            !---------------------------------------------------------------------------------------------------------            
            !(4) Reading events
            IF (en>0) THEN  
               DO i=1,en,1 
                  READ(1,*,IOSTAT=ierror)  (charB(l), l=1,8)  
                  !1st column, time 
                  IF ((TRIM(charB(1))=='Prior').OR.(TRIM(charB(1))=='prior').OR.(TRIM(charB(1))=='PRIOR')) THEN   
                     events(i,1)=-6666666
                     !--
                     z=z+1
                     WRITE( char1, '(I12)' )  i
                     headings(z)='Ev_time_'//TRIM(ADJUSTL(char1))
                  ELSE  
                     READ(charB(1),*) events(i,1) 
                  ENDIF  
                  !2nd column, Ne  
                  IF ((TRIM(charB(2))=='Prior').OR.(TRIM(charB(2))=='prior').OR.(TRIM(charB(2))=='PRIOR')) THEN   
                     events(i,2)=-6666666
                     !--
                     z=z+1
                     WRITE( char1, '(I12)' )  i
                     headings(z)='Ev_Ne_'//TRIM(ADJUSTL(char1))
                  ELSEIF ((TRIM(charB(2))=='Match').OR.(TRIM(charB(2))=='match').OR.(TRIM(charB(2))=='MATCH')) THEN   
                     events(i,2)=-7777777
                  ELSE
                     READ(charB(2),*) events(i,2) 
                  ENDIF  
                  !3rd column, growth  
                  IF ((TRIM(charB(3))=='Prior').OR.(TRIM(charB(3))=='prior').OR.(TRIM(charB(3))=='PRIOR')) THEN   
                     events(i,3)=-6666666
                     !--
                     z=z+1
                     WRITE( char1, '(I12)' )  i
                     headings(z)='Ev_growth_'//TRIM(char1)
                  ELSEIF ((TRIM(charB(3))=='Match').OR.(TRIM(charB(3))=='match').OR.(TRIM(charB(3))=='MATCH')) THEN   
                     events(i,3)=-7777777
                  ELSE
                     READ(charB(3),*) events(i,3) 
                  ENDIF  
                  !4th column, 1st deme
                  READ(charB(4),*) events(i,4)  
                  !5th column, 2st deme
                  READ(charB(5),*) events(i,5)                  
                  !6th column, 1st deme contribution  
                  IF ((TRIM(charB(6))=='Prior').OR.(TRIM(charB(6))=='prior').OR.(TRIM(charB(6))=='PRIOR')) THEN   
                     events(i,6)=-6666666
                     !--
                     z=z+1
                     WRITE( char1, '(I12)' )  i
                     headings(z)='1st_prop_lin_'//TRIM(ADJUSTL(char1))
                  ELSE
                     READ(charB(6),*) events(i,6) 
                  ENDIF  
                  !7th column, 2nd deme contribution  
                  IF ((TRIM(charB(7))=='Prior').OR.(TRIM(charB(7))=='prior').OR.(TRIM(charB(7))=='PRIOR')) THEN   
                     events(i,7)=-6666666
                     !--
                     z=z+1
                     WRITE( char1, '(I12)' )  i
                     headings(z)='2nd_prop_lin_'//TRIM(ADJUSTL(char1))
                  ELSE
                     READ(charB(7),*) events(i,7) 
                  ENDIF
                  !8th column, Block id
                  READ(charB(8),*) events(i,8) 
                  !------------------------------- 
                  WRITE(*,*)'Demographic events'
                  WRITE(*,1049) i
                  WRITE(*,*) '___________________________________'
                  WRITE(*,1050) events(i,1)
                  WRITE(*,1051) events(i,2)
                  WRITE(*,1052) events(i,3)
                  if (events(i,4)==events(i,5)) then
                    WRITE(*,1053) events(i,4)
                    WRITE(*,1056) events(i,6)*100
                  else
                    WRITE(*,1054) events(i,4)  
                    WRITE(*,1055) events(i,5)
                    WRITE(*,1057) events(i,6)*100
                    WRITE(*,1058) events(i,7)*100
                  endif  
                  WRITE(*,1060) events(i,8)                  
                  WRITE(*,*) '___________________________________'
                  1049 FORMAT('Event number ',I2)    
                  1050 FORMAT('Age =',F12.1)
                  1051 FORMAT('New Ne =',F10.0)
                  1052 FORMAT('New growth rate =',F12.9)
                  1053 FORMAT('Deme involved =',F10.0)
                  1054 FORMAT('1st deme involved =',F4.0)
                  1055 FORMAT('2nd deme involved =',F4.0)
                  1056 FORMAT('Taking',F9.2,'% of its lineages')
                  1057 FORMAT('Taking',F9.2,'% of lineages from 1st deme')
                  1058 FORMAT('Taking',F9.2,'% of lineages from 2nd deme')                   
                  1060 FORMAT('Number of the coalescing block =',F4.0) 
		       ENDDO 	 
            ENDIF 
            DO i=1,en,1
               IF ((events(i,1)<1).AND.(events(i,1)/=-7777777).AND.(events(i,1)/=-6666666)) THEN
                  CALL ErrorMessage(1019)
               ENDIF
               IF ((events(i,2)<1).AND.(events(i,2)/=-7777777).AND.(events(i,2)/=-6666666)) THEN
                  CALL ErrorMessage(1020)
               ENDIF
               IF ((events(i,4)<1).OR.(events(i,4)>dn+en)) THEN
                  CALL ErrorMessage(1021)
               ENDIF
               IF ((events(i,5)<1).OR.(events(i,5)>dn+en)) THEN
                  CALL ErrorMessage(1021)
               ENDIF
               IF (((events(i,6)<0.0).OR.(events(i,6)>1.0)).AND.(events(i,6)/=-7777777).AND.(events(i,6)/=-6666666)) THEN
                  CALL ErrorMessage(1022)
               ENDIF
               IF (((events(i,7)<0.0).OR.(events(i,7)>1.0)).AND.(events(i,7)/=-7777777).AND.(events(i,7)/=-6666666)) THEN
                  CALL ErrorMessage(1022)
               ENDIF
            ENDDO
            IF ( (ANY(events(:,8)>dn+en)).OR.(ANY(events(:,8)<dn)) ) THEN
               CALL ErrorMessage(1023)
            ENDIF
            DO i=1,en,1
               IF ( (events(i,4)==events(i,8)) .OR. (events(i,5)==events(i,8)) ) THEN
                  CALL ErrorMessage(1025)
               ENDIF
            ENDDO
            
            !---------------------------------------------------------------------------------------------------------
            !(6) Generation time, marker type, sex ratio, and type of marker
            !Generation time
            READ(1,*,IOSTAT=ierror) char1              
            IF ((TRIM(char1)=='Prior').OR.(TRIM(char1)=='prior').OR.(TRIM(char1)=='PRIOR')) THEN   
               gentime=-6666666
               !--
               z=z+1
               headings(z)='gen_time'
            ELSE  
               READ(char1,*) gentime 
            ENDIF             
            !Sex ratio (focused on the molecular marker)
            READ(1,*,IOSTAT=ierror) char1              
            IF ((TRIM(char1)=='Prior').OR.(TRIM(char1)=='prior').OR.(TRIM(char1)=='PRIOR')) THEN   
               sexratio=-6666666
               !--
               z=z+1
               headings(z)='sex_ratio'
            ELSE  
               READ(char1,*) sexratio 
            ENDIF 
            !Type of molecular marker
            READ(1,*,IOSTAT=ierror) char1              
            IF ((TRIM(char1)=='Prior').OR.(TRIM(char1)=='prior').OR.(TRIM(char1)=='PRIOR')) THEN   
               marker=-6666666
               !--
               z=z+1
               headings(z)='marker'
            ELSE  
               READ(char1,*) marker 
			ENDIF 
            
            !---------------------------------------------------------------------------------------------------------
            !(7) Reading coalescent times and labels
			!Read the empirical coalescent times
            READ(1,*,IOSTAT=ierror) (empCoal(i), i=1,n-1) 
			!Reading the labels of individuals
            READ(1,*,IOSTAT=ierror) (labels(i), i=1,n)
            
            !---------------------------------------------------------------------------------------------------------
            !(8) Reading Priors
            readloop2: DO i=1,r,1
               READ(1,*,IOSTAT=ierror) charPrInfo(1),charPrInfo(2),charPrInfo(3),charPrInfo(4),charPrInfo(5),charPrInfo(6),charPrInfo(7),charPrInfo(8),PrInfo(i,9) 
               !(1) Retrieving the code of the PDF------------
               IF ((TRIM(charPrInfo(1))=='UNIFORM').OR.(TRIM(charPrInfo(1))=='Uniform').OR.(TRIM(charPrInfo(1))=='uniform')) THEN
                  PrInfo(i,1)=1 
               ELSEIF ((TRIM(charPrInfo(1))=='LOG-UNIFORM').OR.(TRIM(charPrInfo(1))=='Log-Uniform').OR.(TRIM(charPrInfo(1))=='Log-uniform').OR.(TRIM(charPrInfo(1))=='log-Uniform').OR.(TRIM(charPrInfo(1))=='log-uniform').OR.(TRIM(charPrInfo(1))=='LOGUNIFORM').OR.(TRIM(charPrInfo(1))=='LogUniform').OR.(TRIM(charPrInfo(1))=='Loguniform').OR.(TRIM(charPrInfo(1))=='logUniform').OR.(TRIM(charPrInfo(1))=='loguniform')) THEN
                  PrInfo(i,1)=2        
               ELSEIF ((TRIM(charPrInfo(1))=='NORMAL').or.(TRIM(charPrInfo(1))=='Normal').or.(TRIM(charPrInfo(1))=='normal')) THEN 
                  PrInfo(i,1)=3
               ELSEIF ((TRIM(charPrInfo(1))=='LOG-NORMAL').OR.(TRIM(charPrInfo(1))=='Log-Normal').OR.(TRIM(charPrInfo(1))=='Log-normal').OR.(TRIM(charPrInfo(1))=='log-Normal').OR.(TRIM(charPrInfo(1))=='log-normal').OR.(TRIM(charPrInfo(1))=='LOGNORMAL').OR.(TRIM(charPrInfo(1))=='LogNormal').OR.(TRIM(charPrInfo(1))=='Lognormal').OR.(TRIM(charPrInfo(1))=='logNormal').OR.(TRIM(charPrInfo(1))=='lognormal')) THEN
                  PrInfo(i,1)=4
               ELSEIF ((TRIM(charPrInfo(1))=='EXPONENTIAL').or.(TRIM(charPrInfo(1))=='Exponential').or.(TRIM(charPrInfo(1))=='exponential')) THEN  
                  PrInfo(i,1)=5 
               ELSEIF ((TRIM(charPrInfo(1))=='GAMMA').or.(TRIM(charPrInfo(1))=='Gamma').or.(TRIM(charPrInfo(1))=='gamma')) THEN  
                  PrInfo(i,1)=6
               ELSEIF ((TRIM(charPrInfo(1))=='GEOMETRIC').or.(TRIM(charPrInfo(1))=='Geometric').or.(TRIM(charPrInfo(1))=='geometric')) THEN  
                  PrInfo(i,1)=7
               ELSEIF ((TRIM(charPrInfo(1))=='BINOMIAL').or.(TRIM(charPrInfo(1))=='Binomial').or.(TRIM(charPrInfo(1))=='binomial')) THEN  
                  PrInfo(i,1)=8 
               ELSEIF ((TRIM(charPrInfo(1))=='BETA').or.(TRIM(charPrInfo(1))=='Beta').or.(TRIM(charPrInfo(1))=='beta')) THEN  
                  PrInfo(i,1)=9
               ELSEIF ((TRIM(charPrInfo(1))=='POISSON').or.(TRIM(charPrInfo(1))=='Poisson').or.(TRIM(charPrInfo(1))=='poisson')) THEN  
                  PrInfo(i,1)=10
               ELSE
                  CALL ErrorMessage(1006) 
               ENDIF 
               !(2) Retrieving the sign-----------------------
               IF (TRIM(charPrinfo(2))=='-') THEN
                  PrInfo(i,2)=-1.0
               ELSEIF (TRIM(charPrinfo(2))=='+') THEN
                  PrInfo(i,2)=1.0 
               ELSE
                  CALL ErrorMessage(1015) 
               ENDIF   
               !(3-7) Retrieving the parameters---------------
               DO j=3,7,1 
                  IF ((charPrInfo(j)(1:5)=='Prior').or.(charPrInfo(j)(1:5)=='prior').or.(charPrInfo(j)(1:5)=='PRIOR')) THEN
                     READ(charPrInfo(j)(6:12),*) x 
                     PrInfo(i,j)=-5555555000-x   
                  ELSE
                     READ(charPrInfo(j),*) PrInfo(i,j)  
                  ENDIF
               ENDDO               
               
               !(8) Retrieving the type of prior--------------
               IF (TRIM(charPrInfo(8))=='NOVAR') THEN
                  PrInfo(i,8)=0.0
               ELSEIF (TRIM(charPrInfo(8))=='PAR') THEN   
                  PrInfo(i,8)=1.0
               ELSE
                  CALL ErrorMessage(10061) 
               ENDIF  
            ENDDO readloop2
            !---------------------------------------------------------------------------------------------------------   
            !(9) MCMC numbers
            READ(1,*,IOSTAT=ierror) chainLength  
            READ(1,*,IOSTAT=ierror) int1  			
			READ(1,*,IOSTAT=ierror) seed1  
            
        ELSE
		    CALL ErrorMessage(1005)           
        ENDIF openif
        
        CLOSE(UNIT=1)
        
        
        RETURN  
           
	endsubroutine read_info
    
	
	
	
	
    SUBROUTINE get_priors(num,r,gr,samplinfo,dn,Ne,growth,en,events,gentime,sexratio,PrInfo,refPriors,Priors)
   
       implicit none
   
       !Externals 
       integer, intent(in):: num                                      !Nr of simulation: when num=1 priors are sampled directly, otherwise they[re sampled from the jump distr. 
       integer, intent(in):: r                                        !# of priors
       integer, intent(in):: gr                                       !# of sampling groups 
       real(8), dimension(3,gr), intent(inout):: samplinfo            !Age, deme and # of sampling groups
       integer, intent(in):: dn                                       !# of demes
       real(8), dimension(dn), intent(inout):: Ne                     !Ne of all the demes
       real(8), dimension(dn), intent(inout):: growth                 !Growth rate of all the demes 
       integer, intent(in):: en                                       !# of events/blocks
       real(8), dimension(en,8), intent(inout):: events               !Events/blocks info
       real, intent(inout) :: gentime                                 !Generation time
       real, intent(inout):: sexratio                                 !Sex ratio focused on the marker (e.g. 0.25 for mtDNA means 25% females)
       real(8), dimension(r,9), intent(inout):: PrInfo                !Array with the info of the priors
       real(8), dimension(r), intent(inout):: refPriors               !Array with the random values simulated from previous priors (they provide the mean of jump density in mcmc)
       real(8), dimension(r), intent(out):: Priors                    !Array with the random values simulated from the priors
       !Internals
       integer:: i, j, k, l, z, q, a
       real(8):: x                                                    !X gets the number of cross-referred prior  
       logical, dimension(r):: done                                   !Indicator of a prior already sampled (if needed by another one) 
       logical, dimension(5):: y                                      !'y' register which of the 5 parameters of each prior are ready
       logical:: mcmc                                                 !Indicates if sampling directly from priors or from the jump distribution 
       !--------------------------------------------------------------------------------------------------
       !--------------------------------------------------------------------------------------------------
       !We get the priors values first 
       !---------------------------------------------------------------------
       if (num==1) then
          mcmc=.false.
       else
          mcmc=.true.
       endif 
       !--------------
       q=r
       done=.false.
       do
          do i=1,r,1
             if (.NOT.done(i)) then
                if ((ANY(INT(PrInfo(i,:)/10000)==-5555555)).or.(ANY(INT(PrInfo(i,:)/1000)==-5555555))&
                    &.or.(ANY(INT(PrInfo(i,:)/100)==-5555555)).or.(ANY(INT(PrInfo(i,:)/10)==-5555555))) then                 
                    !----This is only for getting the bloody reference of prior (number)
                    do j=3,7,1
                       if  ((INT(PrInfo(i,j)/10)==(-5555555)).or.(INT(PrInfo(i,j)/100)==(-5555555))&
                           &.or.(INT(PrInfo(i,j)/1000)==(-5555555)).or.(INT(PrInfo(i,j)/10000)==(-5555555))) then 
                          a=log10(ABS(PrInfo(i,j)))
                          a=int(a-6) !a stores now the number of cifers in the prior reference
                          x=real(int(PrInfo(i,j)/(10**a)))
                          x=(10**a)*x
                          x=ABS(PrInfo(i,j))+x
                          if (done(int(x))) then
                             PrInfo(i,j)=Priors(int(x)) 
                             y(j-2)=.true. 
                          else
                             y(j-2)=.false. 
                          endif
                       else
                          y(j-2)=.true. 
                       endif    
                    enddo
                    if (ALL(y)) then  !i.e. if the 5 possible parameters required by i-th prior are ready
                       Priors(i)=xrandom(PrInfo(i,:),mcmc,refPriors(i))
                       done(i)=.true.
                    endif    
                else    
                    Priors(i)=xrandom(PrInfo(i,:),mcmc,refPriors(i)) 
                    done(i)=.true.
                endif
             endif    
          enddo 
          q=COUNT(done)
          if (q==r) exit
       enddo
       !-------------------------------------------------- 
       !..and then we assign the priors: 
       !--------------------------------------------------
       !Checking samplinfo
       z=1
       do i=1,gr,1
          do j=1,3,1    
             if (samplinfo(j,i)==-6666666) then
                samplinfo(j,i)=Priors(z) 
                z=z+1
             endif    
          enddo          
       enddo    
       !--------------------------------------------------
       !Checking Ne and growth
       do i=1,dn,1
          if (Ne(i)==-6666666) then
             Ne(i)=Priors(z) 
             z=z+1
          endif
          if (growth(i)==-6666666) then
             growth(i)=Priors(z) 
             z=z+1
          endif
       enddo  
       !--------------------------------------------------
       !Checking events
       do i=1,en,1
          do j=1,8,1    
             if (events(i,j)==-6666666) then                 
                events(i,j)=Priors(z)                 
                if ((j==6).or.(j==7)) then  !If the Prior is in the proportion of lineages we have to adjust the proportion in the block that taes the remaining
                   do k=1,en,1              !So we check block by block...
                      do l=4,5,1            !Until we get the one that has the same entring block
                         if ((k==i).and.(events(k,l)==events(i,j-2))) then
                            events(k,l+2)=1-Priors(z)       
                         endif    
                      enddo
                   enddo    
                endif  
                z=z+1
             endif    
          enddo          
       enddo 

       !Substitution-related and other parameters
       if (gentime==-6666666) then 
          gentime=Priors(z) 
          z=z+1
       endif       
       if (sexratio==-6666666) then 
          sexratio=Priors(z) 
          z=z+1
       endif
              
       if (num==1) then
          refPriors=Priors
       endif 
       
       
       return
          
	end subroutine get_priors
        
	
	
	
	
    REAL(8) FUNCTION xrandom(info,mcmc,ref)
	
	
       USE irandom
       USE luxury
   
       IMPLICIT NONE
     
       !Variables
       REAL(8), DIMENSION(9), INTENT(IN):: info       !Info of the probability distribution of priors
       LOGICAL, INTENT(IN):: mcmc                     !If we are performing a MCMC without likelihoods then we sample from the jump density
       REAL(8), INTENT(IN):: ref                      !The previous prior value obtained (used as the mean for a normal jump distribution with MCMC)
     
       INTEGER:: i
       REAL(8):: x, y, sd
       REAL, DIMENSION(1):: RVEC
       
       !----------------------------------
       IF ((mcmc).AND.(info(8)==1.0)) THEN
           sd=info(9)
           SELECT CASE(int(info(1)))
           CASE(1) !Uniform
              x=random_normal()*sd+ref    !Value is sampled from a normal jump density
              x=MAX(info(6),MIN(info(7),x))   !Checking boundaries 
           CASE(2) !Log-uniform
              y=LOG(ref) 
              x=random_normal()*sd+y    !Value is sampled from a normal jump density              
              x=MAX(info(6),MIN(info(7),EXP(x))) !Checking boundaries    
           CASE(3) !Normal
              x=random_normal()*sd+ref    !Value is sampled from a normal jump density
              x=MAX(info(6),MIN(info(7),x))   !Checking boundaries 
           CASE(4) !Log-normal
              y=LOG(ref) 
              x=random_normal()*sd+y    !Value is sampled from a normal jump density              
              x=MAX(info(6),MIN(info(7),EXP(x))) !Checking boundaries  
           CASE(5) !Exponential
              x=random_normal()*sd+ref    !Value is sampled from a normal jump density
              x=MAX(info(6),MIN(info(7),x))   !Checking boundaries 
           CASE(6) !Gamma
              x=random_normal()*sd+ref    !Value is sampled from a normal jump density
              x=MAX(info(6),MIN(info(7),x))   !Checking boundaries
           CASE(7) !Geometric
              x=random_normal()*sd+ref    !Value is sampled from a normal jump density
              x=MAX(info(6),MIN(info(7),x))   !Checking boundaries
           CASE(8) !Binomial   
              x=random_normal()*sd+ref    !Value is sampled from a normal jump density
              x=MAX(info(6),MIN(info(7),x))   !Checking boundaries
           CASE(9) !Beta
              x=random_normal()*sd+ref    !Value is sampled from a normal jump density
              x=MAX(info(6),MIN(info(7),x))   !Checking boundaries
           CASE(10) !Poisson
              x=random_normal()*sd+ref    !Value is sampled from a normal jump density
              x=MAX(info(6),MIN(info(7),x))   !Checking boundaries
           ENDSELECT   
           !-----------------------------------------------------------
       ELSE
           SELECT CASE(int(info(1)))
           CASE(1) !Uniform
              CALL RANLUX (RVEC,1)
              x=RVEC(1)*(info(4)-info(3))+info(3) + info(5) !info(5)=offset
           CASE(2) !Log-uniform
              CALL RANLUX (RVEC,1)
              x=EXP(RVEC(1)*(info(4)-info(3))+info(3)) + info(5) !info(5)=offset   
           CASE(3) !Normal
              x=random_normal()*info(4)+info(3) +info(5) !info(3)=mean, info(4)=standar deviation, info(5)=offset
           CASE(4) !Log-normal
              x=EXP(random_normal()*info(4)+info(3)) + info(5) !info(3)=ln-mean, info(4)=ln-standard deviation, info(5)=offset  
           CASE(5) !Exponential
              x=random_exponential()*info(3) + info(5)   !we use lambda parametrized as the rate (rate=1/lambda), info(5)=offset
           CASE(6) !Gamma
              x=random_gamma(real(info(3)),.true.)
              x=x*info(4)+info(5) !info(3)=shape,info(4)=scale,info(5)=offset         
           CASE(7) !Geometric
              i=0
              do
                 i=i+1 
                 CALL RANLUX (RVEC,1)
                 if (RVEC(1)<info(3)) exit    
              enddo    
              x=i + info(5) !info(3)=p,info(5)=offset 
           CASE(8) !Binomial 
              x=random_binomial2(int(info(3)),real(info(4)),.true.) + info(5) !info(3)=n,info(4)=p,info(5)=offset 
              !x=random_binomial1(int(info(3)),real(info(4)),.true.) + info(5) !info(3)=n,info(4)=p,info(5)=offset 
           CASE(9) !Beta
              x=random_beta(real(info(3)),real(info(4)),.true.) + info(5) !info(3)=alpha,info(4)=beta,info(5)=offset
           CASE(10) !Poisson
              x=random_Poisson(real(info(3)),.true.) + info(5) !info(3)=lambda,info(5)=offset
           ENDSELECT
       ENDIF    
       x=x*info(2) !+ or -
       xrandom = MAX(info(6),MIN(info(7),x))   !Last checking to boundaries. Here always in natural scale
       !-----------------------------     
           
       
       RETURN          
                
    END FUNCTION xrandom  
    
	
	
	
	
    subroutine get_matches(dn,Ne,growth,en,events)
   
       implicit none
   
       !Variables
       integer, intent(in)::dn
       real(8), dimension(dn), intent(in)::Ne
       real(8), dimension(dn), intent(inout):: growth
       integer, intent(in):: en
       real(8), dimension(en,8):: events
       
       integer:: i, j
       integer, dimension(1):: pos
       integer:: a, b
       real(8):: Ne1, Ne2
       real(8):: g1, g2
       real(8):: t, t1, t2
       real(8):: factor
       !Start
       !First we check in the growth column
       do i=1,dn,1
          if (growth(i)==-7777777) then
             a=0
             b=0
             do j=1,en !This instruction is only to get the id's of the blocks that have our current block(i-th) as an entring block
                if ((events(j,4)==i).or.(events(j,5)==i)) then
                   if (a==0) then    
                      a=j
                   else
                      b=j
                   endif
                endif   
             enddo 
             if (events(a,4)/=events(a,5)) then !In case the ancient a-th block takesl ineages also from another block/population, we adjust its proportion of Ne
                if (events(a,4)==i) then !Before getting the Ne, we need what porportion of that Ne is contribution of our entring block 
                   factor=events(a,6) !Which is in column 6 if our entring block id is in cloumn 4
                else
                   factor=events(a,7) !..or it is in column 7 if our entring block id was in col 5
                endif
                factor=factor/(events(a,6)+events(a,7))
             else
                factor=1
             endif
             Ne1=events(a,2)*factor !The Ne is in events(:,2)
             if (b/=0) then
                if (events(b,4)/=events(b,5)) then !In case the ancient a-th block takesl ineages also from another block/population, we adjust its proportion of Ne
                   if (events(b,4)==i) then !Before getting the Ne, we need what porportion of that Ne is contribution of our entring block 
                      factor=events(b,6) !Which is in column 6 if our entring block id is in cloumn 4
                   else
                      factor=events(b,7) !..or it is in column 7 if our entring block id was in col 5
                   endif
                   factor=factor/(events(b,6)+events(b,7))
                else
                   factor=1
                endif
                Ne2=events(b,2)*factor !The Ne is in events(:,2)
             else
                Ne2=0 
             endif
             Ne2=Ne1+Ne2
             Ne1=Ne(i)
             t=events(a,1) !-events(i,1)
             !Finnally we calculate the bloody number:
             growth(i)=(log(Ne2)-log(Ne1))/t  
          endif
       enddo    
       !-----------------------------------------------------------------------------
       do i=1,en,1  !We check all the events looking for a Match code (-7777);
          if (events(i,2)==-7777777) then  !If the Match is in column 2, it is a Ne so: 
             a=events(i,4)      !a & b gets the id of the blocks whose Ne have to match (the entring blocks of our block)
             b=events(i,5)
             if (a<=dn) then    !if the entring block is a population
                g1=growth(a)    !we get growth, Ne and age of the population
                Ne1=Ne(a)
                t1=events(i,1)
             else               !if the entring block is an internal block, the data is in events   
                pos=MAXLOC(events(:,8), Mask= events(:,8)==a) !The entring block is the 
                g1=events(pos(1),3)
                Ne1=events(pos(1),2)
                t1=events(i,1)-events(pos(1),1)
             endif
             if (b/=a) then    !The same made with 'a' is done for 'b' except if it is = a
                 if (b<=dn) then
                    g2=growth(b)
                    Ne2=Ne(b)
                    t2=events(i,1)
                 else
                    pos=MAXLOC(events(:,8), Mask= events(:,8)==b)
                    g2=events(pos(1),3)
                    Ne2=events(pos(1),2)
                    t2=events(i,1)-events(pos(1),1)
                 endif
             else              !if b=a we destroy the secont part of the sum (cuz is repeated)
                g2=1
                Ne2=1
                t2=1
                factor=0
             endif    
             events(i,2) = events(i,6)*(Ne1*(2.718281828459**(t1*g1))) 
             events(i,2) = events(i,2) + events(i,7)*(Ne2*(2.718281828459**(t2*g2))) 
          endif   
          if (events(i,3)==-7777777) then !if Match is in column 3 it is a growth rate
             a=0
             b=0
             do j=i+1,en !This instruction iso nly to get the id's of the blocks that have our current block(i-th) as an entring block
                if ((events(j,4)==events(i,8)).or.(events(j,5)==events(i,8))) then
                   if (a==0) then    
                      a=j
                   else
                      b=j
                   endif
                endif   
             enddo 
             if (events(a,4)==i) then !Before getting the Ne, we need what porportion of that Ne is contribution of our entring block 
                factor=events(a,6) !Which is in column 6 if our entring block id is in cloumn 4
             else
                factor=events(a,7) !..or it is in column 7 if our entring block id was in col 5
             endif
             Ne1=events(a,2)*factor !The Ne is in events(:,2)
             if (b/=0) then
                if (events(b,4)==i) then !Before getting the Ne, we need what porportion of that Ne is contribution of ou entring block 
                   factor=events(b,6) !Which is in column 6 if our entring block id is in cloumn 4
                else
                   factor=events(b,7) !..or it is in column 7 if our entring block id was in col 5
                endif
                Ne2=events(b,2)*factor !The Ne is in events(:,2)
             else
                Ne2=0 
             endif
             Ne2=Ne1+Ne2
             Ne1=events(i,2)
             t=events(a,1)-events(i,1)
             !Finnally we calculate the bloody number:
             events(i,3)=(log(Ne2)-log(Ne1))/t
          endif    
       enddo
       
       
   
       return
   
    end subroutine get_matches
    
    
    
    
    
    subroutine PriorProbability(r,PrInfo,Priors,priorL)
      
      use irandom
      
      implicit none      
      
      !Variables
      integer, intent(in):: r                                   !Number of priors 
      real(8), dimension(r,9), intent(in):: PrInfo              !Array carrying the information about priors
      real(8), dimension(r), intent(in):: Priors                !Array carrying the sampled values of the priors
      real(16), intent(out):: priorL 
      
      integer:: i                !Counter
      real(16):: x, y, z, a, b   !Multipurpose
      !-----------------------------------------------
      !Start
      priorL=0.0
      do i=1,r,1
         if (PrInfo(i,8)==1) then
            select case(int(PrInfo(i,1)))
            case(1)   !Uniform
               a=Prinfo(i,3)
               b=Prinfo(i,4)
               z=-LOG(b-a)
            case(2)   !Log-uniform
               a=Prinfo(i,3)
               b=Prinfo(i,4)
               z=-LOG(b-a)    
            case(3)   !Normal
               x=Priors(i)*PrInfo(i,2)-PrInfo(i,5)  !We bring back the sampled value (without offset or inversion)
               a=Prinfo(i,3)   !median
               b=Prinfo(i,4)   !s.d.
               y=-LOG(b*2.506628274631000502415765284811)
               z=-0.5*((x-a)/b)**2    
               z=z+y               
            case(4)   !Log-normal
               x=Priors(i)*PrInfo(i,2)-PrInfo(i,5)  !We bring back the sampled value (without offset or inversion)
               x=LOG(x)
               a=Prinfo(i,3)   !median
               b=Prinfo(i,4)   !s.d.
               y=-LOG(b*2.506628274631000502415765284811)
               z=-0.5*((x-a)/b)**2    
               z=z+y
            case(5)   !Exponential
               x=Priors(i)*PrInfo(i,2)-PrInfo(i,5)
               a=Prinfo(i,3)**-1     !recall we use rate (rate=1/lambda)
               z=log(a)-a*x 
            case(6)   !Gamma
               x=Priors(i)*PrInfo(i,2)-PrInfo(i,5) 
               a=Prinfo(i,3)   !shape
               b=Prinfo(i,4)   !scale
               y=lngamma(a)    !gamma function
               z=-y-a*log(b)+(a-1)*log(x)-(x/b)
            case(7)   !Geometric
               x=Priors(i)*PrInfo(i,2)-PrInfo(i,5) 
               a=Prinfo(i,3)   !p
               z=log(a)+(x-1)*log(1-a)
            case(8)    !Binomial   
               x=Priors(i)*PrInfo(i,2)-PrInfo(i,5) 
               a=Prinfo(i,3)   !n
               b=Prinfo(i,4)   !p
               y=lngamma(a+1)-lngamma(x+1)-lngamma(a-x+1)  !This is ln(Comb[n,k])
               z=y+x*log(b)+(a-x)*log(1-b)
            case(9)    !Beta
               x=Priors(i)*PrInfo(i,2)-PrInfo(i,5) 
               a=Prinfo(i,3)   !alpha (shape)
               b=Prinfo(i,4)   !beta (shape)
               y=-lngamma(a)-lngamma(b)+lngamma(a+b)
               z=(a-1)*log(x)+(b-1)*log(1-x)+y
            case(10)   !Poisson
               x=Priors(i)*PrInfo(i,2)-PrInfo(i,5) 
               a=Prinfo(i,3)   !lambda
               z=x*log(a)-a-lngamma(x+1) 
           endselect
           priorL=priorL+z
         endif  
         
      enddo
   
   
      return
      
    end subroutine PriorProbability
      
   
   
    
    
   subroutine CoalLikelihood(n,sa,en,samplinfo,coal_t,NeChanges,iniNe,inigrowth,gentime,LnL)  !Ne,growth,e,events,gen,adjust,eages,coal_t,LnL)
   
      implicit none
      
      !Variables
      integer, intent(in):: n                             !Sample size
      integer, intent(in):: sa                            !Number of sample age groups
      integer, intent(in):: en                            !Number of events
      real(8), dimension(3,sa), intent(in):: samplinfo
      real(8), dimension(n-1), intent(in):: coal_t
      real(8), dimension(en,8), intent(in):: NeChanges    ! NeChanges=ievents in the main part 
      real(8), intent(in):: iniNe, inigrowth, gentime
      real(16), intent(out):: LnL      
      
      integer:: i, j, d, lin                              !Multipurpose
      real(8):: x, y, Ne, N0, g, t, growth, tgen, lambda
      real(8), allocatable, dimension(:):: narr, times    !Provisional arrays to sort times
      real(8), allocatable, dimension(:,:):: narray, marray       !Master array:
                                                          !narray    marray 
                                                          !1st col   1st col:= time start block;
                                                          !          2st col:= time end block;
                                                          !2nd col   3rd col:= number of lineages introduced by the event;
                                                          !          4th col:= initial Ne (youngest);
                                                          !          5th col:= final Ne (oldest); 
                                                          !          6th col:= growth rate
                                                          !3rd col   7th col:= type (0:=ages; 1:=coal; 2:=demi) 
      logical:: done
      !------------------------     
      !(1) First we prepare master array (marr) with all relevant events (coalescent, new sample, Ne change)
      d=(n-1)+sa+en
      allocate(narr(d),times(d),narray(d,3),marray(d-1,7))  
      !Listing the new samples introduced by being ancient
      do i=1,sa,1                   
         narray(i,1)=samplinfo(3,i)  !samplinfo ages
         narray(i,2)=samplinfo(1,i)  !samplinfo subsample size
         narray(i,3)=0  
      enddo   
      !Listing the coalescent events from coal_t
      do i=1,n-1,1
         narray(sa+i,1)=coal_t(i)    !time to the ith coalescent
         narray(sa+i,2)=-1           !each coalescent removes a lineage
         narray(sa+i,3)=1  
      enddo
      !Listing the events of Ne change
      do i=1,en,1
         narray(sa+n-1+i,1)=NeChanges(i,1)  
         narray(sa+n-1+i,2)=0
         narray(sa+n-1+i,3)=2
      enddo
      narr=(/(i,i=1,d)/)
      times=narray(:,1)
      CALL SSORT(times,narr,d,2)
      !Now we can make master array:
      !This is only to find the
      i=1  !i counts the rows in narray and marray
      j=0  !j counts the rows that are coalescents
      lin=0
      do
         marray(i,1)=narray(narr(i),1)
         marray(i,2)=narray(narr(i+1),1)
         lin=lin+narray(narr(i),2) 
         marray(i,3)=lin
         marray(i,7)=narray(narr(i),3)
         if (marray(i,7)==1) then
            j=j+1
         endif   
         if ((j==n-2).and.(narray(narr(i+1),3)==1)) then !i.e. the next event is the last coalescent
            done=.true.
         endif   
         if (done) exit
         i=i+1
      enddo  
      !Now we have to set the Ne values and growth
      Ne=iniNe
      growth=inigrowth
      t=NeChanges(1,1)      
      j=0
      do i=1,d-1,1
         if (marray(i,1)>=t) then  !Update Ne, growth, threshold t             
            j=j+1 
            Ne=NeChanges(j,2)            
            growth=NeChanges(j,3)
            if (j==en) then
               t=MAXVAL(times)+1000000 
            else    
               t=NeChanges(j+1,1)
            endif            
         endif
         marray(i,4)=Ne
         if (growth==0.0) then
            marray(i,5)=Ne
         else  
            tgen=(marray(i,2)-marray(i,1))/gentime
            marray(i,5)=Ne*exp(tgen*growth) 
         endif 
         marray(i,6)=growth
      enddo         
      !+++++++++++++++++++++++++++++++++
      !Finally the likelihood
      LnL=0.0
      do i=1,d-1,1
         if (marray(i,6)==0.0) then !growth rate
            Ne=marray(i,4)    
         else
            N0=marray(i,4)   !this is Ne(0) in N(t)=N(0)e^gt 
            t=(marray(i,2)-marray(i,1))/gentime !this is t (time elapsed) 
            g=marray(i,6)
            y=(N0/g)*(exp(g*t)-1)    !y gets the area under the curve 
            Ne=y/t    !area y=Ne*t
         endif    
         y=marray(i,3)
         lambda=y*(y-1)/2
         x=(marray(i,2)-marray(i,1))/gentime
         x=lambda*x/Ne
         x=log(lambda)-log(Ne)-x 
         LnL=LnL+x
      enddo 
        
      
      return
   
   end subroutine CoalLikelihood
   
   
   
   
 
    SUBROUTINE SSORT (X, Y, N, KFLAG) 
    !***BEGIN PROLOGUE SSORT 
    !***PURPOSE Sort an array and optionally make the same interchanges in 
    ! an auxiliary array. The array may be sorted in increasing 
    ! or decreasing order. A slightly modified QUICKSORT 
    ! algorithm is used. 
    !***LIBRARY SLATEC 
    !***CATEGORY N6A2B 
    !***TYPE SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I) 
    !***KEYWORDS SINGLETON QUICKSORT, SORT, SORTING 
    !***AUTHOR Jones, R. E., (SNLA) C Wisniewski, J. A., (SNLA) 
    !***DESCRIPTION 
    ! 
    ! SSORT sorts array X and optionally makes the same interchanges in 
    ! array Y. The array X may be sorted in increasing order or 
    ! decreasing order. A slightly modified quicksort algorithm is used. 
    ! 
    ! Description of Parameters 
    ! X - array of values to be sorted (usually abscissas) 
    ! Y - array to be (optionally) carried along 
    ! N - number of values in array X to be sorted 
    ! KFLAG - control parameter 
    ! = 2 means sort X in increasing order and carry Y along. 
    ! = 1 means sort X in increasing order (ignoring Y) 
    ! = -1 means sort X in decreasing order (ignoring Y) 
    ! = -2 means sort X in decreasing order and carry Y along. 
    ! 
    !***REFERENCES R. C. Singleton, Algorithm 347, An efficient algorithm 
    ! for sorting with minimal storage, Communications of 
    ! the ACM, 12, 3 (1969), pp. 185-187. 
    !***REVISION HISTORY (YYMMDD) 
    ! 761101 DATE WRITTEN 
    ! 761118 Modified to use the Singleton quicksort algorithm. (JAW) 
    ! 890531 Changed all specific intrinsics to generic. (WRB) 
    ! 890831 Modified array declarations. (WRB) 
    ! 891009 Removed unreferenced statement labels. (WRB) 
    ! 891024 Changed category. (WRB) 
    ! 891024 REVISION DATE from Version 3.2 
    ! 891214 Prologue converted to Version 4.0 format. (BAB) 
    ! 900315 CALLs to XERROR changed to CALLs to XERMSG. (THJ) 
    ! 901012 Declared all variables; changed X,Y to SX,SY. (M. McClain) 
    ! 920501 Reformatted the REFERENCES section. (DWL, WRB) 
    ! 920519 Clarified error messages. (DWL)
    ! 920801 Declarations section rebuilt and code restructured to use 
    ! IF-THEN-ELSE-ENDIF. (RWC, WRB) 
    !***END PROLOGUE SSORT 
    ! .. Scalar Arguments .. 
    INTEGER KFLAG, N 
    ! .. Array Arguments .. 
    REAL(8) X(*), Y(*) 
    ! .. Local Scalars .. 
    REAL(8) R, T, TT, TTY, TY 
    INTEGER I, IJ, J, K, KK, L, M, NN 
    ! .. Local Arrays .. 
    INTEGER IL(21), IU(21) 
    ! .. External Subroutines .. 
    ! None 
    ! .. Intrinsic Functions .. 
    INTRINSIC ABS, INT 
    !***FIRST EXECUTABLE STATEMENT SSORT 
    NN = N 
    IF (NN .LT. 1) THEN
       write(*,*) 'The number of values to be sorted is not positive.' 
       RETURN 
    ENDIF 
    ! 
    KK = ABS(KFLAG) 
    IF (KK.NE.1 .AND. KK.NE.2) THEN
       write(*,*) 'The sort control parameter, K, is not 2, 1, -1, or -2.'
       RETURN 
    ENDIF 
    ! 
    ! Alter array X to get decreasing order if needed 
    ! 
    IF (KFLAG .LE. -1) THEN
       DO 10 I=1,NN
          X(I) = -X(I) 
10     CONTINUE 
    ENDIF 
    ! 
    IF (KK .EQ. 2) GO TO 100  !C Sort X only 
    ! 
    M = 1 
    I = 1 
    J = NN 
    R = 0.375E0 
    ! 
20  IF (I .EQ. J) GO TO 60 
    IF (R .LE. 0.5898437E0) THEN
       R = R+3.90625E-2 
    ELSE 
       R = R-0.21875E0 
    ENDIF 
    ! 
30  K = I   ! C Select a central element of the array and save it in location T 
    ! 
    IJ = I + INT((J-I)*R) 
    T = X(IJ)  ! C If first element of array is greater than T, interchange with T  
    !
    IF (X(I) .GT. T) THEN
       X(IJ) = X(I)
       X(I) = T 
       T = X(IJ) 
    ENDIF 
    L = J    ! C If last element of array is less than than T, interchange with T 
    ! 
    IF (X(J) .LT. T) THEN
       X(IJ) = X(J) 
       X(J) = T 
       T = X(IJ)  ! C If first element of array is greater than T, interchange with T 
       ! 
       IF (X(I) .GT. T) THEN 
          X(IJ) = X(I) 
          X(I) = T 
          T = X(IJ) 
       ENDIF 
    ENDIF  ! C Find an element in the second half of the array which is smaller 
          ! than T 
          ! 
40  L = L-1 
    IF (X(L) .GT. T) GO TO 40    ! C Find an element in the first half of the array which is greater 
                                 ! than T 
                                 ! 
50  K = K+1
    IF (X(K) .LT. T) GO TO 50    ! C Interchange these elements 
                                   ! 
    IF (K .LE. L) THEN
       TT = X(L)
       X(L) = X(K)
       X(K) = TT
       GO TO 40 
    ENDIF     ! C Save upper and lower subscripts of the array yet to be sorted 
                ! 
    IF (L-I .GT. J-K) THEN
       IL(M) = I 
       IU(M) = L 
       I = K 
       M = M+1 
    ELSE 
       IL(M) = K 
       IU(M) = J 
       J = L 
       M = M+1 
    ENDIF 
    GO TO 70        ! C Begin again on another portion of the unsorted array 
                      ! 
60  M = M-1 
    IF (M .EQ. 0) GO TO 190 
    I = IL(M) 
    J = IU(M) 
    ! 
70  IF (J-I .GE. 1) GO TO 30 
    IF (I .EQ. 1) GO TO 20 
    I = I-1
    ! 
80  I = I+1
    IF (I .EQ. J) GO TO 60
    T = X(I+1)
    IF (X(I) .LE. T) GO TO 80 
    K = I 
    ! 
90  X(K+1) = X(K) 
    K = K-1 
    IF (T .LT. X(K)) GO TO 90 
    X(K+1) = T 
    GO TO 80    ! C Sort X and carry Y along 
                             ! 
100 M = 1
    I = 1 
    J = NN
    R = 0.375E0 
    ! 
110 IF (I .EQ. J) GO TO 150 
    IF (R .LE. 0.5898437E0) THEN
       R = R+3.90625E-2 
    ELSE 
       R = R-0.21875E0 
    ENDIF
    ! 
120 K = I    ! C Select a central element of the array and save it in location T 
               ! 
    IJ = I + INT((J-I)*R) 
    T = X(IJ) 
    TY = Y(IJ)   ! C If first element of array is greater than T, interchange with T 
                   ! 
    IF (X(I) .GT. T) THEN
       X(IJ) = X(I) 
       X(I) = T
       T = X(IJ)
       Y(IJ) = Y(I)
       Y(I) = TY
       TY = Y(IJ)
    ENDIF 
    L = J     ! C If last element of array is less than T, interchange with T
                      ! 
    IF (X(J) .LT. T) THEN
       X(IJ) = X(J)
       X(J) = T 
       T = X(IJ)
       Y(IJ) = Y(J)
       Y(J) = TY
       TY = Y(IJ)    ! C If first element of array is greater than T, interchange with T 
                       ! 
       IF (X(I) .GT. T) THEN
          X(IJ) = X(I) 
          X(I) = T 
          T = X(IJ) 
          Y(IJ) = Y(I) 
          Y(I) = TY 
          TY = Y(IJ)
       ENDIF 
    ENDIF      ! C Find an element in the second half of the array which is smaller 
                 ! than T 
                 ! 
130 L = L-1 
    IF (X(L) .GT. T) GO TO 130    ! C Find an element in the first half of the array which is greater 
                                    ! than T 
                                    ! 
140 K = K+1 
    IF (X(K) .LT. T) GO TO 140    ! C Interchange these elements 
                                    ! 
    IF (K .LE. L) THEN
       TT = X(L) 
       X(L) = X(K)
       X(K) = TT
       TTY = Y(L) 
       Y(L) = Y(K) 
       Y(K) = TTY 
       GO TO 130
    ENDIF               ! C Save upper and lower subscripts of the array yet to be sorted 
                          ! 
    IF (L-I .GT. J-K) THEN
       IL(M) = I 
       IU(M) = L 
       I = K 
       M = M+1
    ELSE
       IL(M) = K 
       IU(M) = J 
       J = L 
       M = M+1 
    ENDIF 
    GO TO 160      ! C Begin again on another portion of the unsorted array
                     ! 
150 M = M-1
    IF (M .EQ. 0) GO TO 190
    I = IL(M) 
    J = IU(M) 
    ! 
160 IF (J-I .GE. 1) GO TO 120
    IF (I .EQ. 1) GO TO 110 
    I = I-1 
    ! 
170 I = I+1
    IF (I .EQ. J) GO TO 150 
    T = X(I+1) 
    TY = Y(I+1)
    IF (X(I) .LE. T) GO TO 170
    K = I 
    ! 
180 X(K+1) = X(K) 
    Y(K+1) = Y(K) 
    K = K-1
    IF (T .LT. X(K)) GO TO 180 
    X(K+1) = T 
    Y(K+1) = TY 
    GO TO 170        ! C Clean up 
                       ! 
190 IF (KFLAG .LE. -1) THEN
       DO 200 I=1,NN
          X(I) = -X(I)
200    CONTINUE
    ENDIF
   
    RETURN 

END SUBROUTINE SSORT
  
   
   
   
  
   SUBROUTINE ErrorMessage(code)
    
        !USE ifport
        !USE dflib
        !USE ifqwin

        IMPLICIT NONE
  
        INTEGER, INTENT(IN):: code
        integer:: i, ret, ret2       
        !------------------------------
        SELECT CASE(code)
              CASE(1002)
                 write(*,*) 'Error!: ', 'The simulations program could not read the options provided by the interface. Please, check the names and content of iSuS file.' 
                 !ret2 = CLICKMENUQQ(loc(WINEXIT))   
                 STOP
              CASE(1003)
                 write(*,*) 'Error!: ', 'The simulations program could not open the info file. Please check that a file named iSuS appears in the same directory.'   
                 !ret = CLICKMENUQQ(loc(WINEXIT))   
                 STOP
              CASE(1005)
                 write(*,*) 'Error!: ', 'The input file could not be open or was not in the folder. Please check that input file is a ".dat" file and is located in the folder.' 
                 !ret = CLICKMENUQQ(loc(WINEXIT)) 
                 STOP
              CASE(1006)
                 write(*,*) 'Error!: ', 'The probability distribution of the prior was not recognized. Please, chech the input file and consult the handbook if necessary.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(10061)
                 write(*,*) 'Error!: ', 'The status of the prior (8th column) should be either NOVAR (noise variable) or PAR (parameter). Please, add one of those two commands.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP   
              CASE(1007)
                 write(*,*) 'Error!: ', 'The sample sizes provided at the head of the input file and calculated from the sampling section do not coincide.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
			  CASE(1008)
                 write(*,*) 'Error!: ', 'The generation time cannot be negative. Please, check the input file.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT)) 
                 STOP
              CASE(1009)
                 write(*,*) 'Error!: ', 'The mutation rate cannot be negative. Please, check the input file.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT)) 
                 STOP
              CASE(1010)
                 write(*,*) 'Error!: ', 'Please, check the number of nucleotides (lenght) of the DNA fragment.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))
                 STOP
              CASE(1011)
                 write(*,*) 'Error!: ', 'The proportion of invariant sites should be between zero (inclusive) to one.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT)) 
                 STOP
              CASE(1012)
                 write(*,*) 'Error!: ', 'Shape parameter of (gamma parameter) should be positive.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT)) 
                 STOP
              CASE(1013)
                 write(*,*) 'Error!: ', 'The number of migration matrix should be a non-negative integer.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1014)
                 write(*,*) 'Error!: ', 'The migration rates should be given as probabilities (i.e. numbers between 0.0 and 1.0).' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1015)
                 write(*,*) 'Error!: ', 'The name of the prior density should be followed by a "-" or "+" symbol indicating if the sampled parameter should be mirrored or not.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1016)
                 write(*,*) 'Error!: ', 'The number of populations provided at the head of the input file and gotten from the sampling section do not coincide.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1017)
                 write(*,*) 'Error!: ', 'The labels of the populations should be consecutive integers starting in 1 (1, 2, 3, ...).' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1018)
                 write(*,*) 'Error!: ', 'The effective population sizes should be positive and larger than 1.0.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1019)
                 write(*,*) 'Error!: ', 'The times to the events should be positive (they are considered the age of the event).' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1020)
                 write(*,*) 'Error!: ', 'The effective population sizes of the blocks should be larger than 1.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT)) 
                 STOP
              CASE(1021)
                 write(*,*) 'Error!: ', 'The labels of the blocks/populations that are a source of another block should correspond to a numeration continuing the populations numbers.' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1022)
                 write(*,*) 'Error!: ', 'The columns 6th and 7th of events section should contain the contribution of lineages of the source block/population (it is a proportion).' 
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1023)
                 write(*,*) 'Error!: ', 'The labels of the events/blocks should correspond to a numeration continuing the populations numeration.'
				 !ret = CLICKMENUQQ(loc(WINEXIT))  
                 STOP
              CASE(1024)
                 write(*,*) 'Error!: ', 'The FASTA file was not found. Please, check the name and location of the alignment file.'
				 !ret = CLICKMENUQQ(loc(WINEXIT))
                 STOP
              CASE(1025)
                 write(*,*) 'Error!: ', 'The blocks cannot be source of lineages for themlseves (columns 4th and 5th of blocks/events should not coincide with column 8th)'  
				 !ret = CLICKMENUQQ(loc(WINEXIT)) 
                 STOP
              CASE(1026)
                 write(*,*) 'Error!: ', 'When using the "Random_Ne_Hist" command, the number of random changes should be an integer positive or a prior'  
				 !ret = CLICKMENUQQ(loc(WINEXIT)) 
                 STOP
              CASE(1027)
                 write(*,*) 'Error!: ', 'When using the "Random_Ne_Hist" command, the interval of time should be set by two positive integers in increasing order'  
				 !ret = CLICKMENUQQ(loc(WINEXIT)) 
                 STOP   
              CASE(9999)
               !CALL BEEPQQ(2500,150)
               !CALL BEEPQQ(1000,200)
               !CALL BEEPQQ(2500,150)
               !CALL BEEPQQ(900,200)
               !CALL BEEPQQ(500,200)
               !CALL BEEPQQ(2200,750)
               write(*,*) 'THE ANALYSIS HAS FINISHED SUCCESFULLY! :)'              
               !CALL SLEEPQQ(9000000)  
        END SELECT          
        
        RETURN

   END SUBROUTINE ErrorMessage

   
   
   
    
  
    
       
    
end module accessories
