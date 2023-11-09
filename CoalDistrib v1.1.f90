!  CoalDistribv11.f90 
!
!  FUNCTIONS:
!  CoalDistribv11 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: CoalDistribv11
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program CoalDistrib11
     
     use accessories
     use luxury
      
     implicit none

     !Variables
     character(30):: nameA, nameB, iname                            !Name of scenario file   
     integer:: n, dn, en, gr, r                                     !Overall sample size, nr of demes, nr of events, nr of sample groups, nr of priors 
     integer:: chainLength, int1                                    !Number of iterations of the MCMC  
     real(8), allocatable, dimension(:,:):: samplinfo, isamplinfo   !Sampling info: [1st col.] sample group sizes, [2nd col] deme number, [3rd col] group ages 
     real(8), allocatable, dimension(:):: Ne, iNe                   !Vector of population sizes 
     real(8), allocatable, dimension(:):: growth, igrowth           !Vector of growth rates 
     real(8), allocatable, dimension(:):: bid                       !Identifier of population blocks for Ne/growth
     real(8), allocatable, dimension(:,:):: events, ievents         !Info about the internal blocks
     real :: gen, igen                                              !Generation time 
     real:: sexratio, isexratio, adjust                             !Sexratio, the order male-female depends on the marker focus, e.g. 0.2 mitochondrial means 4/1 (female/male)
     integer:: marker                                               !Marker type: 1)Autosomal haploid; 2)Autosomal diploid; 3)Haplo-diploid/X-linked; 4)Mitochondrial/Y-linked
     character(30), allocatable, dimension(:):: labels              !Labels of organisms (sequences)     
     real(8), allocatable, dimension(:,:):: PrInfo, iPrInfo         !Array carrying the information about priors
     real(8), allocatable, dimension(:):: Priors, refPriors         !Array with the random values simulated from the priors          
     real(8), allocatable, dimension(:):: empCoal                   !Array carrying the coalescent times of the empirical data  
     real(8), allocatable, dimension(:,:):: table                   !The table with parameters and likelihoods
     integer:: seed1                                                !seed of the random number generator  
     character(18), allocatable, dimension(:):: headings            !Headings of the results file     
     real(8), allocatable, dimension(:):: ages                      !Ages of each taxa (individual) in the sample (ordered)
     real(8), allocatable, dimension(:):: coal_t                    !Array storing the coalescent times                                 
     integer:: i, j, k, m  
     integer:: ierror
     real(16):: newL, oldL, priorL, coalL, x                        !For storing likelihoods/probabilities  
     logical:: esta, accepted
     character(8):: text    
     real, dimension(1):: rvec     
     !---------------------------------------------------------------------------------------------------
     !Body of CoalDistribv10-----------------------------------------------------------------------------
     write(*,*) 'Enter input file name'
     read(*,*) nameA
     nameB=TRIM(nameA)//'.dat'
     call preread_input(nameB,gr,dn,en,n,r)
     allocate(samplinfo(3,gr),isamplinfo(3,gr),Ne(dn),iNe(dn),growth(dn),igrowth(dn),bid(dn),events(en,8),ievents(en,8))
     allocate(labels(n),PrInfo(r,9),Priors(r),refPriors(r),headings(r+4),empCoal(n-1))
     call read_info(nameB,dn,gr,en,n,Ne,growth,bid,samplinfo,events,gen,sexratio,marker,labels,r,PrInfo,headings,empCoal,chainLength,int1,seed1)
     !MCMC***********************************
     allocate(table(chainLength+1,r+4))
     call RLUXGO(4,seed1,0,0)
     do i=1,chainLength,1
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         isamplinfo=samplinfo
         iNe=Ne
         igrowth=growth            
         ievents=events
         igen=gen
         isexratio=sexratio
         iPrInfo=PrInfo
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
         !Calling the routine that assigns priors 
         if (r>0) then
            call Get_Priors(i,r,gr,isamplinfo,dn,iNe,igrowth,en,ievents,igen,isexratio,iPrInfo,refPriors,Priors)
         endif 
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Resolving 'Match' commands 
         if ( (ANY(ievents(:,2)==-7777777)).or.(ANY(ievents(:,3)==-7777777)).or.(ANY(igrowth==-7777777)) ) then
            call Get_Matches(dn,iNe,igrowth,en,ievents)
         endif  
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
         !Adjustments for marker type 
         select case (marker) 
         case(1)    
            adjust = isexratio*(1-isexratio)*4     !Autosomal haploid
         case(2)
            adjust = isexratio*(1-isexratio)*8     !Autosomal diploid
         case(3)
            adjust = isexratio*(1-isexratio)*9/(1+isexratio) !Haplo-diploid/X-linked
         case(4)
            adjust = isexratio                     !Mitochondrial/Y-linked
         end select          
         ievents(:,2)=ievents(:,2)*adjust
         iNe=iNe*adjust 
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++         
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! M   C   M   C      
         call PriorProbability(r,PrInfo,Priors,priorL) 
         call CoalLikelihood(n,gr,en,isamplinfo,empCoal,ievents,dble(iNe(1)),igrowth(1),dble(igen),coalL)
         newL=coalL+priorL
         !MCMC jump
         if (i>1) then  
            call RANLUX (rvec,1) 
            x=newL-oldL            ! x gets the ln(f(x')/f(xt))  ; x'=proposal, xt=t-th state of the chain 
            if (x<-32) then        ! e^nL > (e^oL)^32
               accepted=.false.  
            elseif (x>32) then     ! e^nL > (e^oL)^32
               accepted=.true. 
            elseif (exp(x)>=1.0) then        !f(x')/f(xt) >= 1.0
               accepted=.true.
            elseif (exp(x)>=rvec(1)) then    !f(x')/f(xt) >= random[0,1]  (is accepted with probability=f(x')/f(x))
               accepted=.true.
            else   
               accepted=.false. 
            endif               
         endif
         if ((accepted).or.(i==1)) then
            oldL=newL    
            refPriors=Priors
         endif 
         if (i==1) then
            text='initial'
         elseif (accepted) then   
            text='accepted'
         else
            text='rejected'
         endif 
         !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++          
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if ((MOD(i,int1)==0).or.(i==1)) then
            write(*,1000) (real(i)/chainLength)*100,'% ',NewL 
            1000 format(1X,T1,F10.4,T12,A1,T15,F15.7)
            !Now we append the results
            if (i==1) then !This is because we only want to enter here the first time
               iname=TRIM(nameA)//'.txt'
               INQUIRE(FILE=iname, EXIST=esta)   
               if (esta) then
	 	         OPEN(UNIT=1, FILE=iname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ierror) 
 		       else
		         OPEN(UNIT=1, FILE=iname, STATUS='NEW', ACTION='WRITE', IOSTAT=ierror)  
		       endif	   
            else
               OPEN(UNIT=1, FILE=iname, STATUS='OLD', ACTION='WRITE', ACCESS='APPEND', IOSTAT=ierror)
            endif 
            openif0: if (ierror==0) then
               if (i==1) then  !This is because we only want to enter here the first time
                  write(1,1010,iostat=ierror) '# ','status ', 'L ', 'coalescent_L ', 'prior_L ', (headings(j),j=1,r)
               endif
               write(1,1011,iostat=ierror)  i, text, NewL, coalL, priorL, (Priors(j), j=1,r)
               1010 format(1x,1000(A18))  
               1011 format(1x,T1,I15,T18,A8,1000(F14.3))
            else
	           write(*,*)'Failed to write the table'
            endif openif0
            CLOSE(UNIT=1)
         endif   
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
     enddo
     !MCMC***********************************
     deallocate(samplinfo,isamplinfo,Ne,iNe,growth,igrowth,bid,events,ievents)
     deallocate(labels,PrInfo,Priors,refPriors,headings,table,empCoal)
     
     
     
     call ErrorMessage(9999)
     STOP
     
    end program CoalDistrib11

