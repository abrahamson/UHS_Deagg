      program Haz45Post

C This program will compute the equal hazard spectra and as well
C the equal deaggregation results for a given suite of PSHA runs. 
C This version is compatible with the output files from Haz45
c which only outputs the mean hazard code combined over all attenuation models. 

C     Version 45.2 (7/2016)

      implicit none
      include 'max_dims.H'

      integer nMagbins, nDistBins, nEpsBins, nXCostBins, bflag, 
     1        count1, count2, nProb, nAmp(MAXPROB), testnum,
     2        nSpecPer, ncurve(MAXPROB), iProb, iTest, nTest, 
     3        iAtten, nAmpper(MAXPROB), idum, i, nflt, j, k, l,
     4        imag, idis, ieps, iXco, iMagBin, iDistBin, iInten
      integer nskip, iMagbins, iDistbins, m, ii

      real amp(MAXPROB,MAXAMP), nEv(MAXPROB,MAXAMP), x, x1, x2,
     1     mBar(MAXPROB,MAXAMP), dBar(MAXPROB,MAXAMP), y1, y2, 
     2     eBar(MAXPROB,MAXAMP), risk(MAXPROB,MAXAMP,MAXCASE),
     3     long, lat, per(MAXPROB), segModelwt(MAXCASE),
     4     wtoutrisk(100,100,100,100), outrisk(100,100,100,100)       
      real wtmBar(MAXPROB,MAXAMP), wtdBar(MAXPROB,MAXAMP),
     1     wteBar(MAXPROB,MAXAMP), al_segwt(MAXCASE), gmdeag(100),
     2     mindist(MAXCASE), mag(MAXAMP,MAXPROB), magBins(MAXBINS),
     3     dis(MAXAMP,MAXPROB), eps(MAXAMP,MAXPROB), XCostBins(MAXBINS),
     4     ampper(MAXPROB,MAXAMP), DistBins(MAXBINS), EpsBins(MAXBINS) 
      real y(100,100,100) 
      real*8 poisson(MAXAMP), gm(MAXAMP,MAXPROB), test(MAXAMP)
      real*8 wtnEv(MAXPROB,MAXAMP)
      real*8 wtrisk(MAXPROB,MAXAMP,MAXCASE)
      character*80 file1, fname(MAXCASE), filein, fileout, dummy

      write (*,'( 2x,''Enter run file'')')
      read (*,'( a80)') file1
      open (15,file=file1,status='old')

      write (*,*) 'Enter the output filename.'
      read (15,*) fileout
      open (25,file=fileout,status='unknown')
      write (25,*) ' *** Output file from program Haz45-Post *** '
      write (25,*) '              *** Version 45.2 ***           '
      write (25,*) 
      write (25,'(a17,2x,a80)') ' Input filename: ', file1
      write (25,*) 

c     Enter run parameters
      write (*,*) 'Enter number of hazard levels to consider'
      read (15,*) nTest
      write (25,'(i6,2x,a27)') nTest, ' Number of hazard levels = '
      write (25,*)

      write (*,*) 'Enter the hazard levels.'
      read (15,*) (test(i),i=1,nTest)
      write (25,*) 'Hazard Levels:'
      write (25,*) '    AEP         RP (yr)'
      do i=1,nTest
         write (25,'(e12.5,f12.3)') test(i), 1.0/test(i) 
      enddo
      write (25,*)

      write (*,*) 'Enter the number spectral periods.'
      read (15,*) nSpecPer
      write (25,'(i6, 2x,a29)') nSpecPer, 'Number of Spectral Periods = '
      write (25,*) ' Period (sec)    Curve Number     Filename'

C     First compute the equal hazard spectra.
C     Start loop over each spectral period. 
      do 1000 i=1,nSpecPer

c     Specify PSHA Haz45 output files. 
         read (15,'(a80)') filein
         open (10,file=filein,status='old')

         write (*,*)'Enter the corresponding hazard curve number and spectral period.'
         read (15,*) ncurve(i), per(i)

         write (*,*) 'Opening Haz45 ouput file: ', filein
         write (25,'(f8.3,8x,i8,11x,a80)') per(i), ncurve(i), filein 

c     Read in the hazard curves for all periods in given Haz45 ouptut file.
         read (10,*) nflt, testnum 
         read (10,'(a1)') dummy

         do j=1,testnum
            read (10,'(a1)') dummy
         enddo

         do j=1,4
            read (10,'( a1)') dummy
         enddo

         read (10,'( 21x,2f9.3)') long, lat
         read (10,*) nProb

C    Check that requested ncurve is not greater than nProb.
         if (ncurve(i) .gt. nProb) then	
             write (*,*) 'Ncurve value is greater than nProb for given'
             write (*,*) 'Haz45 output file! Check input file parameters.'
             stop 99
         endif  

         do iProb=1,nProb

            read (10,'( 15x,i5)')  iAtten
            read (10,*) nAmp(iProb)
            read (10,'( 61x,30f12.4)') (amp(iProb,k),k=1,nAmp(iProb))

            do l=1,nFlt
               read (10,'( 2x,a38,2f6.3,f8.1,1x,30e12.4)') fname(l),
     1              segModelwt(l),al_segwt(l),
     1              mindist(l),(risk(iProb,k,l),k=1,nAmp(iProb))
            enddo

            read (10,'( 61x,50e12.4)') (nEv(iProb,k),k=1,nAmp(iProb))
            read (10,'(a1)') dummy

            read (10,'( 61x,50e12.3)') (mBar(iProb,k),k=1,nAmp(iProb))
            read (10,'( 61x,50e12.3)') (dBar(iProb,k),k=1,nAmp(iProb))
            read (10,'( 61x,50e12.3)') (eBar(iProb,k),k=1,nAmp(iProb))

            do idum=1,3
               read (10,'( a80)') dummy
            enddo
            

         enddo
         close (10)

C Now set the requested cases for each spectral period. 
C Initialize the arrays to zero for each spectral period.
         do k=1,nFlt
            do l=1,nAmp(ncurve(i))
               wtmBar(i,l) = 0.0
               wtdBar(i,l) = 0.0
               wteBar(i,l) = 0.0
               wtnEv(i,l)  = 0.0
               wtrisk(i,l,k) = 0.0
            enddo
         enddo 

c Set up curves for each given spectral period.
         do l=1,nAmp(ncurve(i))
            wtmBar(i,l)=mbar(ncurve(i),l)
            wtdBar(i,l)=dBar(ncurve(i),l)
            wteBar(i,l)=eBar(ncurve(i),l)
            wtnEv(i,l) = nEv(ncurve(i),l)
            do k=1,nFlt
               wtrisk(i,l,k) = risk(ncurve(i),l,k)
            enddo
            ampper(i,l) = amp(ncurve(i),l)
            nAmpper(i)  = namp(ncurve(i))
         enddo
      
 1000 continue

      
c     Interpolate to desired return period
      do iTest=1,nTest
         do i=1,nSpecPer
		  poisson(1) = 1.-exp(-wtnEv(i,1))
	      do l=2,nAmpper(i)
		     poisson(l) = 1.-exp(-wtnEv(i,l))

c Check for zero values in hazard curve.
                if (poisson(l) .eq. 0. ) then
                   write (*,*) 'Zero Values for hazard curve.'
                   write (*,'( 2x,''Prob = 0'',2x,i5,2f10.4, i5)') l, 
     1                   wtnEv(i,l), test(iTest), nAmpper(i)
                   write (*,*) 'Program will continue but user should check input'
                   write (*,*) 'hazard curve data.'
                   pause
                   write (*,*) 'Hit the Return Key to Continue Program.'
                endif
c Interpolate the hazard curve.
                if ( poisson(l) .lt. test(itest) ) then

                   x = (log(test(itest)) - log(poisson(l-1)))/
     1            (log(poisson(l))-log(poisson(l-1))) 
     1          * (log(ampper(i,l))-log(ampper(i,l-1))) + log(ampper(i,l-1))

                  gm(iTest,i) = exp(x)
C Now interpolate the Mbar, Dbar and epsBar based on the ground motion level.
                mag(iTest,i)=wtmBar(i,l-1) + (wtmBar(i,l)-wtmBar(i,l-1))
     1                      *(gm(itest,i)-ampper(i,l-1))/
     1                       (ampper(i,l)-ampper(i,l-1))
                dis(iTest,i)=wtdBar(i,l-1) + (wtdBar(i,l)-wtdBar(i,l-1))
     1                      *(gm(itest,i)-ampper(i,l-1))/
     1                       (ampper(i,l)-ampper(i,l-1))
                eps(iTest,i)=wteBar(i,l-1) + (wteBar(i,l)-wteBar(i,l-1))
     1                      *(gm(itest,i)-ampper(i,l-1))/
     1                       (ampper(i,l)-ampper(i,l-1))
                   goto 10
                endif
              enddo
 10           continue
         enddo
      enddo

C Write out results
      write (*,*) 
      write (*,'(a19,2x,10f12.6)') ' Testing Levels   :',
     1          (test(iTest),iTest=1,nTest)
      write (*,'(a19,2x,10f12.3)') ' Return Period(yr):',
     1          (1.0/test(iTest),iTest=1,nTest)
      write (*,'(a20)') 'Period (sec)'

      write (25,*)
      write (25,'(a19,2x,10f12.6)') ' Testing Levels   :',
     1          (test(iTest),iTest=1,nTest)
      write (25,'(a19,2x,10f12.3)') ' Return Period(yr):',
     1          (1.0/test(iTest),iTest=1,nTest)
      write (25,'(a20)') 'Period (sec)'
      do i=1,nSpecPer
         write (*,'( 2x,i5,f10.3,4x,10f12.6)') i,per(i),
     1             (gm(iTest,i),iTest=1,nTest)
         write (25,'( 2x,i5,f10.3,4x,10f12.6)') i,per(i),
     1             (gm(iTest,i),iTest=1,nTest)
      enddo

C Now write out the Mbar, DBar, and EpsBar values for each test case. 
      write (25,*) 
      write (25,*) 'Magnitude Bar Results: '
      do i=1,nSpecPer
         write (25,'( 2x,i5,f10.3,4x,10f12.2)') i,per(i),
     1         (mag(iTest,i),iTest=1,nTest)
      enddo
      write (25,*) 
      write (25,*) 'Distance Bar Results: '
      do i=1,nSpecPer
         write (25,'( 2x,i5,f10.3,4x,10f12.3)') i,per(i),
     1         (dis(iTest,i),iTest=1,nTest)
      enddo
      write (25,*) 
      write (25,*) 'Epsilon Bar Results: '
      do i=1,nSpecPer
         write (25,'( 2x,i5,f10.3,4x,10f12.3)') i,per(i),
     1         (eps(iTest,i),iTest=1,nTest)
      enddo
      write (25,*)
      

C     Now compute the interpolated deaggregation matrices based on 
c     computed equal hazard ground motion values from previous steps
c     in the program.       
C     Perform the Equal Deaggregation Analysis
      write (*,*) 
      write (*,*) ' *** Deaggregation Results ***'
      write (*,*) 
      write (25,*) 
      write (25,*) ' *** Deaggregation Results ***'
      write (25,*) 
      write (25,'(a17,2x,a80)') ' Input filename: ', file1
      write (25,*) 
      write (25,*)
      write (25,*) ' Period (sec)    Curve Number     Filename'

C     Read extra spacing line in input file to separate Hazard and Deag output files.
      read (15,'(a80)') dummy
      write (*,'( 2x,''dummy line:'',a80)') dummy
      
C     Start loop over each spectral period. 
      do 2000 i=1,nSpecPer

c     Specify PSHA Haz45 deaggregation output file. 
         read (15,'(a80)') filein
         open (10,file=filein,status='old',err=5000)

c         write (*,*)'Enter the corresponding hazard curve number and spectral period.'
         read (15,*) ncurve(i), per(i)
         write (*,'( ''period:'',f10.3,2x,'' Read hazard out file:'', a80)') per(i), filein

c         write (*,'(a40,2x,a80)') 'Opening Haz45 deaggregation ouput file: ', filein
         write (25,'(f8.3,8x,i8,11x,a40)') per(i), ncurve(i), filein 


c Read in the deaggregation bins.
         read (10,'(i5)') nMagbins
         read (10,'(20f10.3)') (magBins(imag),imag=1,nMagbins)
         read (10,'(i5)') nDistbins
         read (10,'(20f10.3)') (distBins(idis),idis=1,nDistbins)
         read (10,'(i5)') nEpsbins
         read (10,'(20f10.3)') (epsBins(ieps),ieps=1,nepsbins)
         read (10,'(i5)') nXcostbins
         read (10,'(20f10.3)') (xCostBins(iXco),iXco=1,nxCostbins)

	 do idum=1,4
	    read (10,'( a1)') dummy
	 enddo
	 read (10,'( 21x,2f9.3)') long, lat
	 read (10,*) nProb

C Loop over each attenuation relationship.
	 do iProb=1,nProb
	    read (10,*) nAmp(iProb)
	    read (10,'( 46x,30f12.3)') (amp(iProb,k),k=1,nAmp(iProb))

c Skip over bar bins, down to first XCosT values
            do idum=1,11
               read (10,'(a1)') dummy
            enddo
            
c Skip over XCosT values, down to first deaggregation values
            do idum=1,nXcostbins+5
               read (10,'(a1)') dummy
            enddo            

C Read in the deaggregation bins for all epsilon.
            do iMagBin=1,nMagBins-1
               do iDistBin=1,nDistBins-1
                  read (10,'(2x,6f7.1,2x,30(e10.3,2x))') epsBins(1),
     1                 epsBins(nEpsBins),
     1                 magBins(iMagBin),
     1                 magBins(iMagBin+1),
     1                 distBins(iDistBin), distBins(iDistBin+1),
     2                 (outrisk(iMagBin, iDistBin, iProb, iInten),
     2                 iInten=1,nAmp(iProb))
               enddo
            enddo

            read (10,'(a1)') dummy

C Now skip over the bins versus Espilon Bins. 
            nskip = (nMagbins-1)*(nDistBins-1)*(nEpsBins-1)
            do idum=1,nskip+8
               read (10,'(a1)') dummy
            enddo

         enddo
         close (10)

C Now set the requested cases for each spectral period. 

C Initialize the arrays to zero for each spectral period.
c         do iMagbins=1,nMagbins
c            do iDistBins=1,nDistbins
c              do l=1,nAmp(ncurve(i))
c                 wtoutrisk(iMagbins,iDistbins,ncurve(i),l) = 0.0
c              enddo
c            enddo
c         enddo

         do iMagbins=1,nMagbins
            do iDistbins=1,nDistbins
               do l=1,nAmp(ncurve(i))
                  wtoutrisk(iMagbins,iDistbins,i,l)=outrisk(iMagbins,iDistbins,ncurve(i),l)
                  ampper(i,l) = amp(ncurve(i),l)
                  nAmpper(i)  = nAmp(ncurve(i))
               enddo
            enddo
         enddo

 2000 continue
      write (25,*)

C Start loop over each spectral period and perform interpolation of deag bins.
      do itest=1,ntest
         do i=1,nSpecPer

C Enter the ground motion level to analyze.
c 117     write (*,*) 
c         write (*,*) 'Enter the ground motion level for this period.'
c         write (*,*) 'Ground motion value must be between:'
c         write (*,'(2x,a20,f10.3)') 'Min Ground Motion = ', ampper(i,1)
c         write (*,'(2x,a20,f10.3)') 'Max Ground Motion = ', ampper(i,nampper(i))
c         write (*,'(a10,f10.3)') ' Period = ', per(i)
c         read (*,*) gm(i)
            gmdeag(i) = gm(itest,i)

            bflag = 0

C Find the closest two amp levels for interpolation.
         do l=1,nAmpper(i)-1
            if (gmdeag(i) .gt. ampper(i,l) .and. gmdeag(i) .le. ampper(i,l+1) ) then 
               count1 = l
               count2 = l+1
               bflag = 1
               goto 111
            endif
         enddo

         if (bflag .eq. 0 ) then
            write (*,*) 'Ground motion value requested: ', gmdeag(i)
            write (*,*) 'is outside range of amp levels from'
            write (*,*) 'the output file. Ground motion'
            write (*,*) 'level must be between:'
            write (*,'(2f10.3)') ampper(i,1),ampper(i,nampper(i))
            stop 99
         endif

C Now interpolate to the requested ground motion level.
 111     do iMagbins=1,nMagbins
            do iDistbins=1,nDistbins
               y1=wtoutrisk(iMagbins,iDistbins,i,count1)
               y2=wtoutrisk(iMagbins,iDistbins,i,count2)
               x1=ampper(i,count1)
               x2=ampper(i,count2)
               y(iMagbins,iDistbins,i) = y1 + (gmdeag(i)-x1)*(y2-y1)/(x2-x1)

C Set small values to zero.
               if (y(iMagbins,iDistbins,i) .lt. 0.0000001) then
                   y(iMagbins,iDistbins,i) = 0.0
               endif
            enddo
         enddo

C Now write out the results for this period.
         write (25,*) 'Deaggregation Results for Spectral Period:'
         write (25,'(10x,f10.3)') per(i)
         write (25,*) 
         write (25,*) 'Epsilon Bin Range: '
         write (25,'(10x,2f10.3)') epsBins(1),epsBins(nEpsBins)
         write (25,*) 

         write (25,*) '      Hazard Level:'
         write (25,*) '    AEP         RP (yr)'
         write (25,'(f12.8,f12.3)') test(itest), 1.0/test(itest)
         write (25,*)
         write (25,'(a22,f10.5)') 'Ground Motion Level = ', gmdeag(i)
         write (25,*) 
         write (25,'(20x,20(2x,f5.2,'' -'',f5.2))') (Magbins(m),magbins(m+1),
     1           m=1,nMagbins-1)
         do ii=1,nDistbins-1
            write (25,'(2f10.3,20e14.4)') Distbins(ii),Distbins(ii+1),
     1          (y(iMagbins,ii,i),iMagbins=1,nMagbins-1)
         enddo
         write (25,*)
      enddo

      enddo

      write (25,*) 
      write (25,*) 
      write (25,*) ' *** Normal Completion of Haz45-Post Program ***'
      close (25)

      stop 

 5000   write (*,'( 2x,''bad deagg file name'')')
        write (*,'( a80)') filein
       stop

      end
	  
 
	  

