Program TimingAccuracyGLQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will test the accuracy of the spherical harmonic GLQ tranformation
!	routines by expanding a field in the space domain, transforming this
!	to spherical harmonics, and comparing the relative error of the coefficients
!
!
!	Written by Mark Wieczorek (July 2006)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS
	
	implicit none	
		
	integer, parameter ::	maxdeg = 2800
	character*80 ::		outfile1, outfile2, outfile3, outfile4, outfile
	real*8 ::		cilm(2,maxdeg+1,maxdeg+1), cilm2(2,maxdeg+1,maxdeg+1), &
				zero(maxdeg+1), gridglq(maxdeg+1,2*maxdeg+1), w(maxdeg+1), &
				huge8, maxerror, err1, err2, beta, rms, timein(3), timeout(3)
	integer ::		lmax, l, m, seed
	
	print*, "Value of beta for power law Sff = l^(-beta) > "
	read(*,*) beta
	
	print*, "output file name > "
	read(*,*) outfile
	
	outfile1 = trim(outfile) // ".timef"
	outfile2 = trim(outfile) // ".timei"
	outfile3 = trim(outfile) // ".maxerror"
	outfile4 = trim(outfile) // ".rmserror"
		
	huge8 = huge(maxerror)
	
	seed = -1053253
	
	cilm = 0.0d0
	
	do l=1, maxdeg
		
		do m=0, l
			if (m==0) then
				cilm(1,l+1,m+1) = RandomGaussian(seed)
			else
				cilm(1,l+1,m+1) = RandomGaussian(seed)
				cilm(2,l+1,m+1) = RandomGaussian(seed)
			endif
		enddo
		cilm(1:2, l+1, 1:l+1) = cilm(1:2, l+1, 1:l+1)*sqrt(dble(l)**beta)/sqrt(2.0d0*l+1)
	enddo	
	
	print*, "Lmax, Maximum rel. error of Cilm, RMS relative error, Precompute time (sec), Time inverse (sec), Time forward (sec)"

	do lmax=100, maxdeg, 100
		
		if (lmax==100) then
			open(12,file=outfile1, status="replace")
			open(13,file=outfile2, status="replace")
			open(14,file=outfile3, status="replace")
			open(15,file=outfile4, status="replace")
		else
			open(12,file=outfile1, position="append")
			open(13,file=outfile2, position="append")
			open(14,file=outfile3, position="append")
			open(15,file=outfile4, position="append")
		endif
		
		call cpu_time(timein(1))
		call PreCompute(lmax, zero(1:lmax+1), w(1:lmax+1))
		call cpu_time(timeout(1))
						
		call cpu_time(timein(2))	
		call MakeGridGLQ(gridglq(1:lmax+1, 1:2*lmax+1), cilm(1:2,1:lmax+1, 1:lmax+1), lmax, zero=zero(1:lmax+1))
		call cpu_time(timeout(2))
				
		call cpu_time(timein(3))
		call SHExpandGLQ(cilm2(1:2, 1:lmax+1, 1:lmax+1), lmax, gridglq(1:lmax+1, 1:2*lmax+1), w, zero=zero(1:lmax+1))
		call cpu_time(timeout(3))
		
		maxerror = 0.0d0
		rms = 0.0d0
		
		do l=1, lmax
			
			do m=0, l
				if (m==0) then
					err1 = abs( (cilm(1,l+1,m+1)-cilm2(1,l+1,m+1)) / cilm(1,l+1,m+1) )
					if (err1 >= maxerror) maxerror = err1
					rms = rms + err1**2
				else
					err1 = abs( (cilm(1,l+1,m+1)-cilm2(1,l+1,m+1)) / cilm(1,l+1,m+1) )
					err2 = abs( (cilm(2,l+1,m+1)-cilm2(2,l+1,m+1)) / cilm(2,l+1,m+1) )
					if (err1 >= maxerror) maxerror = err1
					if (err2 >= maxerror) maxerror = err2
					rms = rms + err1**2 + err2**2
				endif
			enddo
		enddo
		rms = sqrt(rms/dble(l+1)**2)
		
		! elasped time in seconds!
		print*, lmax, maxerror, rms, timeout(1)-timein(1), timeout(2)-timein(2), timeout(3)-timein(3)
		write(12,*) lmax, timeout(2)-timein(2)
		write(13,*) lmax, timeout(3)-timein(3)
		write(14,*) lmax, maxerror
		write(15,*) lmax, rms
		
		if (maxerror > huge8) then
			print*, "Overflow problems"
			print*, "Overflow problems"
			close(12)
			close(13)
			close(14)
			close(15)
			stop
		endif
		
		if (maxerror > 10.0d0) then
			close(12)
			close(13)
			close(14)
			close(15)
			stop
		endif
		
		close(12)
		close(13)
		close(14)
		close(15)
			
	enddo

end program TimingAccuracyGLQ
