subroutine SHReadJPL(filename, cilm, lmax, error, gm, formatstring)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will read the file format of spherical harmonic
!	coefficients sometimes used by the JPL group into standard arrays for
!	of Cilm and the error. In order to for this to work, you need to know
!	a priori the maximum spherical harmonic degree of the file.
!	
!	The input file contains
!		1. comment lines starting with "#"
!		2. GM (optional)
!		3. A list of J_l, which is -C(1,l+1,1)
!		4. A list of the cosine and sine terms
!		5. Starting from 2, the same thing over for the errors.
!
!	Calling Parameters
!		IN
!			filename: 	The name of the file.
!			lmax:		Maximum spherical harmonic degree.
!		OUT
!			cilm:		An array of the spherical harmonic coefficients.
!		OPTIONAL
!			error:	An array containing the error coefficients.
!			formatstring:	This is a string containing an I/O specification
!					for the numbers of the spherical harmonic coefficients.
!					For MGNP180U.ODP this is e19.12 (default). For topo_gtdr32_360, 
!					this is e17.10. Note that this assumes that there are 
!					two spaces after the equal sign.
!
!	Written by Mark Wieczorek (September 2005)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	
	character(*), intent(in) ::	filename
	integer, intent(in) ::	lmax
	real*8, intent(out) ::	cilm(:,:,:)
	real*8, intent(out), optional :: error(:,:,:), gm(2)
	character, intent(in), optional :: formatstring*6
	
	real*8 :: gm1, gm2
	logical ::	gmpresent	
	integer l, m, stat, i, ll1, mm1, ll2, mm2, skip
	character :: c*4, s*4, j*4, dumb*14, dumb2*14, js*4, cs*4, ss*4, inum*2
	
	if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- SHReadJPL"
		print*, "CILM must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif 
	
	if (present(error)) then
		if (size(error(:,1,1)) < 2 .or. size(error(1,:,1)) < lmax+1 .or. size(error(1,1,:)) < lmax+1) then
			print*, "Error --- SHReadJPL"
			print*, "error must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", lmax
			print*, "Input array is dimensioned ", size(error(:,1,1)), size(error(1,:,1)), size(error(1,1,:))
			stop
		endif
	endif
	
	
	gmpresent=.false.
	
	c = "OBAC";	cs = "SIGC"
	s = "OBAS";	ss = "SIGS"
	j = "OBAJ";	js = "SIGJ"
	
	if (lmax <= 99) then
		 inum = "i2"
	elseif (lmax >= 100 .and. lmax <=999) then
		inum = "i3"
	else
		inum = "i4"
	endif
	
	
	open(13,file=filename)
	
	! Skip past header lines
	
	skip = 0
	
	do 
		read(13,"(a)", iostat=stat) dumb
		if (stat /= 0) then
			print*, "Error --- SHReadJPL"
			print*, "Problem reading header lines of ", filename
			stop
		endif
		
		if (dumb(1:1) == "#") then
			skip = skip + 1
		else
			exit
		endif
	enddo
	
	! Try to read GM, if present
	
	rewind(13)
	
	do i=1, skip
		read(13, *, iostat=stat)
	enddo
	
	read(13, "(a12, f14.6)", iostat=stat) dumb, gm1

	if (stat/=0) then
		print*, "Error --- SHReadJPL"
		print*, "Problem reading GM(1)"
		stop
	endif
	
	if (dumb(2:3) == "GM") then
		if ( present(gm) ) then
			gm(1) = gm1
			gmpresent=.true.
		endif
	else
		rewind(13)
		do i=1, skip
			read(13,*)
		enddo
	endif
	
	! Start reading J_l coefficients
		
	do i=1, lmax
		if (present(formatstring)) then
			read(13,"(1x, a5, "//inum//", 5x,"//formatstring//")", iostat=stat) dumb, l, cilm(1,i+1,1)	
		else
			read(13,"(1x, a5, "//inum//", 5x, e19.12)", iostat=stat) dumb, l, cilm(1,i+1,1)	
		endif
		
		if (stat /=0) then
			print*, "Error --- SHReadJPL"
			print*, "Problem reading file", i
			stop
		elseif (dumb(1:4) /= j) then
			print*, "Error --- SHReadJPL"
			print*, "Problem with line specifiers in file ", filename
			print*, "Read sring = ", dumb(1:4)
			print*, "Expected string = ", j
			stop
		elseif (i /=l) then
			print*, "Error --- SHReadJPL"
			print*, "Problem with indices in file ", filename
			print*, "Read indice = ", l
			print*, "Expected indice = ", i
			stop
		endif
		
		cilm(1,i+1,1) = -cilm(1,i+1,1)
		
	enddo

	! Read the Clm and Slm coefficients
	
	do l=1, lmax
		do m=1, l
			if (present(formatstring)) then
				read(13, &
					"(1x, a5, "//inum//", 1x, "//inum//", 5x,"//formatstring//", 5x, a5, "//inum//", 1x, "//inum//", 5x,"//formatstring//")", &
					iostat=stat) dumb, ll1, mm1, cilm(1,l+1,m+1), dumb2, ll2, mm2, cilm(2,l+1,m+1)
			else
				read(13,"(1x, a5, "//inum//", 1x, "//inum//", 5x, e19.12, 5x, a5, "//inum//", 1x, "//inum//", 5x, e19.12)", &
					iostat=stat) dumb, ll1, mm1, cilm(1,l+1,m+1), dumb2, ll2, mm2, cilm(2,l+1,m+1)
			endif
			
			if (stat /=0) then
				print*, "Error --- SHReadJPL"
				print*, "Problem reading file", l, m
				stop
			elseif (dumb(1:4) /= c .or. dumb2(1:4) /=s) then
				print*, "Error --- SHReadJPL"
				print*, "Problem with line specifiers in file ", filename
				print*, "Read srings = ", dumb(1:4), dumb2(1:4)
				print*, "Expected strings = ", c, s
				stop
			elseif (ll1/=l .or. ll2/=l .or. mm1/=m .or. mm2/=m ) then
				print*, "Error --- SHReadJPL"
				print*, "Problem with indices in file ", filename
				print*, "Read indices (l1, m1), (l2, m2) = ", ll1, mm1, ll2, mm2
				print*, "Expected indices (l, m) = ", l, m
				stop
			endif
			
		enddo
	enddo
	
	! Next read uncertainties
	
	if (present(error)) then
		if (gmpresent) then
			read(13, "(a12, f17.6)", iostat=stat) dumb, gm2
			gm(2) = gm2
			if (stat/=0) then
				print*, "Error --- SHReadJPL"
				print*, "Problem reading GM(2)"
				stop
			endif
		endif
		
		do i=1, lmax
			if (present(formatstring)) then
				read(13,"(1x, a5, "//inum//", 5x,"//formatstring//")", iostat=stat) dumb, l, error(1,i+1,1)	
			else
				read(13,"(1x, a5, "//inum//", 5x, e19.12)", iostat=stat) dumb, l, error(1,i+1,1)	
			endif

			if (stat /=0) then
				print*, "Error --- SHReadJPL"
				print*, "Problem reading file during J2 error coefficients", i
				stop
			elseif (dumb(1:4) /= js) then
				print*, "Error --- SHReadJPL"
				print*, "Problem with line specifiers in file ", filename
				print*, "Read sring = ", dumb(1:4)
				print*, "Expected string = ", js
				stop
			elseif (i /=l) then
				print*, "Error --- SHReadJPL"
				print*, "Problem with indices in file ", filename
				print*, "Read indice = ", l
				print*, "Expected indice = ", i
				stop
			endif

		enddo
			
		do l=1, lmax
			do m=1, l
			if (present(formatstring)) then
				read(13, &
				"(1x, a5, "//inum//", 1x, "//inum//", 5x,"//formatstring//", 5x, a5, "//inum//", 1x, "//inum//", 5x,"//formatstring//")", &
					iostat=stat) dumb, ll1, mm1, error(1,l+1,m+1), dumb2, ll2, mm2, error(2,l+1,m+1)
			else
				read(13,"(1x, a5, "//inum//", 1x, "//inum//", 5x, e19.12, 5x, a5, "//inum//", 1x, "//inum//", 5x, e19.12)", &
					iostat=stat) dumb, ll1, mm1, error(1,l+1,m+1), dumb2, ll2, mm2, error(2,l+1,m+1)
			endif
			
			if (stat /=0) then
				print*, "Error --- SHReadJPL"
				print*, "Problem reading file during error coefficients", l, m
				stop
			elseif (dumb(1:4) /= cs .or. dumb2(1:4) /=ss) then
				print*, "Error --- SHReadJPL"
				print*, "Problem with line specifiers in file ", filename
				print*, "Read srings = ", dumb(1:4), dumb2(1:4)
				print*, "Expected strings = ", cs, ss
				stop
			elseif (ll1/=l .or. ll2/=l .or. mm1/=m .or. mm2/=m ) then
				print*, "Error --- SHReadJPL"
				print*, "Problem with indices in file ", filename
				print*, "Read indices (l1, m1), (l2, m2) = ", ll1, mm1, ll2, mm2
				print*, "Expected indices (l, m) = ", l, m
				stop
			endif
			
			enddo
		enddo
	endif
	
	close(13)
	
end subroutine SHReadJPL
