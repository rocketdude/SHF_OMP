program SHMag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will read in a set of magnetic spherical harmonic coefficients
!	and will expand these into gridded raster ascii files corresponding
!	to the phi, theta, and radial components, plus the total field. The field is calculated
!	on a speroid with mean radius r and flattening f.
!
!	The included spherical harmonic file FSU_mars90.sh is the martian magnetic
!	field model of Cain et al., 2003.
!
!	Dependencies:	SHRead
!			MakeMagGrid2D
!
!	Written by Mark Wieczorek, February 2004
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS

	implicit none
	integer, parameter ::	maxdeg = 500
	character*120 ::		infile, radf, thetaf, phif, totalf
	real*8 ::		glm(2,maxdeg+1, maxdeg+1), header(4), interval, r0, r, &
				f, mpr, z, timein, timeout
	real*8, allocatable ::	rad(:,:), phi(:,:), theta(:,:), total(:,:)		
	integer ::		lmax, nlong, nlat, i, j, nlat1, nlong1, astat(4)
	
	
	f = 1.0d0/169.864881d0	! Mars flattening = (R_eq - R_p)/R_eq
	mpr = 3389.508d3	! Mean radius of mars
	z = 145.d3		! mean altitude to calculate field at.
	
	infile = "../ExampleDataFiles/FSU_mars90.sh"
	
	call SHRead(infile, glm, lmax, header=header(1:4), skip=1)
	r0 = header(1)*1.d3
	print*, "Reference radius of spherical harmonic coefficients R0 (km) = ", r0/1.d3
	print*, "Lmax of spherical harmonic model = ", lmax
	
	r = mpr + z
	print*, "Field calculated on a flattened ellipsoid with mean radius R (km) = ", r/1.d3
	print*, "Flattening = ", f
		
	print*, "Interval for gridded data (degrees) > "
	read(*,*) interval
	
	nlat1 = 180.0d0 / interval + 1
	nlong1 = 360.0d0 / interval + 1
	
	allocate(rad(nlat1,nlong1), stat = astat(1))
	allocate(theta(nlat1,nlong1), stat = astat(2))
	allocate(phi(nlat1,nlong1), stat = astat(3))
	allocate(total(nlat1,nlong1), stat = astat(4))
	
	if (sum(astat(1:4)) /= 0) then
		print*, "Problem allocating RAD, THETA, PHI and TOTAL", astat(1), astat(2), astat(3), astat(4)
		stop
	endif
		
	radf = "radial_145f.dat"
	thetaf = "theta_145f.dat"
	phif = "phi_145f.dat"
	totalf = "total_145f.dat"
	
	open(12,file=radf)
	open(13,file=phif)
	open(14,file=thetaf)
	open(15,file=totalf)
	
	call cpu_time(timein)
	
	call MakeMagGrid2D(rad, phi, theta, total, glm, r0, r, f, lmax, interval, nlat, nlong)
	
	call cpu_time(timeout)
	
	print*, "Elapsed time (sec) = ", timeout-timein
	
	print*, "Maximum and minimum intensity (nT) = ", maxval(total(1:nlat,1:nlong)), minval(total(1:nlat,1:nlong))
	
	print*, nlat, nlong
	
	! write(12,*) nlat, nlong
	! write(13,*) nlat, nlong
	! write(14,*) nlat, nlong
	! write(15,*) nlat, nlong
	
	do i=1, nlat
		do j=1, nlong
			write(12,*) rad(i,j)
			write(13,*) phi(i,j)
			write(14,*) theta(i,j)
			write(15,*) total(i,j)
		enddo
	enddo
	
	close(12)
	close(13)
	close(14)
	close(15)
	
	deallocate(rad)
	deallocate(theta)
	deallocate(phi)
	deallocate(total)
	
end program SHMag