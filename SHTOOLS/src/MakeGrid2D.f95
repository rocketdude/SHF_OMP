subroutine MakeGrid2D(grid, cilm, lmax, interval, nlat, nlong, norm, csphase, f, a, north, south, east, west)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Given the Spherical Harmonic coefficients cilm, this subroutine
!	will compute a 2D grid with equal latitude and longitude spacings. 
!	Note that since this is NOT done using FFTs, this routine is therefore 
!	relatively SLOW! You might want to instead use MakeGridDH.
!
!	The value at grid(1,1) correspons to 90 degrees latitude
!	and 0 degrees longitude, and the longitude spacing goes from 0 to 360, 
!	with both points being calculated. If the optional
!	parameters NORTH, SOUTH, EAST and WEST are specified, the upper-left and lower right coordinates
!	of the output grid are (NORTH, WEST) and (SOUTH, EAST), respectively.
!
!	Calling Parameters:
!		IN
!			cilm:		Spherical harmonic coefficients.
!			lmax:		Maximum degree of expansions to be performed.
!			interval:	Spacing of output grid in DEGREES.
!		OUT
!			grid:		Gridded expansion of spherical harmonic coefficients.
!			nlat:		Number of latitude points for the grid.
!			nlong:		Number of longitude points for the grid.
!		OPTIONAL (IN)
!			norm		Spherical harmonic normalization:
!						(1) "geodesy" (default)
!						(2) Schmidt
!						(3) unnormalized
!						(4) orthonormalized
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!			f, a:		Flattening and semimajor axis of the function. If included, 
!					an ellipsoid with these parameters will be subtracted from the
!					data.
!			north		Maximum latitude to compute, in degrees.
!			south		Minimum latitude to compute, in degrees.
!			east		Maximum longitude to compute, in degrees.
!			west		Minimum latitude to compute, in degrees.
!			
!
!	Notes:
!		1.	If lmax is greater than the the maximum spherical harmonic
!			degree of the input file, then this file will be ZERO PADDED!
!			(i.e., those degrees after lmax are assumed to be zero).
!
!	Dependencies:		PlmBar, PlBar, PlmSchmidt, PlSchmidt, PLegendreA, PLegendre, 
!				PlmON, CSPHASE_DEFAULT
!
!	Written by Mark Wieczorek (2003)
!
!	August 2005 - 	sines and cosines are now precomputed, leading to a savings in time
!			by a factor of up to 15 (if the code is compiled with normal
!			optimizations.
!	April 12, 2007	Added option to allow calculation of function referenced to a flattened ellipsoid.
!	July 21, 2007 	Reference ellipsoid is now calculated exactly using the flattening and semi-major
!			axis instead of just the second degree Legendre expansion coefficient.
!	August 19, 2007	Added option to calculate subdomains of a raster grid.
!
!	Copyright (c) 2005-2007, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: PlmBar, PlBar, PlmSchmidt, PlSchmidt, PLegendreA, PLegendre, &
				PlmON, PlON, CSPHASE_DEFAULT

	implicit none

	real*8, intent(in) :: 	cilm(:,:,:), interval
	real*8, intent(out) :: 	grid(:,:)
	integer, intent(in) :: 	lmax
	integer, intent(out) :: nlat, nlong
	integer, intent(in), optional :: norm, csphase
	real*8, intent(in), optional :: f, a, north, south, east, west
	integer :: 		l, m, j, k, index, l1, m1, lmax_comp, phase, lnorm, temp, astat(4)
	real*8 :: 		pi, latmax, latmin, longmin, longmax, lat, longitude, &
				x, intervalrad,  r_ref
	real*8, allocatable ::	pl(:), cosm(:, :), sinm(:, :), cilm2(:,:, :)
	
	temp = 0
	if (present(north)) temp = temp + 1
	if (present(south)) temp = temp + 1
	if (present(east)) temp = temp + 1
	if (present(west)) temp = temp + 1
	
	if (temp /=0 .and. temp /=4) then
		print*, "Error --- MakeGrid2d"
		print*, "The optional parameters NORTH, SOUTH, EAST, and WEST must all be specified", &
			present(north), present(south), present(east), present(west)
		stop
	endif
	
	if (temp == 4) then
		latmax = north
		latmin = south
		longmin = west
		longmax = east
		
		if (latmax < latmin) then
			print*, "Error --- MakeGrid2d"
			print*, "NORTH must be larger than SOUTH."
			print*, "NORTH = ", latmax
			print*, "SOUTH = ", latmin
			stop
		endif
		
		if (latmax > 90.0d0 .or. latmin < -90.0d0) then
			print*, "Error --- MakeGrid2d"
			print*, "NORTH and SOUTH must lie between 90 and -90."
			print*, "NORTH = ", latmax
			print*, "SOUTH = ", latmin
			stop
		endif
		
		if (longmin > longmax) longmax = longmax + 360.0d0
		
	else
		latmax = 90.0d0
		latmin = -90.0d0
		longmin = 0.0d0
		longmax = 360.d0
	endif
	
	nlat = (latmax - latmin) / interval + 1
	nlong = (longmax - longmin) / interval + 1
	
	if (present(norm)) then
		if (norm > 4 .or. norm < 1) then
			print*, "Error - MakeGrid2d"
			print*, "Parameter NORM must be 1 (geodesy), 2 (Schmidt), 3 (unnormalized), or 4 (orthonormalized)."
			print*, "Input value is ", norm
			stop
		endif
		
		lnorm = norm
		
	else 
		lnorm = 1
	endif
	
 	if (present(csphase)) then
     		if (csphase /= -1 .and. csphase /= 1) then
     			print*, "MakeGrid2D --- Error"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)"
     			print*, "Input valuse is ", csphase
     			stop
     		else
     			phase = csphase
     		endif
     	else
     		phase = CSPHASE_DEFAULT
     	endif

    	
 	if (size(grid(:,1)) < nlat .or. size(grid(1,:)) < nlong ) then
		print*, "Error --- MakeGrid2D"
		print*, "GRID must be dimensioned ( (LATMAX-LATMIN)/INTERVAL+1, (LONGMAX-LONGMIN)/INTERVAL+1 ) where"
		print*, "INTERVAL = ", interval
		print*, "LATMAX = ", latmax
		print*, "LATMIN = ", latmin
		print*, "LONGMIN = ", longmin
		print*, "LONGMAX = ", longmax
		print*, "Input array is dimensioned ", size(grid(:,1)), size(grid(1,:))
		stop
	elseif (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
		print*, "Error --- MakeGrid2D"
		print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	if ((present(f) .and. .not. present(a)) .or. (present(a) .and. .not. present(f)) ) then
		print*, "Error --- MakeGrid2D"
		print*, "Both F and A must be present"
		print*, "F ", present(f)
		print*, "A ", present(a)
		stop
	endif
	
	allocate(pl((lmax+1)*(lmax+2)/2), stat = astat(1))
	allocate(cosm(int(360./interval + 1), lmax+1), stat = astat(2))
	allocate(sinm(int(360./interval + 1), lmax+1), stat = astat(3))
	allocate(cilm2(2,lmax+1,lmax+1), stat = astat(4))
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /=0 .or. astat(4) /= 0) then
		print*, "Error --- MakeGrid2D"
		print*, "Problem allocating arrays PL, COSM, SINM, and CILM2", astat(1), astat(2), astat(3), astat(4)
		stop
	endif
	
	pi = acos(-1.0d0)
	grid = 0.0d0
	
	lmax_comp = min(lmax, size(cilm(1,1,:))-1)

	intervalrad = interval*pi/180.0d0
	
	cilm2(1:2,1:lmax+1,1:lmax+1) = cilm(1:2,1:lmax+1,1:lmax+1)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Precomputing sines and cosines leads to an increase in speed by a factor
	!	of almost 4 with no optimization, and by a factor of about 15 with normal 
	!	optimizations.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do k=1, nlong
	
		longitude = longmin*pi/180.0d0 + dble(k-1)*intervalrad
		do m=0, lmax
			cosm(k, m+1) = cos(m*longitude)
			sinm(k, m+1) = sin(m*longitude)
		enddo
		
	enddo

	do j=1, nlat
	
		lat = latmax - dble(j-1)*interval
		x = sin(lat*pi/180.0d0)
		
		if (lat == 90.0d0 .or. lat == -90.0d0) then
		
			select case(lnorm)
				case(1); call PlBar(pl, lmax_comp, x)
				case(2); call PlSchmidt(pl, lmax_comp, x)
				case(3); call PLegendre(pl, lmax_comp, x)
				case(4); call PlON(pl, lmax_comp, x)
			end select
	
			do l = lmax_comp, 0, -1
				l1 = l+1
				grid(j,1) = grid(j,1) + cilm2(1,l1,1) * pl(l1)
			enddo
			
			grid(j, 2:nlong) = grid(j,1)
	
			if (present(f)) grid(j,1:nlong) = grid(j,1:nlong) - a*(1.0d0-f)
		
		else
		
			select case(lnorm)
				case(1); call PlmBar(pl, lmax_comp, x, csphase = phase)
				case(2); call PlmSchmidt(pl, lmax_comp, x, csphase = phase)
				case(3); call PLegendreA(pl, lmax_comp, x, csphase = phase)
				case(4); call PlmON(pl, lmax_comp, x, csphase = phase)
			end select
			
			do k = 1, nlong
			
				do l = lmax_comp, 0, -1
			
					l1 = l+1
					index = (l+1)*l/2 + 1
					grid(j,k) = grid(j,k) + cilm2(1,l1,1) * pl(index)
				
					do m = 1, l, 1
						m1 = m+1
						index = (l+1)*l/2 + m + 1
						grid(j,k) = grid(j,k) + ( cilm2(1,l1,m1)*cosm(k,m+1) + &
							cilm2(2,l1,m1)*sinm(k,m+1) ) * pl(index)
					enddo
				enddo
			
			enddo
		
			if (present(f)) then
				r_ref = a**2 * (1.0d0 + tan(lat*pi/180.0d0)**2) / &
					(1.0d0  + tan(lat*pi/180.0d0)**2 / (1.0d0 - f)**2 )
				r_ref = sqrt(r_ref)
				grid(j,1:nlong) = grid(j,1:nlong) - r_ref
			endif
			
		endif
		
	enddo
		
	
	! deallocate memory
	select case(lnorm)
		case(1); call PlmBar(pl, -1, x, csphase = phase)
		case(2); call PlmSchmidt(pl, -1, x, csphase = phase)
		case(4); call PlmON(pl, -1, x, csphase = phase)
	end select
	deallocate(pl)
	deallocate(cosm)
	deallocate(sinm)
	deallocate(cilm2)
	
end subroutine MakeGrid2D


