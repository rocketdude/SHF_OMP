subroutine MakeMagGridDH(cilm, lmax, r0, a, f, rad_grid, theta_grid, phi_grid, total_grid, n, sampling, lmax_calc, pot_grid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Given the Schmidt normalized magnetic potential spherical harmonic coefficients CILM, this subroutine
!	will compute a 2D Driscoll and Healy sampled grid of the three components and magnitude field. 
!	The output grids are cacluated on a flattened ellipsoid with semi-major axis A and flattening F. 
!	The output grids contain N samples in latitude and longitude by default, but if the optional parameter 
!	SAMPLING is set to 2, the grids will contain N samples in latitude and 2N samples in longitude. 
!	In order to calculate the entire gravitational acceleration, it is necessary that the degree-0 term 
!	be set equal to 1. The radial magnetic field component is assumed to be positive when directed UPWARDS. 
!	The magnetic potential is defined according to
!	
!		V = R0 Sum_{l=1}^LMAX (R0/r)^{l+1} Sum_{m=-l}^l C_{lm} Y_{lm}.
!
!	and
!
!		B = - Grad V
! 
!	The first latitudinal band of the grid corresponds to 90 N, the latitudinal band for 90 S is not calculated, 
!	and the latitudinal sampling interval is 180/N degrees. The first longitudinal band is 0 E, the longitudinal 
!	band for 360 E is not calculated, and the longitudinal sampling interval is 360/N for equally sampled and 180/N 
!	for equally spaced grids, respectively.
!
!	This routine assumes that the spherical harmonic coefficients of geodesy 4-pi normalized.
!
!	Calling Parameters:
!		IN
!			cilm:		Magnetic potential spherical harmonic coefficients.
!			lmax:		The maximum spherical harmonic degree of the function, used to 
!					determines the number of samples N.
!			r0:		Reference radius of potential coefficients.
!			a:		The semimajor axis of the flattened ellipsoid.
!			f:		Flattening of the planet.
!		IN, OPTIONAL
!			sampling	(1) Grid is N latitudes by N longitudes (default).
!					(2) Grid is N by 2N. The higher frequencies resulting
!					from this oversampling in longitude are discarded, and hence not
!					aliased into lower frequencies.
!                       lmax_calc	The maximum spherical harmonic degree to evaluate
!					the coefficients up to.
!		OUT
!			rad_grid	Gridded expansion of the radial component of the 
!					magnetic field.
!			theta_grid	Gridded expansaion of the theta component of the 
!					gravitational field.
!			phi_grid	Gridded expansaion of the phi component of the 
!					gravitational field.
!			total_grid	Gridded expansaion of the the magnitude of the 
!					gravitational field.
!			N		Number of samples in latitude. Number of samples in longitude is N when sampling
!					is 1 (default), and is 2N when sampling is 2.
!		OUT, OPTIONAL
!			pot_grid	Magnetic potential on the ellipsoid in SI units.
!
!	Notes:
!		1.	If lmax is greater than the the maximum spherical harmonic
!			degree of the input file, Cilm will be ZERO PADDED!
!			(i.e., those degrees after lmax are assumed to be zero).
!		2. 	Latitude is geocentric latitude.
!
!	Dependencies:	FFTW3
!
!	Written by Mark Wieczorek 2012.
!
!	Copyright (c) 2012, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use FFTW3
	
	implicit none
	
	real*8, intent(in) :: 	cilm(:,:,:), r0, a, f
	real*8, intent(out) ::	rad_grid(:,:), theta_grid(:,:), phi_grid(:,:), total_grid(:,:)
	real*8, intent(out), optional :: pot_grid(:,:)
	integer, intent(in) :: 	lmax
	integer, intent(out) ::	n
	integer, intent(in), optional :: sampling, lmax_calc
	integer :: 		l, m, i, l1, m1, lmax_comp, i_eq, i_s, astat(4), nlong
	real*8 :: 		grid(4*lmax+4), pi, theta, scalef, rescalem, u, p, dp, pmm, &
				pm1, pm2, z, tempr, r_ex, lat, prefactor(lmax), coefr0, coefu0, coefrs0, &
				coeft0, coefts0, coefp0, coefps0, coefus0
	complex*16 ::		coef(2*lmax+3), coefr(2*lmax+3), coefrs(2*lmax+3), &
				coeft(2*lmax+3), coefts(2*lmax+3), coefp(2*lmax+3), coefps(2*lmax+3), &
				coefu(2*lmax+3), coefus(2*lmax+3), tempc
	integer*8 ::		plan
	real*8, save, allocatable ::	ff1(:,:), ff2(:,:), sqr(:)
	integer*1, save, allocatable ::	fsymsign(:,:)
	integer, save ::	lmax_old=0
	logical ::		calcu
	
	n = 2*lmax+2
	
	if (present(sampling)) then
		if (sampling /= 1 .and. sampling /=2) then
			print*, "Error --- MakeMagGridDH"
			print*, "Optional parameter SAMPLING must be 1 (N by N) or 2 (N by 2N)."
			print*, "Input value is ", sampling
			stop
		endif
	endif
	
	if (size(cilm(:,1,1)) < 2) then
		print*, "Error --- MakeMagGridDH"
		print*, "CILM must be dimensioned as (2, *, *)."
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif 
	
	if (present(sampling)) then
		if (sampling == 1) then
			nlong = n
		else
			nlong = 2*n
		endif
	else
		nlong = n
	endif
	
	if (size(rad_grid(:,1)) < n .or. size(rad_grid(1,:)) < nlong .or. size(theta_grid(:,1)) < n .or. size(theta_grid(1,:)) < nlong &
		.or. size(phi_grid(:,1)) < n .or. size(phi_grid(1,:)) < nlong .or. size(total_grid(:,1)) < n .or. &
		size(total_grid(1,:)) < nlong) then
		print*, "Error --- MakeMagGridDH"
		if (present(sampling)) then
			if (sampling == 1) then
				print*, "RAD_GRID, THETA_GRID, PHI_GRID, and TOTAL_GRID must be dimensioned as (N, N) where N is ", n
			elseif (sampling == 2) then
				print*, "RAD_GRID, THETA_GRID, PHI_GRID, and TOTAL_GRID must be dimensioned as (N, 2N) where N is ", n
			endif
		else
			print*, "RAD_GRID, THETA_GRID, PHI_GRID, and TOTAL_GRID must be dimensioned as (N, N) where N is ", n
		endif
		print*, "Input dimensions are ", size(rad_grid(:,1)), size(rad_grid(1,:)), size(theta_grid(:,1)), size(theta_grid(1,:)), &
			size(phi_grid(:,1)), size(phi_grid(1,:)), size(total_grid(:,1)), size(total_grid(1,:))
		stop
	endif
	
	if (present(pot_grid)) then
		calcu = .true.
		if (size(pot_grid(:,1)) < n .or. size(pot_grid(1,:)) < nlong) then
			print*, "Error --- MakeMagGridDH"
			if (present(sampling)) then
				if (sampling == 1) then
					print*, "POT_GRID must be dimensioned as (N, N) where N is ", n
				elseif (sampling == 2) then
					print*, "POT_GRID must be dimensioned as (N, 2N) where N is ", n
				endif
			else
				print*, "POT_GRID must be dimensioned as (N, N) where N is ", n
			endif
			print*, "Input dimensions are ", size(pot_grid(:,1)), size(pot_grid(1,:))
			stop
		endif
	else 
		calcu = .false.
	endif
			
	pi = acos(-1.0d0)
	
	scalef = 1.0d-280
	
	if (present(lmax_calc)) then
		if (lmax_calc > lmax) then
			print*, "Error --- MakeMagGridDH"
			print*, "LMAX_CALC must be less than or equal to LMAX."
			print*, "LMAX = ", lmax
			print*, "LMAX_CALC = ", lmax_calc
			stop
		else
			lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1, lmax_calc)
		endif
	else
		lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1)
	endif
	
	if (lmax_comp == 0) then
		rad_grid(1:n,1:nlong) = 0.0d0
		theta_grid(1:n,1:nlong) = 0.0d0
		phi_grid(1:n,1:nlong) = 0.0d0
		total_grid(1:n,1:nlong) = 0.0d0
		if (calcu) pot_grid(1:n,1:nlong) = 0.0d0
		return
	endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate recursion constants used in computing the Legendre polynomials
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (lmax_comp /= lmax_old) then
		
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(ff1)) deallocate(ff1)
		if (allocated(ff2)) deallocate(ff2)
		if (allocated(fsymsign)) deallocate(fsymsign)
		
		allocate(sqr(2*lmax_comp+1), stat=astat(1))
		allocate(ff1(lmax_comp+1,lmax_comp+1), stat=astat(2))
		allocate(ff2(lmax_comp+1,lmax_comp+1), stat=astat(3))
		allocate(fsymsign(lmax_comp+1,lmax_comp+1), stat=astat(4))
		
		if (sum(astat(1:4)) /= 0) then
			print*, "MakeMagGridDH --- Error"
			print*, "Problem allocating arrays SQR, FF1, FF2, or FSYMSIGN", astat(1), astat(2), astat(3), astat(4)
			stop
		endif
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Calculate signs used for symmetry of Legendre functions about equator. For 
		!	the first derivative in theta, these signs are reversed.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		do l = 0, lmax_comp, 1
			do m = 0, l, 1
				if (mod(l-m,2) == 0) then
					fsymsign(l+1,m+1) = 1
				else
					fsymsign(l+1,m+1) = -1
				endif
			enddo
		enddo
			
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		!	Precompute square roots of integers that are used several times.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		do l=1, 2*lmax_comp+1
			sqr(l) = sqrt(dble(l))
		enddo

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Precompute multiplicative factors used in recursion relationships
		! 		P(l,m) = x*f1(l,m)*P(l-1,m) - P(l-2,m)*f2(l,m)
		!		k = l*(l+1)/2 + m + 1
		!	Note that prefactors are not used for the case when m=l as a different 
		!	recursion is used. Furthermore, for m=l-1, Plmbar(l-2,m) is assumed to be zero.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		if (lmax_comp /= 0) then
			ff1(2,1) = 1.0d0
			ff2(2,1) = 0.0d0
		endif
				
		do l=2, lmax_comp, 1
			ff1(l+1,1) = dble(2*l-1) /dble(l)
			ff2(l+1,1) = dble(l-1) /dble(l)
			do m=1, l-2, 1
				ff1(l+1,m+1) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                  		ff2(l+1,m+1) = sqr(l-m-1) * sqr(l+m-1) / sqr(l+m) / sqr(l-m)
			enddo
			ff1(l+1,l)= dble(2*l-1) / sqr(l+m) / sqr(l-m)
                  	ff2(l+1,l) = 0.0d0
		enddo
	
		lmax_old = lmax_comp
	
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Determine Clms one l at a time by intergrating over latitude.
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	call dfftw_plan_dft_c2r_1d(plan, nlong, coef(1:nlong/2+1), grid(1:nlong), FFTW_MEASURE)

	i_eq = n/2 + 1	! Index correspondong to zero latitude

	do i=1, i_eq - 1, 1
	
		i_s = 2*i_eq -i
	
		theta = pi * dble(i-1)/dble(n)
		z = cos(theta)
		u = sqrt( (1.0d0-z) * (1.0d0+z) )
		
		lat = pi/2.0d0 - theta
		
		if (i==1) then		! Reference ellipsoid radius
			r_ex = a * (1.0d0 - f)	
		else
			r_ex = (1.0d0 + tan(lat)**2) / &
				(1.0d0  + tan(lat)**2 / (1.0d0 - f)**2 )
			r_ex = a * sqrt(r_ex)
		endif

		coefr(1:lmax+1) = dcmplx(0.0d0,0.0d0)
		coefr0 = 0.0d0
		coefrs(1:lmax+1) = dcmplx(0.0d0,0.0d0)
		coefrs0 = 0.0d0
		
		coeft(1:lmax+1) = dcmplx(0.0d0,0.0d0)
		coeft0 = 0.0d0
		coefts(1:lmax+1) = dcmplx(0.0d0,0.0d0)
		coefts0 = 0.0d0
		
		coefp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
		coefp0 = 0.0d0
		coefps(1:lmax+1) = dcmplx(0.0d0,0.0d0)
		coefps0 = 0.0d0
		
		if (calcu) then
			coefu(1:lmax+1) = dcmplx(0.0d0,0.0d0)
			coefu0 = 0.0d0
			coefus(1:lmax+1) = dcmplx(0.0d0,0.0d0)
			coefus0 = 0.0d0
		endif
		
		pm2 = 1.0d0

		! l = 0 terms are zero
				
		prefactor(1) = (r0/r_ex)**2
		do l=2, lmax_comp, 1
			prefactor(l) = prefactor(l-1)*r0/r_ex
		enddo
			
		pm1 =  ff1(2,1) * z	
		tempr = cilm(1,2,1) * pm1 * (-2) * prefactor(1)	! -2 = (l+1) prefactor
		coefr0 = coefr0 + tempr
		coefrs0 = coefrs0 - tempr 	! fsymsign = -1
			
		if (calcu) then
			tempr = cilm(1,2,1) * pm1 * prefactor(1)
			coefu0 = coefu0 + tempr
			coefus0 = coefus0 - tempr 	! fsymsign = -1
		endif
			
		dp = ff1(2,1)
		tempr = cilm(1,2,1) * dp * prefactor(1)
		coeft0 = coeft0 + tempr
		coefts0 = coefts0 + tempr 	! reverse fsymsign
				
		do l=2, lmax_comp, 1
			l1 = l+1
			p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
			tempr = cilm(1,l1,1) * p * (-l1) * prefactor(l)
			coefr0 = coefr0 + tempr
			coefrs0 = coefrs0 + tempr * fsymsign(l1,1)
			
			if (calcu) then
				tempr = cilm(1,l1,1) * p * prefactor(l)
				coefu0 = coefu0 + tempr
				coefus0 = coefus0 + tempr * fsymsign(l1,1)
			endif
			
			dp = l * (pm1 - z * p) / u**2

			tempr = cilm(1,l1,1) * dp * prefactor(l)
			coeft0 = coeft0 + tempr
			coefts0 = coefts0 - tempr * fsymsign(l1,1)	! reverse fsymsign
			
			pm2 = pm1
			pm1 = p
		enddo
				
		pmm = sqr(2) * scalef
				
		rescalem = 1.0d0/scalef
			
		do m = 1, lmax_comp-1, 1
			
			m1 = m+1
			rescalem = rescalem * u
					
			pmm = pmm * sqr(2*m+1) / sqr(2*m)
			pm2 = pmm / sqr(2*m+1)
			
			tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 * (-m-1) * prefactor(m)	! (m,m)
			coefr(m1) = coefr(m1) + tempc
			coefrs(m1) = coefrs(m1) + tempc
			! fsymsign = 1
			
			if (calcu) then
				tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 * prefactor(m)	! (m,m)
				coefu(m1) = coefu(m1) + tempc
				coefus(m1) = coefus(m1) + tempc
				! fsymsign = 1
			endif
			
			tempc = dcmplx(cilm(2,m1,m1), cilm(1,m1,m1)) * pm2 * prefactor(m) * m ! (m,m)
			coefp(m1) = coefp(m1) + tempc
			coefps(m1) = coefps(m1) + tempc
			! fsymsign = 1
			
			dp = -m * z * pm2 / u**2
			tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * dp * prefactor(m)	! (m,m)
			coeft(m1) = coeft(m1) + tempc
			coefts(m1) = coefts(m1) - tempc	! reverse fsymsign
										
	   		pm1 = z * ff1(m1+1,m1) * pm2		
	   		tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * pm1 * (-m-2) * prefactor(m+1)	! (m+1,m)
	   		coefr(m1) = coefr(m1) + tempc	
			coefrs(m1) = coefrs(m1) - tempc
			! fsymsign = -1
			
			if (calcu) then
				tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * pm1 * prefactor(m+1)	! (m+1,m)
	   			coefu(m1) = coefu(m1) + tempc	
				coefus(m1) = coefus(m1) - tempc
			endif
			
			tempc = dcmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1)) * pm1 * prefactor(m+1) * m	! (m+1,m)
	   		coefp(m1) = coefp(m1) + tempc	
			coefps(m1) = coefps(m1) - tempc
			! fsymsign = -1
			
			dp = (pm2 * sqr(2*m+1) - z * (m+1) * pm1) / u**2

			tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * dp * prefactor(m+1)	! (m+1,m)
	   		coeft(m1) = coeft(m1) + tempc	
			coefts(m1) = coefts(m1) + tempc	! reverse fsymsign
	   				
			do l = m+2, lmax_comp, 1
				l1 = l+1
				p = z * ff1(l1,m1) * pm1 - ff2(l1,m1) * pm2
				tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p * (-l1) * prefactor(l)
				coefr(m1) = coefr(m1) + tempc
				coefrs(m1) = coefrs(m1) + tempc * fsymsign(l1,m1)
				
				if (calcu) then
					tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p * prefactor(l)
					coefu(m1) = coefu(m1) + tempc
					coefus(m1) = coefus(m1) + tempc * fsymsign(l1,m1)
				endif
				
				tempc = dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) * p * prefactor(l) * m
				coefp(m1) = coefp(m1) + tempc
				coefps(m1) = coefps(m1) + tempc * fsymsign(l1,m1)
				
				dp = ( sqr(l+m) * sqr(l-m) * pm1 - l * z * p ) / u**2
      				tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * dp * prefactor(l)
				coeft(m1) = coeft(m1) + tempc
				coefts(m1) = coefts(m1) - tempc * fsymsign(l1,m1) ! reverse fsymsign
      		
				pm2 = pm1
       				pm1 = p
       				
			enddo
					
			coefr(m1) = coefr(m1) * rescalem
			coefrs(m1) = coefrs(m1) * rescalem
			
			if (calcu) then
				coefu(m1) = coefu(m1) * rescalem
				coefus(m1) = coefus(m1) * rescalem
			endif
			
			coeft(m1) = coeft(m1) * rescalem
			coefts(m1) = coefts(m1) * rescalem
			
			coefp(m1) = coefp(m1) * rescalem
			coefps(m1) = coefps(m1) * rescalem
					
		enddo			
								
		rescalem = rescalem * u
				
         	pmm = pmm / sqr(2*lmax_comp) * rescalem         			
        	tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), - cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
			* (-lmax_comp-1) * prefactor(lmax_comp)
        	coefr(lmax_comp+1) = coefr(lmax_comp+1) + tempc
		coefrs(lmax_comp+1) = coefrs(lmax_comp+1) + tempc
		! fsymsign = 1
		
		if (calcu) then
			tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), - cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
				* prefactor(lmax_comp)
        		coefu(lmax_comp+1) = coefu(lmax_comp+1) + tempc
			coefus(lmax_comp+1) = coefus(lmax_comp+1) + tempc
		endif
		
		tempc = dcmplx(cilm(2,lmax_comp+1,lmax_comp+1), cilm(1,lmax_comp+1,lmax_comp+1)) * pmm &
			 * prefactor(lmax_comp) * lmax_comp
        	coefp(lmax_comp+1) = coefp(lmax_comp+1) + tempc
		coefps(lmax_comp+1) = coefps(lmax_comp+1) + tempc
		! fsymsign = 1
		
		dp = -lmax_comp * z * pmm / u**2
		tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), - cilm(2,lmax_comp+1,lmax_comp+1)) * dp &
			* prefactor(lmax_comp)
        	coeft(lmax_comp+1) = coeft(lmax_comp+1) + tempc
		coefts(lmax_comp+1) = coefts(lmax_comp+1) - tempc	! reverse fsymsign
		
		
		coef(1) = dcmplx(coefr0,0.0d0)
		coef(2:lmax+1) = coefr(2:lmax+1)/2.0d0
		if (present(sampling)) then
			if (sampling == 2) then
				coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
			endif
		endif		
               	call dfftw_execute(plan)	! take fourier transform
               	rad_grid(i,1:nlong) = - grid(1:nlong) * r0/r_ex
               	
               	if (calcu) then
               		coef(1) = dcmplx(coefu0,0.0d0)
			coef(2:lmax+1) = coefu(2:lmax+1)/2.0d0
			if (present(sampling)) then
				if (sampling == 2) then
					coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
				endif
			endif		
               		call dfftw_execute(plan)	! take fourier transform
               		pot_grid(i,1:nlong) = grid(1:nlong) * r0
               	endif
               	
               	if (i==1) then
               		theta_grid(1,1:nlong) = 0.0d0	! This derivative is undefined at the pole
               		phi_grid(1,1:nlong) = 0.0d0	! This derivative is undefined at the pole
               	else
              	 	coef(1) = dcmplx(coeft0,0.0d0)
			coef(2:lmax+1) = coeft(2:lmax+1)/2.0d0
			if (present(sampling)) then
				if (sampling == 2) then
					coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
				endif
			endif		
              	 	call dfftw_execute(plan)	! take fourier transform
               		theta_grid(i,1:nlong) = sin(theta)*grid(1:nlong) * r0/r_ex 
               		
               		coef(1) = dcmplx(coefp0,0.0d0)
			coef(2:lmax+1) = coefp(2:lmax+1)/2.0d0
			if (present(sampling)) then
				if (sampling == 2) then
					coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
				endif
			endif		
              	 	call dfftw_execute(plan)	! take fourier transform
               		phi_grid(i,1:nlong) = - grid(1:nlong) * (r0/r_ex) / sin(theta)
               	endif
               	
               	if (i /= 1) then	! don't compute value for south pole.
               		coef(1) = dcmplx(coefrs0,0.0d0)
			coef(2:lmax+1) = coefrs(2:lmax+1)/2.0d0
			if (present(sampling)) then
				if (sampling == 2) then
					coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
				endif
			endif
               		call dfftw_execute(plan)	! take fourier transform
               		rad_grid(i_s,1:nlong) = - grid(1:nlong) * r0/r_ex
               		
               		if (calcu) then
               			coef(1) = dcmplx(coefus0,0.0d0)
				coef(2:lmax+1) = coefus(2:lmax+1)/2.0d0
				if (present(sampling)) then
					if (sampling == 2) then
						coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
					endif
				endif
               			call dfftw_execute(plan)	! take fourier transform
               			pot_grid(i_s,1:nlong) = grid(1:nlong) * r0
               		endif
               		
               		coef(1) = dcmplx(coefts0,0.0d0)
			coef(2:lmax+1) = coefts(2:lmax+1)/2.0d0
			if (present(sampling)) then
				if (sampling == 2) then
					coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
				endif
			endif
               		call dfftw_execute(plan)	! take fourier transform
               		theta_grid(i_s,1:nlong) = sin(theta)*grid(1:nlong) * r0/r_ex
               		
               		coef(1) = dcmplx(coefps0,0.0d0)
			coef(2:lmax+1) = coefps(2:lmax+1)/2.0d0
			if (present(sampling)) then
				if (sampling == 2) then
					coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
				endif
			endif
               		call dfftw_execute(plan)	! take fourier transform
               		phi_grid(i_s,1:nlong) = - grid(1:nlong) * (r0/r_ex) / sin(theta)
               		
               	endif
			
	enddo
	
	! Finally, do equator
	
	r_ex = a
	theta = pi/2.0d0
	z = 0.0d0
	u = 1.0d0

	coefr(1:lmax+1) = dcmplx(0.0d0,0.0d0)
	coefr0 = 0.0d0
	
	if (calcu) then
		coefu(1:lmax+1) = dcmplx(0.0d0,0.0d0)
		coefu0 = 0.0d0
	endif
		
	coeft(1:lmax+1) = dcmplx(0.0d0,0.0d0)
	coeft0 = 0.0d0
		
	coefp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
	coefp0 = 0.0d0
		
	pm2 = 1.0d0
			
	prefactor(1) = (r0/r_ex)**2
	do l=2, lmax_comp,1
		prefactor(l) = prefactor(l-1)*r0/r_ex
	enddo
			
	pm1 =  0.0d0	
		
	if (calcu) coefu0 = coefu0 + cilm(1,2,1) * pm1 * prefactor(1)
	
	dp = ff1(2,1)
	coeft0 = coeft0 + cilm(1,2,1) * dp * prefactor(1)
				
	do l=2, lmax_comp, 1
		l1 = l+1
		p = - ff2(l1,1) * pm2
		coefr0 = coefr0 + cilm(1,l1,1) * p * (-l1) * prefactor(l)
		
		if (calcu) coefu0 = coefu0 + cilm(1,l1,1) * p * prefactor(l)

		dp = l * (pm1) 
		coeft0 = coeft0 + cilm(1,l1,1) * dp * prefactor(l)
			
		pm2 = pm1
		pm1 = p
	enddo
				
	pmm = sqr(2) * scalef
				
	rescalem = 1.0d0/scalef
			
	do m = 1, lmax_comp-1, 1
			
		m1 = m+1
					
		pmm = pmm * sqr(2*m+1) / sqr(2*m)
		pm2 = pmm / sqr(2*m+1)
			
		coefr(m1) = coefr(m1) + dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 * (-m-1) * prefactor(m)
		
		if (calcu) coefu(m1) = coefu(m1) + dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 * prefactor(m)
	
		coefp(m1) = coefp(m1) + dcmplx(cilm(2,m1,m1), cilm(1,m1,m1)) * pm2 * prefactor(m) * m
										
	   	pm1 = 0.0d0	
	
		dp =  (pm2 * sqr(2*m+1))
	   	coeft(m1) = coeft(m1) + dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * dp * prefactor(m+1)
					
		do l = m+2, lmax_comp, 1
			l1 = l+1
			p =  - ff2(l1,m1) * pm2
			coefr(m1) = coefr(m1) + dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p * (-l1) * prefactor(l)
			
			if (calcu) coefu(m1) = coefu(m1) + dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p * prefactor(l)
		
			coefp(m1) = coefp(m1) + dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) * p * prefactor(l) * m
				
			dp = ( sqr(l+m) * sqr(l-m) * pm1)
			coeft(m1) = coeft(m1) + dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * dp * prefactor(l)
      	
			pm2 = pm1
       			pm1 = p
       				
		enddo
					
		coefr(m1) = coefr(m1) * rescalem
		
		if (calcu) coefu(m1) = coefu(m1) * rescalem
		
		coefp(m1) = coefp(m1) * rescalem
			
		coeft(m1) = coeft(m1) * rescalem
					
	enddo
				
        pmm = pmm / sqr(2*lmax_comp) * rescalem			
	coefr(lmax_comp+1) = coefr(lmax_comp+1) + dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), - cilm(2,lmax_comp+1,lmax_comp+1)) &
		* pmm * (-lmax_comp-1) * prefactor(lmax_comp)
		
	if (calcu) coefu(lmax_comp+1) = coefu(lmax_comp+1) + dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), - cilm(2,lmax_comp+1,lmax_comp+1)) &
		* pmm * prefactor(lmax_comp)
	
        coefp(lmax_comp+1) = coefp(lmax_comp+1) + dcmplx(cilm(2,lmax_comp+1,lmax_comp+1), cilm(1,lmax_comp+1,lmax_comp+1)) &
		* pmm * prefactor(lmax_comp) * lmax_comp
	
	dp = -lmax_comp * z * pmm / u**2
       	 coeft(lmax_comp+1) = coeft(lmax_comp+1) + dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), - cilm(2,lmax_comp+1,lmax_comp+1)) &
		* dp * prefactor(lmax_comp)

	coef(1) = dcmplx(coefr0,0.0d0)
	coef(2:lmax+1) = coefr(2:lmax+1)/2.0d0
	if (present(sampling)) then
		if (sampling == 2) then
			coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
		endif
	endif		
        call dfftw_execute(plan)	! take fourier transform
        rad_grid(i_eq,1:nlong) = - grid(1:nlong) * r0/r_ex
        
        if (calcu) then
        	coef(1) = dcmplx(coefu0,0.0d0)
		coef(2:lmax+1) = coefu(2:lmax+1)/2.0d0
		if (present(sampling)) then
			if (sampling == 2) then
				coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
			endif
		endif		
        	call dfftw_execute(plan)	! take fourier transform
        	pot_grid(i_eq,1:nlong) = grid(1:nlong) * r0
        endif
               	
        coef(1) = dcmplx(coeft0,0.0d0)
	coef(2:lmax+1) = coeft(2:lmax+1)/2.0d0
	if (present(sampling)) then
		if (sampling == 2) then
			coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
		endif
	endif		
       	call dfftw_execute(plan)	! take fourier transform
	theta_grid(i_eq,1:nlong) = sin(theta)*grid(1:nlong) * r0/r_ex 
               		
       	coef(1) = dcmplx(coefp0,0.0d0)
	coef(2:lmax+1) = coefp(2:lmax+1)/2.0d0
	if (present(sampling)) then
		if (sampling == 2) then
			coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0)	
		endif
	endif		
       	call dfftw_execute(plan)	! take fourier transform
       	phi_grid(i_eq,1:nlong) = - grid(1:nlong) * (r0/r_ex) / sin(theta)

	call dfftw_destroy_plan(plan)
	
	total_grid(1:n, 1:nlong) = sqrt(rad_grid(1:n,1:nlong)**2 + phi_grid(1:n,1:nlong)**2 + theta_grid(1:n,1:nlong)**2)
				
end subroutine MakeMagGridDH

