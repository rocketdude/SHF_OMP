subroutine MakeGridGLQC(gridglq, cilm, lmax, plx, zero, norm, csphase, lmax_calc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Given the complex spherical harmonic coefficients CILM, this subroutine
!	will evalate these coefficients on a grid with equal spacing in
!	longitude and with latitude points appropriate for Gauss-Legendre
!	quadrature integrations. Note that this is done using FFTs for each degree of
!	each latitude band. The grid spacing is determined by the spherical harmonic
!	bandwidth LMAX. Nevertheless, the coefficients can be evaluated up to
!	smaller spherical harmonic degree by specifying the optional parameter
!	LMAX_CALC.
!
!	The optional array PLX contains precomputed associated legendre functions
!	evaluated on the Gauss-Legendre quadrature nodes (obtained from PreCompute with CNORM=1)
!	and should not be precomputed when memory is an issue (i.e., lmax>360).
!	If PLX is not present, the Legendre functions are computed on the fly
!	using the scaling methodolgy presented in Holmes and Featherston (2002). 
!	When NORM=1,2 or 4, these are accurate to degree 2800. When NORM=3, the
!	routine is only stable to about degree 15!
!
!	Calling Parameters
!		IN
!			cilm		Input spherical harmonic coefficients with 
!					dimensions (2, LMAX+1, LMAX+1).
!			lmax		Maximum spherical harmonic degree used in the expansion.
!					This value determines the grid spacing of the output function.
!		OUT
!			gridglq		Gridded data of the spherical harmonic
!					coefficients CILM with dimensions (LMAX+1 , 2*LMAX+1). 
!					The first index (latitude) corresponds to the 
!					Gauss points, and the second index corresponds to 
!					360*(k-1)/nlong = 360*(k-1)/(2*LMAX +1). 
!		OPTIONAL (IN)
!			plx:		Input array of Associated Legendre Polnomials computed
!					at the Gauss points (determined from a call to
!					PreCompute). If this is not included, then the optional
!					array ZERO MUST be inlcuded. PLX must be computed using CNORM=1.
!			zero		Array of dimension lmax+1 that contains the latitudinal
!					gridpoints used in the Gauss-Legendre quadrature integration
!					scheme. Only needed if PLX is not included.
!			norm		Normalization to be used when calculating Legendre functions
!						(1) "geodesy" (default)
!						(2) Schmidt
!						(3) unnormalized
!						(4) orthonormalized
!			csphase	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!			lmax_calx	The maximum spherical harmonic degree to evaluate the coefficients
!					up to.
!
!	Notes:
!		1.	If lmax is greater than the the maximum spherical harmonic
!			degree of the input file, then this file will be ZERO PADDED!
!			(i.e., those degrees after lmax are assumed to be zero).
!		2.	Latitudes are geocentric latitude.
!
!	Dependencies:	FFTW3, NGLQSH, CSPHASE_DEFAULT
!
!	Written by Mark Wieczorek (May 2008).
!
!	November 21, 2011	Fixed problem where saved variables used in Plm recursion were not recalculated
!				if NORM changed from one call to the next (with the same value of N).
!
!	Copyright (c) 2005-2011, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use FFTW3
	use SHTOOLS, only: NGLQSH, CSPHASE_DEFAULT
	
	implicit none
	
	complex*16, intent(in) :: 	cilm(:,:,:)
	real*8, intent(in), optional :: plx(:,:), zero(:)
	complex*16, intent(out) ::	gridglq(:,:)
	integer, intent(in) :: 	lmax
	integer, intent(in), optional :: norm, csphase, lmax_calc
	integer :: 		l, m, i, nlat, nlong, l1, m1, lmax_comp, phase, i_s, astat(4), lnorm, &
				k, kstart
	real*8 :: 		pi, scalef, rescalem, u, p, pmm, pm1, pm2, z
	complex*16 ::		coef(2*lmax+1), coefs(2*lmax+1), grid(2*lmax+1)
	integer*8 ::		plan
	real*8, save, allocatable ::	f1(:), f2(:), sqr(:), symsign(:)
	integer, save ::	lmax_old=0, norm_old = 0


	if (size(cilm(:,1,1)) < 2) then
		print*, "Error --- MakeGridGLQC"
		print*, "CILM must be dimensioned as (2, *, *)."
		print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	elseif (size(gridglq(1,:)) < 2*lmax+1 .or. size(gridglq(:,1)) < lmax+1 ) then
		print*, "Error --- MakeGridGLQC"
		print*, "GRIDGLQ must be dimensioned as (LMAX+1, 2*LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(gridglq(:,1)), size(gridglq(1,:))
		stop
	endif

	if (present(plx)) then
		if (size(plx(:,1)) < lmax+1 .or. size(plx(1,:)) < (lmax+1)*(lmax+2)/2) then
			print*, "Error --- MakeGridGLQC"
			print*, "PLX must be dimensioned as (LMAX+1, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
			print*, "Input array is dimensioned as ", size(plx(:,1)), size(plx(1,:))
			stop
		endif
	elseif (present(zero)) then
		if (size(zero) < lmax+1) then
			print*, "Error --- MakeGridGLQC"
			print*, "ZERO must be dimensioned as (LMAX+1) where LMAX is ", lmax
			print*, "Input array is dimensioned ", size(zero)
			stop
		endif
	else
		print*, "Error --- MakeGridGLQC"
		print*, "Either PLX or ZERO must be specified."
		stop
	endif
	
	if (present(norm)) then
		if (norm > 4 .or. norm < 1) then
			print*, "Error - MakeGridGLQC"
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
     			print*, "Error --- MakeGridGLQC"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			stop
     		else
     			phase = csphase
     		endif
     	else
     		phase = CSPHASE_DEFAULT
     	endif
     	
	
	pi = acos(-1.0d0)
	
	nlong = 2*lmax + 1
	
	nlat = NGLQSH(lmax)	! nlat = lmax+1
	
	scalef = 1.0d-280
	
	if (present(lmax_calc)) then
		lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1, lmax_calc)
	else
		lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1)
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate recursion constants used in computing the Legendre polynomials
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if ( (lmax_comp /= lmax_old .or. lnorm /= norm_old) .and. .not. present(plx)) then
		
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		if (allocated(symsign)) deallocate(symsign)
		
		allocate(sqr(2*lmax_comp+1), stat=astat(1))
		allocate(f1((lmax_comp+1)*(lmax_comp+2)/2), stat=astat(2))
		allocate(f2((lmax_comp+1)*(lmax_comp+2)/2), stat=astat(3))
		allocate(symsign((lmax_comp+1)*(lmax_comp+2)/2), stat=astat(4))
		
		if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 .or. astat(4) /= 0) then
			print*, "MakeGridGLQ --- Error"
			print*, "Problem allocating arrays SQR, F1, F2, or SYMSIGN", astat(1), astat(2), astat(3), astat(4)
			stop
		endif
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Calculate signs used for symmetry of Legendre functions about equator
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		k = 0
		do l = 0, lmax_comp, 1
			do m = 0, l, 1
				k = k + 1
				if (mod(l-m,2) == 0) then
					symsign(k) = 1.0d0
				else
					symsign(k) = -1.0d0
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
	
		k = 1
		
		select case(lnorm)
		
			case(1,4)
				if (lmax_comp /= 0) then
					k = k + 1
					f1(k) = sqr(3)
					f2(k) = 0.0d0
					k = k + 1
				endif
				
				do l=2, lmax_comp, 1
					k = k + 1
					f1(k) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
					f2(k) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)
					do m=1, l-2, 1
						k = k+1
						f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                				f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                  				 	/ sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
					enddo
					k = k+1
					f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                			f2(k) = 0.0d0
					k = k + 1
				enddo
			
			case(2)
				
				if (lmax_comp /= 0) then
					k = k + 1
					f1(k) = 1.0d0
					f2(k) = 0.0d0
					k = k + 1
				endif
				
				do l=2, lmax_comp, 1
					k = k + 1
					f1(k) = dble(2*l-1) /dble(l)
					f2(k) = dble(l-1) /dble(l)
					do m=1, l-2, 1
						k = k+1
						f1(k) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                  				f2(k) = sqr(l-m-1) * sqr(l+m-1) / sqr(l+m) / sqr(l-m)
					enddo
					k = k+1
					f1(k) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                  			f2(k) = 0.0d0
					k = k + 1
				enddo
			
			case(3)		
			
				do l=1, lmax_comp, 1
					k = k + 1
					f1(k) = dble(2*l-1) /dble(l)
					f2(k) = dble(l-1) /dble(l)
					do m=1, l-1, 1
						k = k+1
						f1(k) = dble(2*l-1) / dble(l-m)
                  				f2(k) = dble(l+m-1) / dble(l-m)
					enddo
					k = k + 1
				enddo

		end select
	
		lmax_old = lmax_comp
		norm_old = lnorm
	
	endif	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Do special case of lmax_comp = 0
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (lmax_comp == 0) then
	
		select case(lnorm)
			case(1,2,3);	pm2 = 1.0d0
			case(4);	pm2 = 1.0d0 / sqrt(4.0d0*pi)
		end select
		
		gridglq(1:nlat, 1:nlong) = cilm(1,1,1) * pm2
	
		return
	
	endif

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Determine Cilms, one l at a time I by integrating over all 
	!	latitudes using Gauss-Legendre Quadrature. When PLX is not 
	!	present, the Legendre functions are computed on the fly 
	!	during the summations over l and m. These are scaled using
	!	the methodology of Holmesand Featherstone (2002), with the
	!	exception of the m=0 terms that do not need to be scaled 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	call dfftw_plan_dft_1d(plan, nlong, coef(1:nlong), grid(1:nlong), FFTW_BACKWARD, FFTW_MEASURE)	! create generic plan
	
	if (present(plx)) then
	
		do i = 1, nlat
			coef = dcmplx(0.0d0,0.0d0)
			! This summation order is intended to add the smallest terms first
			do l = lmax_comp, 0, -1
				l1 = l+1
				k =  l1*l/2 + 1	! m=0
				coef(1) = coef(1) + cilm(1,l1,1) * plx(i,k)
				do m = 1, l, 1
					m1 = m+1
					k = l1*l/2 + m1
					coef(m1) = coef(m1) + cilm(1,l1,m1) * plx(i,k)
					coef(nlong-(m-1)) = coef(nlong-(m-1)) + &
						cilm(2,l1,m1) * plx(i,k) * dble((-1)**mod(m,2))
				enddo
			enddo
			
               		call dfftw_execute(plan)	! take fourier transform
                
			gridglq(i,1:nlong) = grid(1:nlong)
		enddo
		
	else

		do i=1, (nlat+1)/2
			
			coef = dcmplx(0.0d0,0.0d0)
		
			if ( i==(nlat+1)/2 .and. mod(nlat,2) /=0) then		! This latitude is the equator; z=0, u=1
				
				u = 1.0d0
				
				select case(lnorm)
					case(1,2,3);	pm2 = 1.0d0
					case(4);	pm2 = 1.0d0 / sqrt(4.0d0*pi)
				end select

				coef(1) = coef(1) + cilm(1,1,1) * pm2
				
				k = 2
				
				do l=2, lmax_comp, 2
					l1 = l+1
					k = k + l
					p = - f2(k) * pm2
					pm2 = p
					coef(1) = coef(1) + cilm(1,l1,1) * p
					k = k + l + 1
				enddo
				
				select case(lnorm)
					case(1,2,3);	pmm = scalef
					case(4);	pmm = scalef / sqrt(4.0d0*pi)
				end select
								
				rescalem = 1.0d0/scalef
				kstart = 1
			
				do m = 1, lmax_comp-1, 1
				
					m1 = m + 1
					kstart = kstart + m + 1
					
					select case(lnorm)
						case(1,4)
							pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
							pm2 = pmm
						case(2)
							pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
							pm2 = pmm / sqr(2*m+1)
						case(3)
							pmm = phase * pmm * dble(2*m-1)
							pm2 = pmm
					end select
					
					coef(m1) = coef(m1) + cilm(1,m1,m1) * pm2
					coef(nlong-(m-1)) = coef(nlong-(m-1)) + cilm(2,m1,m1) * pm2 
										
					k = kstart+m+1
	   				
					do l = m+2, lmax_comp, 2
						l1 = l+1
						k = k + l
						p = - f2(k) * pm2
						coef(m1) = coef(m1) + cilm(1,l1,m1) * p
						coef(nlong-(m-1)) = coef(nlong-(m-1)) + cilm(2,l1,m1) * p 
						pm2 = p
						k = k + l + 1
					enddo
					
					coef(m1) = coef(m1) * rescalem
					coef(nlong-(m-1)) = coef(nlong-(m-1)) * rescalem * dble((-1)**mod(m,2))
										
				enddo			
				
           			select case(lnorm)
            				case(1,4);	pmm = phase * pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem
            				case(2);	pmm = phase * pmm / sqr(2*lmax_comp) * rescalem
            				case(3);	pmm = phase * pmm * dble(2*lmax_comp-1) * rescalem
        			end select
          			
        			coef(lmax_comp+1) = coef(lmax_comp+1) + cilm(1,lmax_comp+1,lmax_comp+1) * pmm
        			coef(nlong-(lmax_comp-1)) = coef(nlong-(lmax_comp-1)) + &
        				cilm(2,lmax_comp+1,lmax_comp+1) * pmm  * dble((-1)**mod(lmax_comp,2))
		
               			call dfftw_execute(plan)	! take fourier transform
                
				gridglq(i,1:nlong) = grid(1:nlong)
			
			else
			
				z = zero(i)
				u = sqrt( (1.0d0-z) * (1.0d0+z) )
				
				i_s = nlat+1-i

				coefs = dcmplx(0.0d0,0.0d0)
				
				select case(lnorm)
					case(1,2,3);	pm2 = 1.0d0
					case(4);	pm2 = 1.0d0 / sqrt(4.0d0*pi)
				end select

				coef(1) = coef(1) + cilm(1,1,1) * pm2
				coefs(1) = coefs(1) + cilm(1,1,1) * pm2 	! symsign is always 1 for l=m=0
				
				k = 2
				pm1 =  f1(k) * z * pm2
				coef(1) = coef(1) + cilm(1,2,1) * pm1
				coefs(1) = coefs(1) + cilm(1,2,1) * pm1 * symsign(k)
				
				do l=2, lmax_comp, 1
					l1 = l+1
					k = k+l
					p = f1(k) * z * pm1 - f2(k) * pm2
					coef(1) = coef(1) + cilm(1,l1,1) * p
					coefs(1) = coefs(1) + cilm(1,l1,1) * p * symsign(k)
					pm2 = pm1
					pm1 = p
				enddo
				
				select case(lnorm)
					case(1,2,3);	pmm = scalef
					case(4);	pmm = scalef / sqrt(4.0d0*pi)
				end select
				
				rescalem = 1.0d0/scalef
				kstart = 1
			
				do m = 1, lmax_comp-1, 1
				
					m1 = m+1
					rescalem = rescalem * u
					kstart = kstart+m+1
					
					select case(lnorm)
						case(1,4)
							pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
							pm2 = pmm
						case(2)
							pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
							pm2 = pmm / sqr(2*m+1)
						case(3)
							pmm = phase * pmm * dble(2*m-1)
							pm2 = pmm
					end select
					
					coef(m1) = coef(m1) + cilm(1,m1,m1) * pm2
					coef(nlong-(m-1)) = coef(nlong-(m-1)) + cilm(2,m1,m1) * pm2
					coefs(m1) = coefs(m1) + cilm(1,m1,m1) * pm2
					coefs(nlong-(m-1)) = coefs(nlong-(m-1)) + cilm(2,m1,m1) * pm2
					! symsign(kstart) = 1
										
					k = kstart+m+1
	   				pm1 = z * f1(k) * pm2
	   				
	   				coef(m1) = coef(m1) + cilm(1,m1+1,m1) * pm1
	   				coef(nlong-(m-1)) = coef(nlong-(m-1)) + cilm(2,m1+1,m1) * pm1	
					coefs(m1) = coefs(m1) - cilm(1,m1+1,m1) * pm1
					coefs(nlong-(m-1)) = coefs(nlong-(m-1)) - cilm(2,m1+1,m1) * pm1
					! symsign = -1
	   				
					do l = m+2, lmax_comp, 1
						l1 = l+1
						k = k + l
						p = z * f1(k) * pm1 - f2(k) * pm2
						pm2 = pm1
                  				pm1 = p
						coef(m1) = coef(m1) + cilm(1,l1,m1) * p
						coef(nlong-(m-1)) = coef(nlong-(m-1)) + cilm(2,l1,m1) * p
						coefs(m1) = coefs(m1) + cilm(1,l1,m1) * p * symsign(k)
						coefs(nlong-(m-1)) = coefs(nlong-(m-1)) + cilm(2,l1,m1) * p * symsign(k)
					enddo
					
					coef(m1) = coef(m1) * rescalem
					coef(nlong-(m-1)) = coef(nlong-(m-1)) * rescalem * dble((-1)**mod(m,2))
					coefs(m1) = coefs(m1) * rescalem
					coefs(nlong-(m-1)) = coefs(nlong-(m-1)) * rescalem * dble((-1)**mod(m,2))
					
				enddo			
								
				rescalem = rescalem * u
				
          			select case(lnorm)
            				case(1,4);	pmm = phase * pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem
            				case(2);	pmm = phase * pmm / sqr(2*lmax_comp) * rescalem
            				case(3);	pmm = phase * pmm * dble(2*lmax_comp-1) * rescalem
        			end select
         			
        			coef(lmax_comp+1) = coef(lmax_comp+1) + cilm(1,lmax_comp+1,lmax_comp+1) * pmm
        			coef(nlong-(lmax_comp-1)) = coef(nlong-(lmax_comp-1)) + cilm(2,lmax_comp+1,lmax_comp+1) &
        				* pmm * dble((-1)**mod(lmax_comp,2))
				coefs(lmax_comp+1) = coefs(lmax_comp+1) + cilm(1,lmax_comp+1,lmax_comp+1) * pmm
				coefs(nlong-(lmax_comp-1)) = coefs(nlong-(lmax_comp-1)) + cilm(2,lmax_comp+1,lmax_comp+1) &
					* pmm * dble((-1)**mod(lmax_comp,2))
				! symsign = 1
               			
               			call dfftw_execute(plan)	! take fourier transform
               			gridglq(i,1:nlong) = grid(1:nlong)
               			
				coef(1:nlong) = coefs(1:nlong)
               			
               			call dfftw_execute(plan)	! take fourier transform
               			gridglq(i_s,1:nlong) = grid(1:nlong)
 

			endif
			
		enddo
		
	endif
	
	call dfftw_destroy_plan(plan)
				
end subroutine MakeGridGLQC

