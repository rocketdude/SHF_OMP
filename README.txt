Event Horizon Finder using Pseudospectral methods.
(for previous versions with finite-differencing, see LVLST)

03/12/2012
-	SHF: an event horizon finder that uses pseudospectral method.
	Data in SHF are mostly vectorized to reduce computational time.
	
03/25/2012
-	Uses cubic spline interpolation subroutines from 
	http://orion.math.iastate.edu/burkardt/f_src/spline/spline.html
	written by John Burkardt

04/07/2012
-	Parallelized the code with OpenMP

04/30/2012
-	Uses LAPACK calls to calculate the Moore-Penrose pseudoinverse.
-	Uses the pseudoinverse to calculate evolve the data.