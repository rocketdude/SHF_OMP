!----------------------------------------------------------!
!     MODULE Dynamic Metric Array                          !
!----------------------------------------------------------!

    MODULE DynMetricArray

    !This module contains the dynamically allocated metric array
    !used to read metric data from hdf5 files

    INTEGER*4                               :: metricdims(3)
    REAL*8, DIMENSION(:,:,:), ALLOCATABLE   :: hdfmetric
    
    !NOTE       !hdfmetric will behave like a global array to all subroutines
                !that uses DynMetricArray. Remember to deallocate after using.
                !Otherwise, might get segfault.


    END MODULE

