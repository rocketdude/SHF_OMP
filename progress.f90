   SUBROUTINE progress(it, reinit, Maxit, t, Uave)
 
     implicit none
     integer(kind=4) :: it, reinit, Maxit, j, k
     real(kind=8) :: t, Uave
     character(len=98) :: bar="???% |                                                  | t=             U=                       "
     character(len=98) :: bar2="???% |                                                  | t=             U=             ==REINIT=="

                              
     j = it*50/Maxit

     IF( mod(it,reinit) .eq. 0 ) THEN

        ! updates the fraction of calculation done
        write(unit=bar2(1:3),fmt='(i3)') 2*j
        do k = 1, j
           bar2(6+k:6+k)="*"
        enddo

        write(unit=bar2(61:72), fmt='(f12.9)') t
        write(unit=bar2(76:87), fmt='(f12.9)') Uave

        write(*,fmt='(a1,a98)',advance='no') char(13), bar2
     ELSE

        ! updates the fraction of calculation done
        write(unit=bar(1:3),fmt='(i3)') 2*j
        do k = 1, j
           bar(6+k:6+k)="*"
        enddo

        write(unit=bar(61:72), fmt='(f12.9)') t
        write(unit=bar(76:87), fmt='(f12.9)') Uave

        ! print the progress bar.
        write(*,fmt='(a1,a98)',advance='no') char(13), bar
     END IF

     RETURN
   END SUBROUTINE progress
   
