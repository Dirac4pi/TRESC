! |-------<THOMAS RELATIVISTIC ELECTRONIC STRUCTURE CALCULATIONS>------|
! |------------------------------<DIRAC4PI>----------------------------|

program main
    use Indicator
    implicit none
    integer :: mainloop_i, mainloop_j, mainloop_k
    call sconf_calc(keep=0, kill=.true.)
    !do mainloop_i = 1, 10
    !    do batch calculation, plot with TRESCplot.py
    !end do
end program