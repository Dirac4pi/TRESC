! |-------<THOMAS RELATIVISTIC ELECTRONIC STRUCTURE CALCULATIONS>------|
! |------------------------------<DIRAC4PI>----------------------------|
! |----------------------------<DEVELOPMENT>---------------------------|
! |----------------------<LIUKUNYU000213@163.COM>----------------------|

program main
    use Indicator
    implicit none
    
    !integer :: mainloop_i, mainloop_j, mainloop_k
    !character(len=2) :: charai, charaj
    !character(len=67) :: command
    call elec_stc_calc()
    ! batch single-point energy calculations that can be used in conjunction with the TRESCplot.py script.
    !maxiter = 100
    !do mainloop_i = 1, 10
    !    damp = real(mainloop_i)/10.0
    !    do mainloop_j = 4, 30
    !        subsp = mainloop_j
    !        nodiis = subsp + 10
    !        write(charai,'(I2.2)') mainloop_i
    !        write(charaj,'(I2.2)') mainloop_j
    !        command = "copy E:\TRESC\H2\H2.tot E:\TRESC\H2\S-0.5ccS-0.5DIIStest-2\"//charai//charaj//".tot"
    !        call elec_stc_calc()
    !        close(60)
    !        call EXECUTE_COMMAND_LINE(command)
    !        CALL SLEEPQQ(500)
    !    end do
    !end do
    pause
    stop
end program