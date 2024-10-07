! |-------<THOMAS RELATIVISTIC ELECTRONIC STRUCTURE CALCULATIONS>------|
! |------------------------------<DIRAC4¦°>-----------------------------|
! |----------------------------<DEVELOPMENT>---------------------------|
! |----------------------<LIUKUNYU000213@163.COM>----------------------|

program main
    use Indicator
    implicit none
    
    integer :: mainloop_i, mainloop_j, mainloop_k
    character(len=2) :: charai, charaj
    character(len=66) :: command
    maxloop = 100
    do mainloop_i = 7, 10
        nudge = real(mainloop_i)/10.0
        do mainloop_j = 4, 30
            subsp = mainloop_j
            nodiis = subsp + 10
            write(charai,'(I2.2)') mainloop_i
            write(charaj,'(I2.2)') mainloop_j
            command = "copy E:\TRESC\H2\H2.tot E:\TRESC\H2\S-0.5rouS-0.5DIIStest\"//charai//charaj//".tot"
            call elec_stc_calc()
            close(60)
            call EXECUTE_COMMAND_LINE(command)
            CALL SLEEPQQ(500)
        end do
        if (mainloop_i == 2) then
            do mainloop_j = 29, 30
                subsp = mainloop_j
                nodiis = subsp + 10
                write(charai,'(I2.2)') mainloop_i
                write(charaj,'(I2.2)') mainloop_j
                command = "copy E:\TRESC\H2\H2.tot E:\TRESC\H2\S-0.5rouS-0.5DIIStest\"//charai//charaj//".tot"
                call elec_stc_calc()
                close(60)
                call EXECUTE_COMMAND_LINE(command)
                CALL SLEEPQQ(500)
            end do
        end if
    end do
    pause
    stop
end program