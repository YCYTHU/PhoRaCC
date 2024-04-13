module module_fitDecayData
    use module_print
    logical :: useCorrected, doCompute
    contains

    recursive subroutine getParameters(kp, kd, phiPF, phiDF, doCompute, mainOpt)
        implicit none
        character (len=2) :: parameterType
        character :: mainOpt
        logical :: useCorrected, doCompute
        real*8 :: kp, kd, phiPF, phiDF, taup, taud, Ap, Ad, phiPLQY, denominator
        doCompute =.true.

        write(*,*) "                 ------------ Select parameter type ------------"
        write(*,*) "-1. Display formulas used for calculations"
        write(*,*) "0. Return"
        write(*,*) "1. Compute Rate Constants with tau_p, tau_d, Ap, Ad and Phi_PLQY"
        write(*,*) "2. Compute Rate Constants with kp, kd, Ap, Ad and Phi_PLQY"
        write(*,*) "3. Compute Rate Constants with tau_p, tau_d, Phi_PF and Phi_DF"
        write(*,*) "4. Compute Rate Constants with kp, kd, Phi_PF and Phi_DF"
        read(*,*) parameterType
        select case (parameterType)
            case ('-1')
                call printFormula(mainOpt)
                call getParameters(kp, kd, phiPF, phiDF, doCompute, mainOpt)
            case ('0')
                doCompute = .false.
                return
            case ('1')
                call getCompomentMethod(useCorrected)
                write(*,*) "Input tau_p (ns), tau_d (us), Ap, Ad and Phi_PLQY"
                read(*,*) taup, taud, Ap, Ad, phiPLQY
                kp = 1.D9 / taup
                kd = 1.D6 / taud
                denominator = Ap * kd + Ad * kp
                if (useCorrected) then
                    phiPF = (Ap + Ad) * kd * phiPLQY / denominator
                    phiDF = Ad * (kp - kd) * phiPLQY / denominator
                else
                    phiPF = Ap * kd * phiPLQY / denominator
                    phiDF = Ad * kp * phiPLQY / denominator
                end if
            case ('2')
                call getCompomentMethod(useCorrected)
                write(*,*) "Input kp (s-1), kd (s-1), Ap, Ad and Phi_PLQY"
                read(*,*) kp, kd, Ap, Ad, phiPLQY
                denominator = Ap * kd + Ad * kp
                if (useCorrected) then
                    phiPF = (Ap + Ad) * kd * phiPLQY / denominator
                    phiDF = Ad * (kp - kd) * phiPLQY / denominator
                else
                    phiPF = Ap * kd * phiPLQY / denominator
                    phiDF = Ad * kp * phiPLQY / denominator
                end if
            case ('3')
                write(*,*) "Input tau_p (ns), tau_d (us), Phi_PF and Phi_DF"
                read(*,*) taup, taud, phiPF, phiDF
                kp = 1.D9 / taup
                kd = 1.D6 / taud
            case ('4')
                write(*,*) "Input kp (s-1), kd (s-1), Phi_PF and Phi_DF"
                read(*,*) kp, kd, phiPF, phiDF
            case default
                stop
        end select
    end subroutine

    subroutine getCompomentMethod(useCorrected)
        implicit none
        logical :: useCorrected
        integer :: compomentMethod
        call printFitMethod()
        read(*,*) compomentMethod
        if (compomentMethod == 1) then
            useCorrected = .false.
        elseif (compomentMethod == 2) then
            useCorrected = .true.
        else
            stop
        end if
    end subroutine
end module
