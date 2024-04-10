program main
  implicit none
  character :: opt = '0'
  real*8 :: kp, kd, phiPF, phiDF
  real*8 :: krs, knrs, kisc, krt, knrt, krisc
  logical :: useCorrected, doCompute
  write(*,*) "PhoRaCC: Compute Photophysical Rate Constants based on various assumptions"
  write(*,*) "Programmed by Yaocy"
  write(*,*) "Release date: 2024-Apr-10"
  write(*,*)
  doCompute =.true.
  do while (.true.)
    write(*,*) "                   ************ Three-State TADF ************"
    !(DOI:10.1021/acs.jpca.1c04056)
    write(*,*) "1. Steady State Assumption & krS + knrS + kISC >> kRISC & Phi_nrS = 0 & krT = 0"
    write(*,*) "2. Steady State Assumption & krS + knrS + kISC >> kRISC & Phi_nrT = 0 & krT = 0"
    write(*,*) "3. krT = 0"
    write(*,*) "4. Phi_DF / Phi_PF > 4 & krT = knrT = 0"
    write(*,*) "5. krT = knrT = 0"
    read(*,*) opt
    select case (opt)
      case ('1')
        call getParameters(kp, kd, phiPF, phiDF, doCompute)
        if (.not. doCompute) then
          cycle
        end if
        call SSA_NOnrs_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printRate(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        write(*,*) "Tsuchiya, Y., Diesing, S., ... & Adachi, C. (2021). The Journal of Physical Chemistry A, 125(36), 8074-8089."
        write(*,*) "(DOI: 10.1021/acs.jpca.1c04056)"
        write(*,*)
        write(*,*)
      case ('2')
        call getParameters(kp, kd, phiPF, phiDF, doCompute)
        if (.not. doCompute) then
          cycle
        end if
        call SSA_NOnrt_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printRate(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        write(*,*) "Tsuchiya, Y., Diesing, S., ... & Adachi, C. (2021). The Journal of Physical Chemistry A, 125(36), 8074-8089."
        write(*,*) "(DOI: 10.1021/acs.jpca.1c04056)"
        write(*,*)
        write(*,*)
      case ('3')
        call getParameters(kp, kd, phiPF, phiDF, doCompute)
        if (.not. doCompute) then
          cycle
        end if
        call NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printRate(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        write(*,*) "Tsuchiya, Y., Diesing, S., ... & Adachi, C. (2021). The Journal of Physical Chemistry A, 125(36), 8074-8089."
        write(*,*) "(DOI: 10.1021/acs.jpca.1c04056)"
        write(*,*)
        write(*,*)
      case ('4')
        call getParameters(kp, kd, phiPF, phiDF, doCompute)
        if (.not. doCompute) then
          cycle
        end if
        call Monkman_LargeDF(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printRate(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        write(*,*) "Dias, F. B., Penfold, T. J., & Monkman, A. P. (2017). Methods and applications in fluorescence, 5(1), 012001."
        write(*,*) "(DOI: 10.1088/2050-6120/aa537e)"
        write(*,*)
        write(*,*)
      case ('5')
        call getParameters(kp, kd, phiPF, phiDF, doCompute)
        if (.not. doCompute) then
          cycle
        end if
        call Kaji_NOnrt_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printRate(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        write(*,*) "Wada, Y., Nakagawa, H., ... & Kaji, H. (2020). Nature Photonics, 14(10), 643-649."
        write(*,*) "(DOI: 10.1038/s41566-020-0667-0)"
        write(*,*)
        write(*,*)
      case ('q')
        exit
      case default
        cycle
    end select
  enddo

  contains

  subroutine printRate(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
    implicit none
    real*8 :: kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc
    write(*,*)
    write(*,*) "============= ============= ============= ============="
    write(*,*) "     kp            kd          Phi_PF        Phi_DF"
    write(*,*) "------------- ------------- ------------- -------------"
    write(*,"(3X,ES9.3,5X,ES9.3,5X,ES9.3,5X,ES9.3)") kp, kd, phiPF, phiDF
    write(*,*) "============= ============= ============= ============="
    write(*,*) "     krS          knrS          kISC           kS"
    write(*,*) "------------- ------------- ------------- -------------"
    write(*,"(3X,ES9.3,5X,ES9.3,5X,ES9.3,5X,ES9.3)") krs, knrs, kisc, krs + knrs + kisc
    write(*,*) "============= ============= ============= ============="
    write(*,*) "     krT          knrT          kRISC          kT"
    write(*,*) "------------- ------------- ------------- -------------"
    write(*,"(3X,ES9.3,5X,ES9.3,5X,ES9.3,5X,ES9.3)") krt, knrt, krisc, krt + knrt + krisc
    write(*,*) "============= ============= ============= ============="
    write(*,*)
  end subroutine

  subroutine getParameters(kp, kd, phiPF, phiDF, doCompute)
    implicit none
    character :: parameterType
    logical :: useCorrected, doCompute
    real*8 :: kp, kd, phiPF, phiDF, taup, taud, Ap, Ad, phiPLQY, denominator

    write(*,*) "                 ------------ Select parameter type ------------"
    write(*,*) "0. Return"
    write(*,*) "1. Compute Rate Constants with tau_p, tau_d, Ap, Ad and Phi_PLQY"
    write(*,*) "2. Compute Rate Constants with kp, kd, Ap, Ad and Phi_PLQY"
    write(*,*) "3. Compute Rate Constants with tau_p, tau_d, Phi_PF and Phi_DF"
    write(*,*) "4. Compute Rate Constants with kp, kd, Phi_PF and Phi_DF"
    read(*,*) parameterType
    select case (parameterType)
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

    write(*,*) "                    ------------ Fitting Method ------------"
    write(*,*)
    write(*,*) "General Method :"
    write(*,*)
    write(*,*) "PF : DF = Ap*kd : Ad*kp"
    write(*,*)
    write(*,*) "    +----------------------------------------------------------------------+   "
    write(*,*) "    |                                                                      |   "
    write(*,*) "    |                                                     PL decay $$$$$$$ |   "
    write(*,*) "    |$                                                          PF ******* |   "
    write(*,*) "    |$                                                          DF ####### |   "
    write(*,*) "    | $$                                                                   |   "
    write(*,*) "    |  *$                                                                  |   "
    write(*,*) "    |   *$$                                                                |   "
    write(*,*) "    |####* $$$                                                             |   "
    write(*,*) "    |   ######$$$$$                                                        |   "
    write(*,*) "    |      *    ##$$$$$$$$                                                 |   "
    write(*,*) "    |       *           ##$$$$$$$$                                         |   "
    write(*,*) "    |        *                   #$$$$$$$$                                 |   "
    write(*,*) "    |         *                          #$$$$$$$$$                        |   "
    write(*,*) "    |          *                                ###$$$$$$$$$               |   "
    write(*,*) "    |           *                                       ###$$$$$$$$$       |   "
    write(*,*) "    |            *                                               ###$$$$$$$|   "
    write(*,*) "    |             *                                                      ##|   "
    write(*,*) "    |              *                                                       |   "
    write(*,*) "    |               *                                                      |   "
    write(*,*) "    +----------------------------------------------------------------------+   "
    write(*,*)
    write(*,*) "Corrected Method (DOI:10.1021/acs.jpca.1c04056):"
    write(*,*)
    write(*,*) "PF : DF = (Ap+Ad)*kd : Ad*(kp-kd)"
    write(*,*)
    write(*,*) "    +----------------------------------------------------------------------+   "
    write(*,*) "    |                                                                      |   "
    write(*,*) "    |                                                     PL decay $$$$$$$ |   "
    write(*,*) "    |$                                                          PF ******* |   "
    write(*,*) "    |$                                                          DF ####### |   "
    write(*,*) "    | $$                                                                   |   "
    write(*,*) "    |  *$                                                                  |   "
    write(*,*) "    |   *$$                                                                |   "
    write(*,*) "    |    * $$$                                                             |   "
    write(*,*) "    |     *   $$$$$                                                        |   "
    write(*,*) "    |      *  ####$$$$$$$$                                                 |   "
    write(*,*) "    |      ###          ##$$$$$$$$                                         |   "
    write(*,*) "    |    ##  *                   #$$$$$$$$                                 |   "
    write(*,*) "    |   #     *                          #$$$$$$$$$                        |   "
    write(*,*) "    |   #      *                                ###$$$$$$$$$               |   "
    write(*,*) "    |  #        *                                       ###$$$$$$$$$       |   "
    write(*,*) "    |  #         *                                               ###$$$$$$$|   "
    write(*,*) "    | #           *                                                      ##|   "
    write(*,*) "    | #            *                                                       |   "
    write(*,*) "    |#              *                                                      |   "
    write(*,*) "    +----------------------------------------------------------------------+   "
    write(*,*)
    write(*,*) "Choose the Method for fitting PF & DF components. Input 1 for General Method, 2 for Corrected Method."
    write(*,*) "tips: You can pre-specify the method in settings.ini to avoid displaying the plots"
    read(*,*) compomentMethod
    if (compomentMethod == 1) then
        useCorrected = .false.
    elseif (compomentMethod == 2) then
        useCorrected = .true.
    else
        stop
    end if
  end subroutine

  subroutine SSA_NOnrs_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
    implicit none
    real*8 :: kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc
    knrs = 0.D0
    krt = 0.D0
    krs = kp * phiPF
    kisc = kp * (1.D0 - phiPF)
    krisc = (kd * phiDF) / (phiPF * (1.D0 - phiPF))
    knrt = kd - krisc * phiPF
  end subroutine

  subroutine SSA_NOnrt_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
    implicit none
    real*8 :: kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc, phiISC
    knrt = 0.D0
    krt = 0.D0
    krs = kp * phiPF
    phiISC = phiDF / (phiPF + phiDF)
    kisc = kp * phiISC
    knrs = kp * (1.D0 - phiPF - phiISC)
    krisc = (kd * phiDF) / (phiPF * phiISC)
  end subroutine

  subroutine NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
    implicit none
    real*8 :: kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc, phiPLQY
    phiPLQY = phiPF + phiDF
    krt = 0.D0
    knrt = 0.D0
    krs = kp * phiPF
    knrs = kp * phiPF * (1.D0 - phiPLQY) / phiPLQY
    kisc = kp * phiDF / phiPLQY - kd * phiDF / phiPF
    krisc = kd * phiPLQY / phiPF
  end subroutine

  subroutine Monkman_LargeDF(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
    implicit none
    real*8 :: kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc, phiPLQY
    krt = 0
    knrt = 0
    krs = kp * phiPF
    kisc = kp * phiDF / (phiPF + phiDF)
    knrs = kp - krs - kisc
    krisc = kd * (phiPF + phiDF) / phiDF
  end subroutine

  subroutine Kaji_NOnrt_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
    implicit none
    real*8 :: kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc, phiPLQY
    krt = 0
    knrt = 0
    krisc = (kp + kd) / 2 - SQRT(((kp + kd) / 2) ** 2 - kp * kd * (1 + phiDF / phiPF))
    kisc = kp * kd * phiDF / (krisc * phiPF)
    krs = kp * kd  * (phiPF + phiDF) / krisc
    knrs = kp * kd  * (1 - phiPF - phiDF) / krisc
  end subroutine
end

