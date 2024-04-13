module module_print
    contains

    subroutine printTADF(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
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

    subroutine printFormula(mainOpt)
        implicit none
        character :: mainOpt
        write(*,*)
        select case (mainOpt)
            case ('1')
                write(*,*) "+-------------------------------------------------------+"
                write(*,*) "|  krS   = kp * Phi_PF                                  |"
                write(*,*) "|  knrS  = 0                                            |"
                write(*,*) "|  kISC  = kp * (1 - Phi_PF)                            |"
                write(*,*) "|  krT   = 0                                            |"
                write(*,*) "|  knrT  = kd * ((1 - Phi_PF - Phi_DF) / (1 - Phi_PF))  |"
                write(*,*) "|  kRISC = (kd * Phi_DF) / (Phi_PF * (1 - Phi_PF))      |"
                write(*,*) "+-------------------------------------------------------+"
            case ('2')
                write(*,*) "+-------------------------------------------------------------------+"
                write(*,*) "|  krS   = kp * Phi_PF                                              |"
                write(*,*) "|  knrS  = kp * Phi_PF * (1 - Phi_PF - Phi_DF) / (Phi_PF + Phi_DF)  |"
                write(*,*) "|  kISC  = kp * Phi_DF / (Phi_PF + Phi_DF)                          |"
                write(*,*) "|  krT   = 0                                                        |"
                write(*,*) "|  knrT  = 0                                                        |"
                write(*,*) "|  kRISC = kd * (Phi_PF + Phi_DF) / Phi_PF                          |"
                write(*,*) "+-------------------------------------------------------------------+"
            case ('3')
                write(*,*) "+----------------------------------------------------------------------+"
                write(*,*) "|  krS   = kp * Phi_PF                                                 |"
                write(*,*) "|  knrS  = kp * Phi_PF * (1 - Phi_PF - Phi_DF) / (Phi_PF + Phi_DF)     |"
                write(*,*) "|  kISC  = (kp * Phi_DF / (Phi_PF + Phi_DF)) - (kd * Phi_DF / Phi_PF)  |"
                write(*,*) "|  krT   = 0                                                           |"
                write(*,*) "|  knrT  = 0                                                           |"
                write(*,*) "|  kRISC = kd * (Phi_PF + Phi_DF) / Phi_PF                             |"
                write(*,*) "+----------------------------------------------------------------------+"
            case ('4')
                write(*,*) "+-------------------------------------------------------------------+"
                write(*,*) "|  krS   = kp * Phi_PF                                              |"
                write(*,*) "|  knrS  = kp * Phi_PF * (1 - Phi_PF - Phi_DF) / (Phi_PF + Phi_DF)  |"
                write(*,*) "|  kISC  = kp * Phi_DF / (Phi_PF + Phi_DF)                          |"
                write(*,*) "|  krT   = 0                                                        |"
                write(*,*) "|  knrT  = 0                                                        |"
                write(*,*) "|  kRISC = kd * (Phi_PF + Phi_DF) / Phi_DF                          |"
                write(*,*) "+-------------------------------------------------------------------+"
            case ('5')
                write(*,*) "+---------------------------------------------------------------------------------------+"
                write(*,*) "|  krS   = kp * kd  * (Phi_PF + Phi_DF) / kRISC                                         |"
                write(*,*) "|  knrS  = kp * kd  * (1 - Phi_PF - Phi_DF) / kRISC                                     |"
                write(*,*) "|  kISC  = kp * kd * Phi_DF / (kRISC * Phi_PF)                                          |"
                write(*,*) "|  krT   = 0                                                                            |"
                write(*,*) "|  knrT  = 0                                                                            |"
                write(*,*) "|  kRISC = (kp + kd) / 2 - SQRT(((kp + kd) / 2) ^ 2 - kp * kd * (1 + Phi_DF / Phi+PF))  |"
                write(*,*) "+---------------------------------------------------------------------------------------+"
        end select
        write(*,*)
    end subroutine

    subroutine printFitMethod()
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
    end subroutine
end module
