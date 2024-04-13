module module_computeRate
    implicit none
    real*8 :: kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc
    contains

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
        krt = 0.D0
        knrt = 0.D0
        krs = kp * phiPF
        kisc = kp * phiDF / (phiPF + phiDF)
        knrs = kp - krs - kisc
        krisc = kd * (phiPF + phiDF) / phiDF
    end subroutine

    subroutine Kaji_NOnrt_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        implicit none
        real*8 :: kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc, phiPLQY
        krt = 0.D0
        knrt = 0.D0
        krisc = (kp + kd) / 2.D0 - SQRT(((kp + kd) / 2.D0) ** 2 - kp * kd * (1.D0 + phiDF / phiPF))
        kisc = kp * kd * phiDF / (krisc * phiPF)
        krs = kp * kd  * (phiPF + phiDF) / krisc
        knrs = kp * kd  * (1.D0 - phiPF - phiDF) / krisc
    end subroutine
end module
