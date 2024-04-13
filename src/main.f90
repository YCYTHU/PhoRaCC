! PhoRaCC: Simplify the calculation of photophysical rate constants
! Programmed by Yaocy

program main
  use module_computeRate
  use module_fitDecayData
  use module_print

  character :: opt = '0'

  write(*,*) "PhoRaCC: Simplify the calculation of photophysical rate constants"
  write(*,*) "Programmed by Yaocy"
  write(*,*) "Release date: 2024-Apr-11"
  write(*,*)

  do while (.true.)
    write(*,*) "                   ************ Three-State TADF ************"
    write(*,*) "1. krS + knrS + kISC >> kRISC & knrS = 0 & krT = 0 (DOI: 10.1021/acs.jpca.1c04056)"
    write(*,*) "2. krS + knrS + kISC >> kRISC & knrT = 0 & krT = 0 (DOI: 10.1021/acs.jpca.1c04056)"
    write(*,*) "3. krT = 0 (DOI: 10.1021/acs.jpca.1c04056)"
    write(*,*) "4. Phi_DF / Phi_PF > 4 & krT = knrT = 0 (DOI: 10.1088/2050-6120/aa537e)"
    write(*,*) "5. krT = knrT = 0 (DOI: 10.1038/s41566-020-0667-0)"
    read(*,*) opt
    select case (opt)
      case ('1')
        write(*,*)
        write(*,*) "Tsuchiya, Y., Diesing, S., ... & Adachi, C. (2021). The Journal of Physical Chemistry A, 125(36), 8074-8089."
        write(*,*) "(DOI: 10.1021/acs.jpca.1c04056)"
        write(*,*)
        call getParameters(kp, kd, phiPF, phiDF, doCompute, opt)
        if (.not. doCompute) then
          cycle
        end if
        call SSA_NOnrs_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printTADF(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
      case ('2')
        write(*,*)
        write(*,*) "Tsuchiya, Y., Diesing, S., ... & Adachi, C. (2021). The Journal of Physical Chemistry A, 125(36), 8074-8089."
        write(*,*) "(DOI: 10.1021/acs.jpca.1c04056)"
        write(*,*)
        call getParameters(kp, kd, phiPF, phiDF, doCompute, opt)
        if (.not. doCompute) then
          cycle
        end if
        call SSA_NOnrt_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printTADF(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
      case ('3')
        write(*,*)
        write(*,*) "Tsuchiya, Y., Diesing, S., ... & Adachi, C. (2021). The Journal of Physical Chemistry A, 125(36), 8074-8089."
        write(*,*) "(DOI: 10.1021/acs.jpca.1c04056)"
        write(*,*)
        call getParameters(kp, kd, phiPF, phiDF, doCompute, opt)
        if (.not. doCompute) then
          cycle
        end if
        call NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printTADF(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
      case ('4')
        write(*,*)
        write(*,*) "Dias, F. B., Penfold, T. J., & Monkman, A. P. (2017). Methods and applications in fluorescence, 5(1), 012001."
        write(*,*) "(DOI: 10.1088/2050-6120/aa537e)"
        write(*,*)
        call getParameters(kp, kd, phiPF, phiDF, doCompute, opt)
        if (.not. doCompute) then
          cycle
        end if
        call Monkman_LargeDF(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printTADF(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
      case ('5')
        write(*,*)
        write(*,*) "Wada, Y., Nakagawa, H., ... & Kaji, H. (2020). Nature Photonics, 14(10), 643-649."
        write(*,*) "(DOI: 10.1038/s41566-020-0667-0)"
        write(*,*)
        call getParameters(kp, kd, phiPF, phiDF, doCompute, opt)
        if (.not. doCompute) then
          cycle
        end if
        call Kaji_NOnrt_NOrt(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
        call printTADF(kp, kd, phiPF, phiDF, krs, knrs, kisc, krt, knrt, krisc)
      case ('q')
        exit
      case default
        cycle
    end select
  enddo

  contains



end

