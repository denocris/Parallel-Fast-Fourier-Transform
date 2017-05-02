
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine errore(where, why, number)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

  character(len=*) :: where, why
  integer :: number

  if(number.eq.0) return

  write(6,*)
  write(6,'(80("-"))')
  write(6,*)

  if(number.lt.0) then
     write(6,'("WARNING: In routine ",a," "&
          &,a," :: ",i5)') where, why, number
     write(6,*)
     write(6,'(80("-"))')
     write(6,*)
  else
     write(6,&
          & '("ERROR: In routine ",a," "&
          &,a," :: ",i5)') where, why, number
     write(6,*)
     write(6,'(80("-"))')
     write(6,*)

     stop
  endif

end subroutine errore

