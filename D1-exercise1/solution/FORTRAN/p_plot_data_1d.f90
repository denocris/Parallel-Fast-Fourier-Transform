!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine plot_data_1d(filename,n1, n2, n3, local_n3, local_n3_offset, idir, data)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none

  include 'mpif.h'

  integer, intent(in) :: n1, n2, n3, idir, local_n3, local_n3_offset
  real(8), intent(in) :: data(n1,n2,local_n3)
  character(len=*), intent(in) :: filename

  character(len=100) :: fnme
  character(len=4) :: fnbr
  integer :: fnumber
  logical :: exst
  real(8) :: local_n3_buffer(local_n3)
  real(8) :: buffer(n3)
  integer :: i1,i2,i3,i,owner, rank, ierror, size 
  integer, dimension(:), allocatable :: sizes, displ

  fnumber = 1
  file_loop: do
     write(fnbr,'(i4.4)') fnumber
     fnme = trim(filename)//trim(fnbr)//".dat"
     inquire(file=fnme, exist=exst)
     if(.not.exst) exit file_loop
     fnumber = fnumber+1
     if(fnumber.ge.10000) call errore('plot_data_1d','no free filename found',fnumber)
  enddo file_loop


  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  owner=size+1;
  if ( (n3/2 .gt. local_n3_offset) .and. (n3/2 .le. local_n3_offset+local_n3)  ) then 
     owner=rank
  end if 

  if(idir.eq.1) then
     i2=n2/2
     i3=n3/2
     if (rank .eq. owner) then 
        open(unit=11,file=fnme, status='new')
        do i1=1,n1
           write(11,fmt='(f14.6,"  ")') data(i1,i2,i3-local_n3_offset)
        enddo
        close(11)
     end if
  else if(idir.eq.2) then
     i1=n1/2
     i3=n3/2
     if (rank .eq. owner) then
        open(unit=11,file=fnme, status='new')
        do i2=1,n2
           write(11,fmt='(f14.6,"  ")') data(i1,i2,i3-local_n3_offset)
        enddo
        close(11)
     end if 
  else if(idir.eq.3) then
     i1=n1/2
     i2=n2/2
     allocate(sizes(size))
     allocate(displ(size))
     call MPI_Gather(local_n3, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)
     if (rank .eq. 0) then
         displ=0
         do i=2,size
            displ(i)=sizes(i-1)+displ(i-1)
         end do
     end if
     local_n3_buffer=data(i1,i2,:)
     call MPI_Gatherv(local_n3_buffer,local_n3, MPI_DOUBLE,buffer,sizes, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierror)
     deallocate(sizes)     
     if (rank .eq. 0) then 
        open(unit=11,file=fnme, status='new')
        do i3=1,n3
           write(11,fmt='(f14.6,"  ")') buffer(i3)
        enddo
        close(11)
     end if 

  else
     call errore('plot_data_1d','wrong value for idir',1)
  endif

  call MPI_Barrier(MPI_COMM_WORLD,ierror) 
 
end subroutine plot_data_1d
