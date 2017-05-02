Module fft_wrapper

  use, intrinsic :: iso_c_binding

  implicit none

  include 'fftw3-mpi.f03'
  include 'mpif.h'

  type(C_PTR) :: fw_plan, bw_plan
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), pointer :: fftw_data
  integer(C_INTPTR_T) :: alloc_local
  integer :: local_n3, local_n3_offset
  type(C_PTR) :: p

contains

  subroutine init_fft( n1, n2, n3 )

    implicit none

    integer, intent(in) :: n1, n2, n3
    integer(C_INTPTR_T) :: n1_c,n2_c,n3_c,local_n3_c, local_n3_offset_c
    n1_c=n1
    n2_c=n2
    n3_c=n3

    ! Allocate aligned memory
    ! See details here:
    ! http://www.fftw.org/fftw3_doc/Allocating-aligned-memory-in-Fortran.html#Allocating-aligned-memory-in-Fortran
    ! HINT: the allocation size of the c and fortran pointer, is given by fftw_mpi_local_size_3d
    alloc_local = fftw_mpi_local_size_3d(n3_c, n2_c, n1_c, MPI_COMM_WORLD, local_n3_c, local_n3_offset_c)
    local_n3=local_n3_c
    local_n3_offset=local_n3_offset_c
    p = fftw_alloc_complex(alloc_local)
    call c_f_pointer(p, fftw_data, [ n1_c, n2_c, local_n3_c ])

    !!http://www.fftw.org/doc/FFTW-MPI-Fortran-Interface.html
    fw_plan = fftw_mpi_plan_dft_3d( n3_c, n2_c, n1_c, fftw_data, fftw_data, MPI_COMM_WORLD,FFTW_FORWARD, FFTW_ESTIMATE )

    !!http://www.fftw.org/doc/FFTW-MPI-Fortran-Interface.html
    bw_plan = fftw_mpi_plan_dft_3d( n3_c, n2_c, n1_c, fftw_data, fftw_data, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE )

  end subroutine init_fft


  subroutine close_fft( )

    implicit none

    call fftw_destroy_plan(bw_plan)
    call fftw_destroy_plan(fw_plan)
    call fftw_free(p)


  end subroutine close_fft


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fft_3d(n1,n2,n3, data_direct, data_rec, direct_to_reciprocal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! This subroutine uses fftw to calculate 3-dimensional discrete FFTs.
  ! The data in direct space is assumed to be real-valued
  ! The data in reciprocal space is complex. 
  ! direct_to_reciprocal indicates in which direction the FFT is to be calculated
  ! 
  ! Note that for real data in direct space (like here), we have
  ! F(N-j) = conj(F(j)) where F is the array in reciprocal space.
  ! Here, we do not make use of this property.
  ! Also, we do not use the special (time-saving) routines of FFTW which
  ! allow one to save time and memory for such real-to-complex transforms.
  !
  ! f: array in direct space
  ! F: array in reciprocal space
  ! 
  ! F(k) = \sum_{l=0}^{N-1} exp(- 2 \pi I k*l/N) f(l)
  ! f(l) = 1/N \sum_{k=0}^{N-1} exp(+ 2 \pi I k*l/N) F(k)
  !
  ! See also: http://www.fftw.org/fftw3_doc/Calling-FFTW-from-Modern-Fortran.html#Calling-FFTW-from-Modern-Fortran
  !

  implicit none
  
  integer, intent(in) :: n1,n2,n3
  real(C_DOUBLE), dimension(n1,n2,local_n3), intent(inout) :: data_direct
  complex(C_DOUBLE_COMPLEX), dimension(n1,n2,local_n3), intent(inout) :: data_rec
  logical, intent(in) :: direct_to_reciprocal

  real(C_DOUBLE) :: fac

  ! Allocate aligned memory
  ! See details here:
  ! http://www.fftw.org/fftw3_doc/Allocating-aligned-memory-in-Fortran.html#Allocating-aligned-memory-in-Fortran

  !
  ! Now distinguish in which direction the FFT is performed
  !
  if(direct_to_reciprocal) then


      fftw_data(:,:,:) = cmplx(data_direct(:,:,:),0.d0,C_DOUBLE_COMPLEX)

      call fftw_mpi_execute_dft(fw_plan, fftw_data, fftw_data)

      data_rec(:,:,:) = fftw_data(:,:,:)

  else

     fftw_data(:,:,:) = data_rec(:,:,:)
     
     call fftw_mpi_execute_dft(bw_plan, fftw_data, fftw_data)

     data_direct(:,:,:) = real(fftw_data(:,:,:), C_DOUBLE)
     
     fac = 1.d0/real(n1*n2*n3,8)

     data_direct(:,:,:) = data_direct(:,:,:) * fac

  endif

end subroutine fft_3d
 
end Module fft_wrapper
                       
