! 
! This routines calcuates the directional derivatives of given data using FFTW
! 
!
! Created by G.P. Brandino, I. Girotto, R. Gebauer
! Last revision: March 2016
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine derivative(n1,n2,n3,L1,L2,L3,ipol,data,deriv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !
  ! Calculate the derivative in direction ipol of the array 'data'
  !
  use fft_wrapper 

  implicit none
 
  integer, intent(in) :: n1,n2,n3
  real(8), intent(in) :: L1,L2,L3
  integer, intent(in) :: ipol
  real(8), dimension(n1, n2, n3), intent(inout) :: data
  real(8), dimension(n1, n2, n3), intent(out) :: deriv


  real(8), parameter :: pi = 3.14159265358979323846d0
  
  logical :: use_3d_fft
  complex(8), dimension(:,:,:), allocatable :: aux
  complex(8), dimension(:), allocatable :: aux_1d
  real(8), dimension(:), allocatable :: aux_real_1d
  integer :: status
  real(8) :: G
  integer :: i,i1,i2,i3


  ! 
  ! This implementation will perform
  ! a full 3D FFT of the data, and then derive.
  ! 

  allocate(aux(n1, n2, n3),stat=status)
  call errore('derivative','error allocating array aux',status)
  !
  ! First get the FFT of data
  !
  call fft_3d(n1, n2, n3, data, aux, .true.)

  if(ipol.eq.1) then
     G = 2.d0*pi/L1
     do i1=1,n1
        i = i1-1
        if(i1-1.gt.n1/2) i = i-n1
        if(i1-1.eq.n1/2) i = 0 
        aux(i1,:,:) = aux(i1,:,:) * cmplx(0.d0,G*i,8)
     enddo
  else if(ipol.eq.2) then
     G = 2.d0*pi/L2
     do i1=1,n2
        i = i1-1
        if(i1-1.gt.n2/2) i = i-n2
        if(i1-1.eq.n2/2) i = 0 
        aux(:,i1,:) = aux(:,i1,:) * cmplx(0.d0,G*i,8)
     enddo
  else if(ipol.eq.3) then
     G = 2.d0*pi/L3
     do i1=1,n3
        i = i1-1
        if(i.gt.n3/2) i = i-n3
        if(i.eq.n3/2) i = 0 
        aux(:,:,i1) = aux(:,:,i1) * cmplx(0.d0,G*i,8)
     enddo
  else
     call errore('derivative','wrong value for ipol',1)
  endif
  !
  ! Now go back to real space
  !
  call fft_3d(n1, n2, n3, deriv, aux, .false.)

end subroutine derivative
  
