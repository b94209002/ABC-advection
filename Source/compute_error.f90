module compute_error_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use ml_restrict_fill_module

  implicit none

  private

  public :: compute_error_on_level, compute_error,compute_error_from_refinement

contains

  subroutine compute_error_on_level(error,phi,phi_exact,dx)

    type(multifab) , intent(inout) :: phi,phi_exact
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(inout) :: error

    ! local
    integer i,ng,dm
    integer :: lo(phi%dim), hi(phi%dim)
    type(multifab) :: res

    real(kind=dp_t), pointer :: dp(:,:,:,:),de(:,:,:,:),dr(:,:,:,:)

    ng = phi%ng
    dm = phi%dim

    call multifab_build(res,phi%la,1,0)

    do i=1,nfabs(phi)
       dp => dataptr(phi,i)
       de => dataptr(phi_exact,i)
       dr => dataptr(res,i)
       lo = lwb(get_box(phi,i))
       hi = upb(get_box(phi,i))
       select case(dm)
       case (2)
          call compute_error_2d(dp(:,:,1,1), de(:,:,1,1), dr(:,:,1,1), ng, lo, hi, dx)
       case (3)
          call compute_error_3d(dp(:,:,:,1), de(:,:,:,1), dr(:,:,:,1), ng, lo, hi, dx)
       end select
    end do
    error = norm_inf(res)

    call multifab_destroy(res)

  end subroutine compute_error_on_level
  
  subroutine compute_error(error, mla,phi,phi_exact,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:),phi_exact(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(inout) :: error
    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n
    type(multifab) :: res

    real(kind=dp_t), pointer :: dp(:,:,:,:),de(:,:,:,:),dr(:,:,:,:)

    ng = phi(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel
! here might have error since nlev = 1 
    call multifab_build(res,phi(nlevs)%la,1,0)

    do n=1,nlevs

       do i=1,nfabs(phi(n))
          dp => dataptr(phi(n),i)
          de => dataptr(phi_exact(n),i)
          dr => dataptr(res,i) 
          lo = lwb(get_box(phi(n),i))
          hi = upb(get_box(phi(n),i))
          select case(dm)
          case (2)
             call compute_error_2d(dp(:,:,1,1), de(:,:,1,1), dr(:,:,1,1), ng, lo, hi, dx(n))
          case (3)
             call compute_error_3d(dp(:,:,:,1), de(:,:,:,1), dr(:,:,:,1), ng, lo, hi, dx(n))
          end select
       end do
       
    end do

    error = norm_inf(res)
    call multifab_destroy(res)

  end subroutine compute_error

  subroutine compute_error_from_refinement(error, mla,phi,dx)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(inout) :: error
    ! local variables
    integer :: lo(mla%dim), hi(mla%dim), lof(mla%dim), hif(mla%dim)
    integer :: nlevs, dm, ng, i, n, it 
    type(multifab) :: res, phi_exact

    real(kind=dp_t), pointer :: dp(:,:,:,:),dpf(:,:,:,:),dr(:,:,:,:),de(:,:,:,:)

    ng = phi(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel
! here might have error since nlev = 1 
    call multifab_build(res,phi(1)%la,1,0)
    call multifab_build(phi_exact,phi(1)%la,1,0) 



    do n=1,1
       do i=1,nfabs(phi(n+1))
          it = (i-1)/2**dm +1
!          it = mod(i-1,4) + 1
          de => dataptr(phi_exact,it)
          dpf => dataptr(phi(n+1),i)
          lo = lwb(get_box(phi(n),it))
          hi = upb(get_box(phi(n),it))
          lof = lwb(get_box(phi(n+1),i))
          hif = upb(get_box(phi(n+1),i))

          select case(dm)
          case (2)
             call compute_exact_from_refinement_2d(de(:,:,1,1), lo, hi, & 
                                                  dpf(:,:,1,1), lof, hif, ng)
          case (3)
             call compute_exact_from_refinement_3d(de(:,:,:,1), lo, hi, & 
                                                  dpf(:,:,:,1), lof, hif, ng)
          end select
       end do

    if (.true.) then
       do i=1,nfabs(phi(n))
          dp => dataptr(phi(n),i)
          de => dataptr(phi_exact,i)
          dr => dataptr(res,i)
          lo = lwb(get_box(phi(n),i))
          hi = upb(get_box(phi(n),i))

          select case(dm)
          case (2)
             call compute_error_2d(dp(:,:,1,1), de(:,:,1,1), dr(:,:,1,1), ng, lo, hi, dx(n))
          case (3)
             call compute_error_3d(dp(:,:,:,1), de(:,:,:,1), dr(:,:,:,1), ng, lo, hi, dx(n))
          end select
       end do

    endif 

    end do

    error = norm_l2(res)
    call multifab_destroy(res)
    call multifab_destroy(phi_exact)


  end subroutine compute_error_from_refinement


  subroutine compute_error_2d(phi, phi_exact, res, ng, lo, hi, dx)

    integer          :: lo(2), hi(2), ng
    double precision :: phi(lo(1)-ng:,lo(2)-ng:), phi_exact(lo(1):,lo(2):), res(lo(1):,lo(2):)
    double precision :: dx
 
    ! local varables
    integer          :: i,j
    double precision :: r1

    r1 = sqrt(dx*dx)

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          ! compute error in 2 norm 
          res(i,j) =r1*( phi(i,j) -phi_exact(i,j))
          ! compute error in inf-norm            
          ! res(i,j) = abs(phi(i,j) - phi_exact(i,j))

          ! if (res(i,j) .gt. 1.d0) Print*,i,j,phi(i,j),phi_exact(i,j)
          !print*,i,j,phi(i,j),phi_exact(i,j)
       end do
    end do

    end subroutine compute_error_2d

    subroutine compute_error_3d(phi, phi_exact, res, ng, lo, hi, dx)

    integer          :: lo(3), hi(3), ng
    double precision :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:),phi_exact(lo(1):,lo(2):,lo(3):), res(lo(1):,lo(2):,lo(3):)
    double precision :: dx
 
    ! local varables
    integer          :: i,j,k
    double precision :: r1

    !$omp parallel do private(i,j,k,r1)
    r1 = sqrt(dx*dx*dx)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             res(i,j,k) = r1*( phi(i,j,k) -phi_exact(i,j,k))

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine compute_error_3d

  subroutine compute_exact_from_refinement_2d(phi_exact, lo, hi, phi2, lo2, hi2, ng)

    integer          :: lo(2), hi(2), lo2(2), hi2(2), ng
    double precision :: phi_exact(lo(1):,lo(2):)
    double precision :: phi2(lo2(1)-ng:,lo2(2)-ng:)

    ! local varables
    integer          :: i,j,it,jt



    do j=lo2(2),hi2(2),2
       jt = j/2
       do i=lo2(1),hi2(1),2
          it = i/2
          ! compute error in 2 norm
          phi_exact(it,jt) = .25d0*(phi2(i,j) + phi2(i+1,j) + phi2(i,j+1) + phi2(i+1,j+1))
       end do
    end do

  end subroutine compute_exact_from_refinement_2d
  subroutine compute_exact_from_refinement_3d(phi_exact, lo, hi, phi2, lo2, hi2, ng)

    integer          :: lo(3), hi(3), lo2(3), hi2(3), ng
    double precision :: phi_exact(lo(1):,lo(2):,lo(3):)
    double precision :: phi2(lo2(1)-ng:,lo2(2)-ng:,lo2(3)-ng:)

    ! local varables
    integer          :: i,j,k,it,jt,kt

    !$omp parallel do private(i,j,k,r1)
    do k=lo2(3),hi2(3),2
       kt = k/2 
       do j=lo2(2),hi2(2),2
          jt = j/2
          do i=lo2(1),hi2(1),2
             it = i/2
             phi_exact(it,jt,kt) =.125d0*(phi2(i,j,k) + phi2(i+1,j,k) + phi2(i,j+1,k) + phi2(i,j,k+1) +  &
                        phi2(i+1,j+1,k) + phi2(i,j+1,k+1) + phi2(i+1,j,k+1) + phi2(i+1,j+1,k+1))  
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine compute_exact_from_refinement_3d

end module compute_error_module
