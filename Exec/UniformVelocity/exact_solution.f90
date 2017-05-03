module exact_sol_module

  use ml_layout_module
  use bl_constants_module

  implicit none

  private

  public :: exact_sol_on_level, exact_sol

contains

  subroutine exact_sol_on_level(exact_phi,prob_lo,prob_hi,dx,time)

    type(multifab) , intent(inout) :: exact_phi(:)
    real(kind=dp_t), intent(in   ) :: dx,time
    real(kind=dp_t), intent(in   ) :: prob_lo(exact_phi(1)%dim), prob_hi(exact_phi(1)%dim)

    ! local
    integer :: i,dm,ng
    integer :: lo(exact_phi(1)%dim), hi(exact_phi(1)%dim)

    real(kind=dp_t), pointer :: df(:,:,:,:)
    
    do i=1,nfabs(exact_phi(1))
       df => dataptr(exact_phi(1),i)
       lo = lwb(get_box(exact_phi(1),i))
       hi = upb(get_box(exact_phi(1),i))
       select case(dm)
       case (2)
          call exact_sol_2d(df(:,:,1,1), ng, lo, hi, prob_lo, prob_hi, dx, time)
       case (3)
          call exact_sol_3d(df(:,:,:,1), ng, lo, hi, prob_lo, prob_hi, dx, time)
       end select
    end do

  end subroutine exact_sol_on_level

  subroutine exact_sol(mla,exact_phi,prob_lo,prob_hi, dx,time)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: exact_phi(:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    real(kind=dp_t), intent(in   ) :: prob_lo(exact_phi(1)%dim), prob_hi(exact_phi(1)%dim)

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n

    real(kind=dp_t), pointer :: df(:,:,:,:)

    ng = exact_phi(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       do i=1,nfabs(exact_phi(n))
          df => dataptr(exact_phi(n),i)
          lo = lwb(get_box(exact_phi(n),i))
          hi = upb(get_box(exact_phi(n),i))
          select case(dm)
          case (2)
          call exact_sol_2d(df(:,:,1,1), ng, lo, hi, prob_lo, prob_hi, dx(n), time)
          case (3)
          call exact_sol_3d(df(:,:,:,1), ng, lo, hi, prob_lo, prob_hi, dx(n), time)
          end select
       end do
    end do

  end subroutine exact_sol
  
  subroutine exact_sol_2d(phi, ng, lo, hi, prob_lo, prob_hi, dx, time)

    integer          :: lo(2), hi(2), ng
    double precision :: phi(lo(1)-ng:,lo(2)-ng:)
    double precision :: prob_lo(2),prob_hi(2)
    double precision :: dx, time 


    integer          :: i,j
    double precision :: pi2,vavg,x,y,length
    pi2 = 2.d0*3.14159265359d0
    vavg = 0.d0
    length = 1.d0


    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j) + .5d0) * dx
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+ .5d0) * dx

! velx set both x and y in the center 
          phi(i,j) = sin(.5d0*pi2*x)!vavg - 2.d0*(cos(pi2*(x-vavg*time)/length)*sin(pi2*(y-vavg*time)/length))

!          print*,i,j,phi(i,j)

       end do
    end do

    end subroutine exact_sol_2d

    subroutine exact_sol_3d(phi, ng, lo, hi, prob_lo,prob_hi, dx,time)

    integer          :: lo(3), hi(3), ng
    double precision :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: prob_lo(3),prob_hi(3)
    double precision :: dx,time

    ! local varables
    integer          :: i,j,k
    double precision :: x,y,z,r1,tmp

    !$omp parallel do private(i,j,k,x,y,z,r1)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+0.5d0) * dx - time  !velocity
       if (z .lt. prob_lo(3)) z = z + (prob_hi(3) - prob_lo(3))
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0) * dx - time !velocity
          if (y .lt. prob_lo(2)) y = y + (prob_hi(2) - prob_lo(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0) * dx - time !velocity
             if (x .lt. prob_lo(1)) x = x + (prob_hi(1) - prob_lo(1))

             r1 = ((x)**2 + (y)**2 + (z)**2) / 0.01d0

             phi(i,j,k) = 1.d0 + exp(-r1)

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine exact_sol_3d



end module exact_sol_module
