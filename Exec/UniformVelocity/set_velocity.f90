module set_velocity_module

  use ml_layout_module
  use bl_constants_module

  implicit none

  private

  public :: set_velocity_on_level, set_velocity, set_forcing

contains

  subroutine set_velocity_on_level(velocity,dx,time,prob_lo)

    type(multifab) , intent(inout) :: velocity(:)
    real(kind=dp_t), intent(in   ) :: dx,time
    real(kind=dp_t), intent(in   ) :: prob_lo(:)
    ! local
    integer :: i,dm,ng
    integer :: lo(velocity(1)%dim), hi(velocity(1)%dim)

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)
    
    do i=1,nfabs(velocity(1))
       dp1 => dataptr(velocity(1),i)
       dp2 => dataptr(velocity(2),i)
       lo = lwb(get_box(velocity(1),i))
       hi = upb(get_box(velocity(1),i))
       select case(dm)
       case (2)
          call set_velocity_2d(dp1(:,:,1,1), dp2(:,:,1,1), ng, &
                                   lo, hi, dx, time, prob_lo)
       case (3)
          dp3 => dataptr(velocity(3),i)
          call set_velocity_3d(dp1(:,:,:,1), dp2(:,:,:,1), dp3(:,:,:,1), ng, &
                                   lo, hi, dx, time, prob_lo)
       end select
    end do

    do i=1,dm
       call multifab_fill_boundary(velocity(i))
    end do

  end subroutine set_velocity_on_level


  subroutine set_velocity(mla,velocity,dx,time,prob_lo)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: velocity(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    real(kind=dp_t), intent(in   ) :: prob_lo(:)

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n

    real(kind=dp_t), pointer :: dp1(:,:,:,:)
    real(kind=dp_t), pointer :: dp2(:,:,:,:)
    real(kind=dp_t), pointer :: dp3(:,:,:,:)

    ng = velocity(1,1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       do i=1,nfabs(velocity(n,1))
          dp1 => dataptr(velocity(n,1),i)
          dp2 => dataptr(velocity(n,2),i)
          lo = lwb(get_box(velocity(n,1),i))
          hi = upb(get_box(velocity(n,1),i))
          select case(dm)
          case (2)
             call set_velocity_2d(dp1(:,:,1,1), dp2(:,:,1,1), ng, &
                                      lo, hi, dx(n), time, prob_lo)
          case (3)
             dp3 => dataptr(velocity(n,3),i)
             call set_velocity_3d(dp1(:,:,:,1), dp2(:,:,:,1), dp3(:,:,:,1), ng, &
                                      lo, hi, dx(n), time, prob_lo)
          end select
       end do
    end do
    
    do n=1,nlevs
       do i=1,dm
          call multifab_fill_boundary(velocity(n,i))
       end do
    end do

  end subroutine set_velocity
  subroutine set_forcing(mla,forcing,dx,time,prob_lo)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: forcing(:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    real(kind=dp_t), intent(in   ) :: prob_lo(:)

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n

    real(kind=dp_t), pointer :: df(:,:,:,:)

    ng = forcing(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       do i=1,nfabs(forcing(n))
          df => dataptr(forcing(n),i)
          lo = lwb(get_box(forcing(n),i))
          hi = upb(get_box(forcing(n),i))
          select case(dm)
          case (2)
             call set_forcing_2d(df(:,:,1,1), ng, lo, hi, dx(n), time, prob_lo)
          case (3)
             call set_forcing_3d(df(:,:,:,1), ng, lo, hi, dx(n), time, prob_lo)
          end select
       end do
    end do

    do n=1,nlevs
       call multifab_fill_boundary(forcing(n))
    end do

  end subroutine set_forcing
  
  subroutine set_velocity_2d(velx, vely, ng, lo, hi, dx, time, prob_lo)

    integer          :: lo(2), hi(2), ng
    double precision :: velx(lo(1)-ng:,lo(2)-ng:)
    double precision :: vely(lo(1)-ng:,lo(2)-ng:)
    double precision :: dx,time
    double precision :: prob_lo(2) 
   
    integer          :: i,j
    double precision :: pi2,vavg,x,y,h,length
    pi2 = 2.d0*3.14159265359d0
    vavg = 0.d0  
    length = 1.d0

    h = .5d0*dx

    do j=lo(2)-1,hi(2)+1
       y = prob_lo(2) + (dble(j) + .5d0) * dx
       do i=lo(1)-1,hi(1)+2
          x = prob_lo(1) + (dble(i)) * dx

! velx set x in the edge and y in the center 
          velx(i,j) = vavg - 2.d0*cos(pi2*(x-vavg*time)/length)*sin(pi2*(y-vavg*time)/length)
!          print*,i,j,x,y,velx(i,j)
       end do
    end do
    do j=lo(2)-1,hi(2)+2
       y = prob_lo(2) + (dble(j)) * dx
       do i=lo(1)-1,hi(1)+1
          x = prob_lo(1) + (dble(i) + .5d0) * dx

! velx set x in the edge and y in the center 
          vely(i,j) = vavg + 2.d0*sin(pi2*(x-vavg*time)/length)*cos(pi2*(y-vavg*time)/length)
!          print*,i,j,x,y,velx(i,j)
       end do
    end do

!  velx =  1.d0
!  vely = -1.d0 

  end subroutine set_velocity_2d

  subroutine set_velocity_3d(velx, vely, velz, ng, lo, hi, dx, time, prob_lo)

    integer          :: lo(3), hi(3), ng
    double precision :: velx(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: vely(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: velz(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: dx,time, prob_lo(3)

    ! Constant velocity field
    velx = 1.d0
    vely = 1.d0
    velz = 1.d0

  end subroutine set_velocity_3d

  subroutine set_forcing_2d(f, ng, lo, hi, dx, time, prob_lo)

    integer          :: lo(2), hi(2), ng
    double precision :: f(lo(1)-ng:,lo(2)-ng:)
    double precision :: dx,time,prob_lo(2)

    integer          :: i,j
    double precision :: pi4,vavg,x,y,h,length
    pi4 = 4.d0*3.14159265359d0
    vavg = 0.d0
    length = 1.d0


    do j=lo(2)-1,hi(2)+1
       y = prob_lo(2) + (dble(j)+ .5d0) * dx
       do i=lo(1)-1,hi(1)+1
          x = prob_lo(1) + (dble(i)+ .5d0) * dx

          f(i,j) =  -1.d0*pi4/length*sin(pi4*(x-vavg*time)/length)
          !print*,i,j,f(i,j)
       end do
    end do

  end subroutine set_forcing_2d

  subroutine set_forcing_3d(f, ng, lo, hi, dx, time,prob_lo)

    integer          :: lo(3), hi(3), ng
    double precision :: f(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    double precision :: dx,time,prob_lo(3)

    ! Constant forcing field
    f = 0.d0

  end subroutine set_forcing_3d


end module set_velocity_module
