
module update_phi_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use ml_restrict_fill_module
  use bndry_reg_module

  implicit none

  private

  public :: update_phi, update_phi_single_level, update_bds_phi, update_bds_phi_single_level

contains

  subroutine update_phi(mla,phi_old,phi_new,flux,source,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:)
    type(multifab) , intent(inout) :: phi_new(:)
    type(multifab) , intent(in   ) :: flux(:,:)
    type(multifab) , intent(in   ) :: source(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: n, nlevs 

    nlevs = mla%nlevel

    do n=1,nlevs

       call update_phi_single_level(mla,phi_old(n),phi_new(n),flux(n,:), source(n),dx(n),dt(n))

    end do

    ! restrict the multi-level data, and
    ! fill all boundaries: same-level, coarse-fine, periodic, and domain boundaries 
    call ml_restrict_and_fill(nlevs, phi_new, mla%mba%rr, the_bc_tower%bc_tower_array)

  end subroutine update_phi

  subroutine update_phi_single_level(mla,phi_old,phi_new,flux,source,dx,dt)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi_old
    type(multifab) , intent(inout) :: phi_new
    type(multifab) , intent(in   ) :: flux(:)
    type(multifab) , intent(in   ) :: source
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: dt

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, ng_s, i

    real(kind=dp_t), pointer ::  po(:,:,:,:)
    real(kind=dp_t), pointer ::  pn(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)
    real(kind=dp_t), pointer :: src(:,:,:,:)

    dm    = mla%dim

    ng_p = phi_old%ng
    ng_f = flux(1)%ng
    ng_s = source%ng

    do i=1,nfabs(phi_old)
       po  => dataptr(phi_old,i)
       pn  => dataptr(phi_new,i)
       fxp => dataptr(flux(1),i)
       fyp => dataptr(flux(2),i)
       src => dataptr(source,i)
       lo = lwb(get_box(phi_old,i))
       hi = upb(get_box(phi_old,i))
       select case(dm)
       case (2)
          call update_phi_2d(po(:,:,1,1), pn(:,:,1,1), ng_p, &
                             fxp(:,:,1,1), fyp(:,:,1,1), ng_f, &
                             src(:,:,1,1), ng_s, lo, hi, dx, dt)
       case (3)
          fzp => dataptr(flux(3),i)
          call update_phi_3d(po(:,:,:,1), pn(:,:,:,1), ng_p, &
                             fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                             src(:,:,:,1), ng_s, lo, hi, dx, dt)
       end select
    end do

  end subroutine update_phi_single_level

  subroutine update_bds_phi(mla,phi_old,phi_new,flux,velocity,source,dx,dt,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: phi_old(:)
    type(multifab) , intent(inout) :: phi_new(:)
    type(multifab) , intent(in   ) :: flux(:,:)
    type(multifab) , intent(in   ) :: velocity(:,:)
    type(multifab) , intent(in   ) :: source(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: dt(:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: n, nlevs

    nlevs = mla%nlevel

    do n=1,nlevs

       call update_bds_phi_single_level(mla,phi_old(n),phi_new(n),flux(n,:),velocity(n,:),source(n), dx(n),dt(n))

    end do

    ! restrict the multi-level data, and
    ! fill all boundaries: same-level, coarse-fine, periodic, and domain
    ! boundaries 
    call ml_restrict_and_fill(nlevs, phi_new, mla%mba%rr, the_bc_tower%bc_tower_array)

  end subroutine update_bds_phi

  subroutine update_bds_phi_single_level(mla,phi_old,phi_new,flux,velocity,source,dx,dt)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi_old
    type(multifab) , intent(inout) :: phi_new
    type(multifab) , intent(in   ) :: flux(:)
    type(multifab) , intent(in   ) :: velocity(:)
    type(multifab) , intent(in   ) :: source
    real(kind=dp_t), intent(in   ) :: dx
    real(kind=dp_t), intent(in   ) :: dt

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: dm, ng_p, ng_f, ng_v, ng_s, i

    real(kind=dp_t), pointer ::  po(:,:,:,:)
    real(kind=dp_t), pointer ::  pn(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)
    real(kind=dp_t), pointer ::  up(:,:,:,:)
    real(kind=dp_t), pointer ::  vp(:,:,:,:)
    real(kind=dp_t), pointer ::  wp(:,:,:,:)
    real(kind=dp_t), pointer :: src(:,:,:,:)

    dm    = mla%dim

    ng_p = phi_old%ng
    ng_f = flux(1)%ng
    ng_v = velocity(1)%ng

    do i=1,nfabs(phi_old)
       po  => dataptr(phi_old,i)
       pn  => dataptr(phi_new,i)
       fxp => dataptr(flux(1),i)
       fyp => dataptr(flux(2),i)
       up  => dataptr(velocity(1),i) 
       vp  => dataptr(velocity(2),i)
       src => dataptr(source,i)
       lo = lwb(get_box(phi_old,i))
       hi = upb(get_box(phi_old,i))
       select case(dm)
       case (2)
          call update_bds_phi_2d(po(:,:,1,1), pn(:,:,1,1), ng_p, &
                             fxp(:,:,1,1), fyp(:,:,1,1), ng_f, &
                             up(:,:,1,1), vp(:,:,1,1), ng_v, &
                             src(:,:,1,1), ng_s, lo, hi, dx, dt)
       case (3)
          fzp => dataptr(flux(3),i)
          wp  => dataptr(velocity(3),i) 
          call update_bds_phi_3d(po(:,:,:,1), pn(:,:,:,1), ng_p, &
                             fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_f, &
                             up(:,:,:,1), vp(:,:,:,1), wp(:,:,:,1), ng_v, &
                             src(:,:,:,1), ng_s, lo, hi, dx, dt)
       end select
    end do

  end subroutine update_bds_phi_single_level

  subroutine update_phi_2d(phi_old, phi_new, ng_p, fluxx, fluxy, ng_f, source, ng_s, lo, hi, dx, dt)

    integer          :: lo(2), hi(2), ng_p, ng_f, ng_s
    double precision :: phi_old(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: phi_new(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision ::   fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision ::   fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision ::  source(lo(1)-ng_s:,lo(2)-ng_s:)
    double precision :: dx, dt

    ! local variables
    integer i,j
    double precision :: dtdx

    dtdx = dt/dx

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          phi_new(i,j) = phi_old(i,j) - &
               ( fluxx(i+1,j)-fluxx(i,j) &
                +fluxy(i,j+1)-fluxy(i,j) ) * dtdx 
          phi_new(i,j) = phi_new(i,j) + source(i,j)*dt
       end do
    end do

  end subroutine update_phi_2d

  subroutine update_phi_3d(phi_old, phi_new, ng_p, fluxx, fluxy, fluxz, ng_f, source, ng_s, lo, hi, dx, dt)

    integer          :: lo(3), hi(3), ng_p, ng_f, ng_s
    double precision :: phi_old(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: phi_new(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision ::   fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision ::   fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision ::   fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision ::  source(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    double precision :: dx, dt

    ! local variables
    integer i,j,k
    double precision :: dtdx

    dtdx = dt/dx

    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phi_new(i,j,k) = phi_old(i,j,k) - &
                  ( fluxx(i+1,j,k)-fluxx(i,j,k) &
                   +fluxy(i,j+1,k)-fluxy(i,j,k) &
                   +fluxz(i,j,k+1)-fluxz(i,j,k) ) * dtdx 
             phi_new(i,j,k) = phi_new(i,j,k) + source(i,j,k)*dt

          end do
       end do
    end do
    !$omp end parallel do

  end subroutine update_phi_3d

  subroutine update_bds_phi_2d(phi_old, phi_new, ng_p, fluxx, fluxy, ng_f, &
                               velx, vely, ng_v, source, ng_s, lo, hi, dx, dt)

    integer          :: lo(2), hi(2), ng_p, ng_f, ng_v, ng_s
    double precision :: phi_old(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision :: phi_new(lo(1)-ng_p:,lo(2)-ng_p:)
    double precision ::   fluxx(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision ::   fluxy(lo(1)-ng_f:,lo(2)-ng_f:)
    double precision ::    velx(lo(1)-ng_v:,lo(2)-ng_v:)
    double precision ::    vely(lo(1)-ng_v:,lo(2)-ng_v:)
    double precision ::  source(lo(1)-ng_s:,lo(2)-ng_s:)
    double precision :: dx, dt

    ! local variables
    integer i,j
    double precision :: dtdx

    dtdx = dt/dx

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          phi_new(i,j) = phi_old(i,j) - &
               ( velx(i+1,j)*fluxx(i+1,j)-velx(i,j)*fluxx(i,j) &
                +vely(i,j+1)*fluxy(i,j+1)-vely(i,j)*fluxy(i,j) ) * dtdx
          phi_new(i,j) = phi_new(i,j) + dt * source(i,j)
       end do
    end do

  end subroutine update_bds_phi_2d

  subroutine update_bds_phi_3d(phi_old, phi_new, ng_p, fluxx, fluxy, fluxz, ng_f, &
                               velx, vely, velz, ng_v, source, ng_s, lo, hi, dx, dt)

    integer          :: lo(3), hi(3), ng_p, ng_f, ng_v, ng_s
    double precision :: phi_old(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision :: phi_new(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    double precision ::   fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision ::   fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision ::   fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    double precision ::    velx(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    double precision ::    vely(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    double precision ::    velz(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    double precision ::  source(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    double precision :: dx, dt

    ! local variables
    integer i,j,k
    double precision :: dtdx

    dtdx = dt/dx

    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             phi_new(i,j,k) = phi_old(i,j,k) - &
                  ( velx(i+1,j,k)*fluxx(i+1,j,k)-velx(i,j,k)*fluxx(i,j,k) &
                   +vely(i,j+1,k)*fluxy(i,j+1,k)-vely(i,j,k)*fluxy(i,j,k) &
                   +velz(i,j,k+1)*fluxz(i,j,k+1)-velz(i,j,k)*fluxz(i,j,k) ) * dtdx
             phi_new(i,j,k) = phi_new(i,j,k) + dt * source(i,j,k)
          end do
       end do
    end do
    !$omp end parallel do

  end subroutine update_bds_phi_3d



end module update_phi_module
