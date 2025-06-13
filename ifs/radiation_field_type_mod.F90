! radiation_field_type_mod.F90 - FIELD API backed data structures
! that own the memory used by ecrads derived types
!
! (C) Copyright 2022- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! Author:  Johan Ericsson
! Email:   johan.ericsson@ecmwf.int
!


module radiation_field_type_module

  use parkind1,                 only : jprb, jpim ! Working precision, integer type

  ! FIELD API imports
  use field_module, only: field_2rb, field_3rb, field_2im
  use field_factory_module

  ! radiation imports
  use radiation_single_level,   only: single_level_type

  implicit none

  public

  type single_level_field_type

    integer :: nblocks, ncol, nalbedobands, nemissbands, n_bands_sw
    ! single level type access pointers
    real(jprb), pointer, dimension(:) :: &
         &   cos_sza=>null(), &          ! (ncol) Cosine of solar zenith angle
         &   skin_temperature=>null()    ! (ncol) Skin temperature (K)
    real(jprb), pointer, dimension(:,:) :: &
         &   sw_albedo=>null(), &        ! (ncol,nalbedobands)
         &   sw_albedo_direct=>null()    ! (ncol,nalbedobands)
    real(jprb), pointer, dimension(:,:) :: &
         &   lw_emissivity=>null()       ! (ncol,nemissbands) If
    real(jprb), pointer, dimension(:,:) :: &
         &   lw_emission=>null()         ! (ncol,nemissbands)
    real(jprb), pointer, dimension(:) :: &
         &   spectral_solar_scaling=>null() ! (n_bands_sw)
    integer, pointer, dimension(:) :: iseed=>null() ! (ncol)

    ! single level type fields
    class(field_2rb), pointer :: &
         &   f_cos_sza=>null(), &          ! (ncol,nblocks) Cosine of solar zenith angle
         &   f_skin_temperature=>null()    ! (ncol,nblocks) Skin temperature (K)
    class(field_3rb), pointer :: &
         &   f_sw_albedo=>null(), &        ! (ncol,nalbedobands,nblocks)
         &   f_sw_albedo_direct=>null()    ! (ncol,nalbedobands,nblocks)
    class(field_3rb), pointer :: &
         &   f_lw_emissivity=>null()       ! (ncol,nemissbands,nblocks) If
    class(field_3rb), pointer :: &
         &   f_lw_emission=>null()         ! (ncol,nemissbands,nblocks)
    class(field_2rb), pointer :: &
         &   f_spectral_solar_scaling=>null() ! (n_bands_sw, nblocks)
    class(field_2im), pointer :: f_iseed=>null() ! (ncol, nblocks)

    logical                   :: is_simple_surface = .true.
  contains

    procedure :: init => single_level_field_init
    procedure :: final => single_level_field_final
    procedure :: update_view => single_level_field_update_view
    procedure :: update_single_level => single_level_field_update_single_level

  end type single_level_field_type


contains


  !---------------------------------------------------------------------
  ! Initialise single_field_type
  subroutine single_level_field_init(this, nblocks, ncol, nalbedobands, nemisbands, &
       &                           use_sw_albedo_direct, is_simple_surface)

    use yomhook, only : lhook, dr_hook, jphook

    class(single_level_field_type), intent(inout) :: this
    integer,                  intent(in)    :: nblocks, ncol, nalbedobands, nemisbands
    logical,        optional, intent(in)    :: use_sw_albedo_direct
    logical,        optional, intent(in)    :: is_simple_surface

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type_module:single_level_field_init',0,hook_handle)

    if (present(is_simple_surface)) this%is_simple_surface = is_simple_surface

    call field_new(this%f_cos_sza, ubounds=[ncol, nblocks], persistent=.true.)

    if (this%is_simple_surface) then
      call field_new(this%f_skin_temperature, ubounds=[ncol, nblocks], persistent=.true.)
    else
      call field_new(this%f_lw_emission, ubounds=[ncol, nemisbands, nblocks], persistent=.true.)
    end if
    call field_new(this%f_lw_emissivity, ubounds=[ncol, nemisbands, nblocks], persistent=.true.)

    call field_new(this%f_sw_albedo, ubounds=[ncol, nalbedobands, nblocks], persistent=.true.)

    if (present(use_sw_albedo_direct)) then
      if (use_sw_albedo_direct) then
        call field_new(this%f_sw_albedo_direct, ubounds=[ncol, nalbedobands, nblocks], persistent=.true.)
      end if
    end if

    call field_new(this%f_iseed, ubounds=[ncol, nblocks], persistent=.true.)

    if (lhook) call dr_hook('radiation_field_type_module:single_level_field_init',1,hook_handle)
  end subroutine single_level_field_init


  !---------------------------------------------------------------------
  ! Finalise single_level_field_type
  subroutine single_level_field_final(this)

    use yomhook, only : lhook, dr_hook, jphook

    class(single_level_field_type), intent(inout) :: this

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:single_level_field_final',0,hook_handle)

    if (associated(this%f_cos_sza)) then
      call field_delete(this%f_cos_sza)
    end if
    this%cos_sza => null()

    if (associated(this%f_skin_temperature)) then
      call field_delete(this%f_skin_temperature)
    end if
    this%skin_temperature => null()

    if (associated(this%f_sw_albedo)) then
      call field_delete(this%f_sw_albedo)
    end if
    this%sw_albedo => null()

    if (associated(this%f_sw_albedo_direct)) then
      call field_delete(this%f_sw_albedo_direct)
    end if
    this%sw_albedo_direct => null()

    if (associated(this%f_lw_emissivity)) then
      call field_delete(this%f_lw_emissivity)
    end if
    this%lw_emissivity => null()

    if (associated(this%f_lw_emission)) then
      call field_delete(this%f_lw_emission)
    end if
    this%lw_emission => null()

    if (associated(this%f_spectral_solar_scaling)) then
      call field_delete(this%f_spectral_solar_scaling)
    end if
    this%spectral_solar_scaling => null()

    if (associated(this%f_iseed)) then
      call field_delete(this%f_iseed)
    end if
    this%iseed => null()

    if (lhook) call dr_hook('radiation_field_type:single_level_field_final',1,hook_handle)

  end subroutine single_level_field_final

  !---------------------------------------------------------------------
  ! Update view pointers with block level views
  subroutine single_level_field_update_view(this, block_index)

    use yomhook, only : lhook, dr_hook, jphook

    class(single_level_field_type), intent(inout) :: this
    integer,                        intent(in)    :: block_index
    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:single_level_field_update_view',0,hook_handle)

    if (associated(this%f_cos_sza)) this%cos_sza => this%f_cos_sza%get_view(block_index)

    if (associated(this%f_skin_temperature)) this%skin_temperature => this%f_skin_temperature%get_view(block_index)

    if (associated(this%f_sw_albedo)) this%sw_albedo => this%f_sw_albedo%get_view(block_index)

    if (associated(this%f_sw_albedo_direct)) this%sw_albedo_direct => this%f_sw_albedo_direct%get_view(block_index)

    if (associated(this%f_lw_emissivity)) this%lw_emissivity => this%f_lw_emissivity%get_view(block_index)

    if (associated(this%f_lw_emission)) this%lw_emission => this%f_lw_emission%get_view(block_index)

    if (associated(this%f_spectral_solar_scaling)) this%spectral_solar_scaling => this%f_spectral_solar_scaling%get_view(block_index)

    if (associated(this%f_iseed)) this%iseed => this%f_iseed%get_view(block_index)

    if (lhook) call dr_hook('radiation_field_type:single_level_field_update_view',1,hook_handle)

  end subroutine single_level_field_update_view

  !---------------------------------------------------------------------
  ! Update single_level_type with the view pointers of this object
  subroutine single_level_field_update_single_level(this, single_level)

    use yomhook, only : lhook, dr_hook, jphook

    class(single_level_field_type), intent(inout) :: this
    class(single_level_type),       intent(inout) :: single_level

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_field_type:field_update_single_level',0,hook_handle)

    if (associated(this%cos_sza)) single_level%cos_sza => this%cos_sza

    if (associated(this%skin_temperature)) single_level%skin_temperature => this%skin_temperature

    if (associated(this%sw_albedo)) single_level%sw_albedo => this%sw_albedo

    if (associated(this%sw_albedo_direct)) single_level%sw_albedo_direct => this%sw_albedo_direct

    if (associated(this%lw_emissivity)) single_level%lw_emissivity => this%lw_emissivity

    if (associated(this%lw_emission)) single_level%lw_emission => this%lw_emission

    if (associated(this%spectral_solar_scaling)) single_level%spectral_solar_scaling => this%spectral_solar_scaling

    single_level%is_simple_surface = this%is_simple_surface

    if (lhook) call dr_hook('radiation_field_type:field_update_single_level',1,hook_handle)

  end subroutine single_level_field_update_single_level

end module radiation_field_type_module
