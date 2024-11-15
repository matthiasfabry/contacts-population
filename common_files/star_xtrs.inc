   use star_lib
   use star_def
   use const_def
   use chem_def
   use num_lib
   use binary_def
   use binary_lib, only : binary_eval_rlobe
   use utils_lib, only : is_bad
   use interp_2d_lib_db
   use auto_diff

   implicit none


   real(dp), parameter :: expected_runtime = 120 ! minutes

   real(dp), pointer :: xvals(:), yvals(:), yvals_gtr_than_1(:), fpfunc1d(:), ftfunc1d(:), &
         irotfunc1d(:), otherrfunc1d(:), afunc1d(:)
   logical :: inter_ok = .false.
   integer :: num_xpts, num_ypts, num_ypts_gtr_than_1
   character(len=strlen) :: upstairs = '../../../../../../../data_tables/'  ! where fp/ft data lives

   ! logical params
   integer, parameter :: ilx_depleted_H = 1
   ! real params
   integer, parameter :: ix_ET_RL = 1, ix_ET_surf = 2, &
         ix_ET_RL_old = 3, ix_ET_surf_old = 4, &  ! s% xtra_old isnt saved in photos, so we do it manually
         ix_depth = 7, ix_wind = 8
   integer, parameter :: ilx_pre_ms = 1, ilx_had_ET = 7  ! indexes logical xtras of binary_info!

   real(dp), parameter :: nudge = 1d-4, &  ! dim-less step to be taken when evaluating numerical derivatives of the interpolating tables
         et_scale = 1d-2

contains

   subroutine extras_controls(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: xtest, ytest, testval
      type (auto_diff_real_1var_order1) :: test, A
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      ! this is the place to set any procedure pointers you want to change
      ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
      s% extras_startup => extras_startup
      s% extras_check_model => extras_check_model
      s% extras_start_step => extras_start_step
      s% extras_finish_step => extras_finish_step
      s% extras_after_evolve => extras_after_evolve

      s% how_many_extra_history_columns => how_many_extra_history_columns
      s% data_for_extra_history_columns => data_for_extra_history_columns
      s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
      s% data_for_extra_profile_header_items => data_for_extra_profile_header_items
      s% how_many_extra_profile_columns => how_many_extra_profile_columns
      s% data_for_extra_profile_columns => data_for_extra_profile_columns

      s% other_wind => brott_wind

      if(s% use_other_energy) then
         write(*, *) "use energy transfer!"
         s% other_energy => energy_transfer
      else
         write(*, *) "no energy transfer!"
      end if

      if(s% use_other_eval_fp_ft .and. s% use_other_eval_i_rot) then
         write(*, *) "use custom deformation!"
         s% other_eval_fp_ft => my_fp_ft
         s% other_eval_i_rot => my_irot
      else
         write(*, *) "use default deformation"
      end if

      if(s% use_other_surface_pt) then
         write(*, *) "use new BC!"
         s% other_surface_PT => my_atm
      else
         write(*, *) "use default BC!"
      end if

      s% job% warn_run_star_extras = .false.

      if (.not. inter_ok) then
         write(*, *) 'starting interpolator setup'
         call setup_interpolator(trim(upstairs) // 'fp_data.txt', xvals, num_xpts, yvals, &
               num_ypts, fpfunc1d, ierr)
         call setup_interpolator(trim(upstairs) // 'ft_data.txt', xvals, num_xpts, yvals, &
               num_ypts, ftfunc1d, ierr)
         call setup_interpolator(trim(upstairs) // 'irot_data.txt', xvals, num_xpts, &
               yvals, num_ypts, irotfunc1d, ierr)
         call setup_interpolator(trim(upstairs) // 'area_data.txt', xvals, num_xpts, &
               yvals, num_ypts, afunc1d, ierr)
         call setup_interpolator(trim(upstairs) // 'other_r_data.txt', xvals, num_xpts, &
               yvals, num_ypts, otherrfunc1d, ierr)
         xtest = -0.5
         ytest = 1.35
         ! test fp interpolator
         write(*, *) num_xpts, num_ypts, 'grid size'
         write(*, *) 'setup interpolators succesful,'

         call interp_evbipm_db(xtest, ytest, xvals, num_xpts, yvals, num_ypts,&
               fpfunc1d, num_xpts, testval, ierr)
         write(*, *) 'fp   test gave', testval, 'should be close to 0.6'
         call interp_evbipm_db(xtest, ytest, xvals, num_xpts, yvals, num_ypts,&
               ftfunc1d, num_xpts, testval, ierr)
         write(*, *) 'ft   test gave', testval, 'should be close to 0.8'
         call interp_evbipm_db(xtest, ytest, xvals, num_xpts, yvals, num_ypts,&
               irotfunc1d, num_xpts, testval, ierr)
         write(*, *) 'irot test gave', testval, 'should be close to 0.4'
         inter_ok = .true.
      end if

   end subroutine extras_controls

   subroutine extras_startup(id, restart, ierr)
      integer, intent(in) :: id
      logical, intent(in) :: restart
      integer, intent(out) :: ierr
      type (star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      if (.not. restart) then
         s% xtra(ix_ET_RL) = 0d0
         s% xtra(ix_ET_RL_old) = 0d0
         s% xtra(ix_ET_surf) = 0d0
         s% xtra(ix_ET_surf_old) = 0d0
         s% xtra(ix_depth) = 0d0
         s% xtra(ix_wind) = 0d0
         s% lxtra(ilx_depleted_H) = .false.
      end if

   end subroutine extras_startup

   integer function extras_start_step(id)  ! doesn't fire in binary runs!
      integer, intent(in) :: id

      extras_start_step = keep_going
   end function extras_start_step

   integer function extras_check_model(id)
      integer, intent(in) :: id
      type (star_info), pointer :: s
      integer :: ierr
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      extras_check_model = keep_going
   end function extras_check_model

   integer function extras_finish_step(id)
      integer, intent(in) :: id
      integer :: ierr, i
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)

      extras_finish_step = keep_going

      if (s% center_h1 < 1d-5 .and. .not. s% lxtra(ilx_depleted_H)) then
         write(*, *) "H depletion of star ", s% id
         s% lxtra(ilx_depleted_H) = .true.
         extras_finish_step = terminate
      end if

      s% xtra(ix_ET_RL_old) = s% xtra(ix_ET_RL)
      s% xtra(ix_ET_surf_old) = s% xtra(ix_ET_surf)

   end function extras_finish_step

   subroutine extras_after_evolve(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      real(dp) :: dt
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

   end subroutine extras_after_evolve

   integer function how_many_extra_history_columns(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_history_columns = 5

   end function how_many_extra_history_columns

   subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(n), thermal_eq_difference
      integer :: this_star, other_star, k, j
      integer, intent(out) :: ierr
      type (star_info), pointer :: s, so
      type (binary_info), pointer :: b
      ierr = 0

      call star_ptr(id, s, ierr)
      call binary_ptr(s% binary_id, b, ierr)

      call assign_stars(id, this_star, other_star, ierr)
      call star_ptr(b% star_ids(other_star), so, ierr)
      if (ierr /= 0) return

      names(1) = 'delta_shu_top'
      names(2) = 'delta_shu_mid'
      names(3) = "frac_ET_thickness"
      names(4) = 'thickness_bernouilli'
      names(5) = "thermal_eq_difference"
      vals(:) = 0d0

      if (s% xtra(ix_ET_RL) /= 0d0) then
         k = 1  ! find top of ET zone
         do while(s% r(k) > (1d0 + et_scale) * b% rl(this_star) .and. k < s% nz)
            k = k + 1
         end do
         vals(1) = abs(s% xtra(ix_ET_RL)) / (b% separation) ** 2 / &
               ((s% rho(k) * s% energy(k) + s% Peos(k)) * s% csound(k))

         do while(s% r(k) > b% rl(this_star))  ! find mid of ET zone
            k = k + 1
         end do
         vals(2) = abs(s% xtra(ix_ET_RL)) / (b% separation) ** 2 / &
               ((s% rho(k) * s% energy(k) + s% Peos(k)) * s% csound(k))
         vals(3) = b% separation / s% r(1) * pow(vals(1), 0.4d0)

         do while(s% r(k) > (1d0 - et_scale) * b% rl(this_star))  ! find bottom of ET zone
            k = k + 1
         end do

         j=1
         do while(so% r(j) > b% rl(other_star) .and. j < so% nz)
            j = j + 1
         end do
         vals(4) = pow(abs(s% xtra(ix_ET_RL)) / &
                 ((s% Cp(k) + so% Cp(j))*abs(s% T(k) - so% T(j))* s% rho(k)* &
                         pow(abs(s% Peos(k) - so% Peos(j))*(1/s% rho(k) + 1/so% rho(j)), 0.5d0)), 0.5d0) / &
                 b% separation
      end if

      thermal_eq_difference = 0d0
      do k = 1, s% nz
         if (s% q(k) > 0.99d0) cycle
         thermal_eq_difference = max(thermal_eq_difference, abs((s% dL_dm(k) - (s% eps_nuc(k) + s% extra_heat(k)% val))/(s% L(1)/s% m(1))))
      end do
      vals(5) = thermal_eq_difference

   end subroutine data_for_extra_history_columns

   integer function how_many_extra_profile_header_items(id)
      integer, intent(in) :: id
      integer :: ierr
      type (star_info), pointer :: s
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      how_many_extra_profile_header_items = 2
   end function how_many_extra_profile_header_items

   subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
      integer, intent(in) :: id, n
      character (len=maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(n)
      type(star_info), pointer :: s
      type(binary_info), pointer :: b
      integer, intent(out) :: ierr
      integer :: this_star, other_star
      ierr = 0
      call star_ptr(id,s,ierr)
      if(ierr/=0) return
      call binary_ptr(s% binary_id, b, ierr)

      call assign_stars(id, this_star, other_star, ierr)
      names(1) = 'rl'
      vals(1) = b% rl(this_star)
      names(2) = 'other_rl'
      vals(2) = b% rl(other_star)
   end subroutine data_for_extra_profile_header_items

   integer function how_many_extra_profile_columns(id)
      integer, intent(in) :: id
      how_many_extra_profile_columns = 5
   end function how_many_extra_profile_columns

   subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
      integer, intent(in) :: id, n, nz
      character (len = maxlen_profile_column_name) :: names(n)
      real(dp) :: vals(nz, n)
      real(dp) :: r_roche, a, m1, m2, lq
      real(dp) :: Area(nz)
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: this_star, other_star, i
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call binary_ptr(s% binary_id, b, ierr)

      call assign_stars(id, this_star, other_star, ierr)

      m1 = b% m(this_star)
      m2 = b% m(other_star)
      a = b% separation
      r_roche = b% rl(this_star)
      ! if masses or sep are not yet set at the binary level, use these
      if (m1 <= 0) m1 = b% m1 * Msun
      if (m2 <= 0) m2 = b% m2 * Msun
      if (a <= 0) a = pow(standard_cgrav * (m1 + m2) * &
            pow(b% initial_period_in_days * 86400, 2) / (4 * pi2), one_third)

      if (r_roche <= 0) r_roche = binary_eval_rlobe(m1, m2, a)
      lq = log10(m2 / m1)

      do i=1, s% nz
         Area(i) = eval_area(lq, s% r(i) / r_roche, ierr) * a * a
      end do

      names(1) = 'Teff_local'
      names(2) = 'r_tilde'
      names(3) = 'fp/ft'
      names(4) = 'dL_dm'
      names(5) = 'other_r_tilde'
      vals(:, 1) = pow(s% L(1:s% nz) / (boltz_sigma * Area(:)), 0.25d0)
      vals(:, 2) = (s% r(1:s% nz) - r_roche) / &
              (b% r(this_star) - r_roche)
      vals(:, 3) = s% fp_rot(1:s% nz)% val / s% ft_rot(1:s% nz)% val
      vals(:, 4) = s% dL_dm(1:s% nz)
      do i = 1, s% nz  ! other r tilde for this equipotential
         vals(i, 5) = (eval_other_r(lq, s% r(i) / r_roche, ierr) - 1) * b% rl(other_star) / &
                 (b% r(other_star) - b% rl(other_star))
      end do

   end subroutine data_for_extra_profile_columns

   ! *** custom routines ***
   ! interpolator data reading
   subroutine setup_interpolator(filename, &
         xs, num_xs, ys, num_ys, func1d, ierr)
      integer, intent(out) :: ierr, num_xs, num_ys
      character(LEN = *) :: filename
      integer :: k, iounit
      real(dp), pointer, intent(out) :: xs(:), ys(:), func1d(:)
      real(dp), pointer :: func(:,:,:)

      write(*, *) 'loading ', filename
      ! open data to interpolate
      open(newunit = iounit, file = filename, status = 'old', action = 'read',&
            iostat = ierr)

      read(iounit, *, iostat = ierr) num_xs
      allocate(xs(num_xs))
      do k = 1, num_xs
         read(iounit, *, iostat = ierr) xs(k)
      end do

      read(iounit, *, iostat = ierr) num_ys
      allocate(ys(num_ys))
      do k = 1, num_ys
         read(iounit, *, iostat = ierr) ys(k)
      end do

      ! create a 1d array with all the data, point func to it
      allocate(func1d(4 * num_xs * num_ys))
      func(1:4, 1:num_xs, 1:num_ys) => func1d(1:4 * num_xs * num_ys)
      do k = 1, num_xs
         read(iounit, *, iostat = ierr) func(1, k, :)
      end do

      if (ierr /= 0) then
         close(iounit)
      end if
      ! create interpolator
      call interp_mkbipm_db(xs, num_xs, ys, num_ys, func1d, num_xs, ierr)
   end subroutine setup_interpolator

   ! interpolator evaluators
   real(dp) function eval_other_r(lq, ar, ierr) result(other_r)
      ! evaluates, at same equipotential, r_other/r_rl_other as function of
      ! lq = log10(q == m(other)/m(this)) and ar = r_this/r_rl_this
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr
      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 1d-4
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts, &
         otherrfunc1d, num_xpts, other_r, ierr)
      if (ierr/=0) then
         write(*, *) "error in eval other r", ar, other_r
      end if
   end function eval_other_r

   real(dp) function eval_area(lq, ar, ierr) result(area)
      ! evaluates the area of the equipotential shell with equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), in units of separation^2
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr
      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 1d-4
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts, afunc1d,&
            num_xpts, area, ierr)
      if (ierr/=0) then
         write(*, *) "error in eval area", ar, area
      end if
   end function eval_area

   real(dp) function eval_fp(lq, ar, ierr) result(fp)
      ! evaluates fp of the equipotential shell with equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr
      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            fpfunc1d, num_xpts, fp, ierr)
      if (ierr/=0) then
         write(*, *) "error in eval fp", ar, fp
      end if
   end function eval_fp

   real(dp) function eval_ft(lq, ar, ierr) result(ft)
      ! evaluates ft of the equipotential shell with equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr
      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            ftfunc1d, num_xpts, ft, ierr)
      if (ierr/=0) then
         write(*, *) "error in eval ft", ar, ft
      end if
   end function eval_ft

   real(dp) function eval_irot(lq, ar, ierr) result(irot)
      ! evaluates moment of inertia irot of the equipotential shell with equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), in units of separation^2
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr
      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            irotfunc1d, num_xpts, irot, ierr)
      if (ierr/=0) then
         write(*, *) "error in eval irot", ar, irot
      end if
   end function eval_irot

   ! deformation
   subroutine my_fp_ft(id, nz, xm, r, rho, aw, ft, fp, r_polar, r_equatorial, report_ierr, ierr)
      integer, intent(in) :: id, nz
      real(dp), intent(in) :: aw(:), r(:), rho(:), xm(:) ! (nz)
      logical, intent(in) :: report_ierr
      real(dp), intent(inout) :: r_polar(:), r_equatorial(:) ! (nz)
      type (auto_diff_real_star_order1), intent(out) :: ft(:), fp(:)
      integer, intent(out) :: ierr

      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: j, this_star=0, other_star=0
      real(dp) :: r_roche, lq, m1, m2, a, ar

      include 'formats'

      if (.not. inter_ok) then
         write(*, *) "interpolators not setup, should happen in startup"
         stop
      end if

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr /= 0) return
      call assign_stars(id, this_star, other_star, ierr)
      if (ierr /= 0) return

      m1 = b% m(this_star)
      m2 = b% m(other_star)
      a = b% separation
      r_roche = b% rl(this_star)
      ! if masses or sep are not yet set at the binary level, use these
      if (m1 <= 0) m1 = b% m1 * Msun
      if (m2 <= 0) m2 = b% m2 * Msun
      if (a <= 0) a = pow(standard_cgrav * (m1 + m2) * &
            pow((b% initial_period_in_days) * 86400, 2) / (4 * pi2), one_third)
      if (r_roche <= 0) r_roche = binary_eval_rlobe(m1, m2, a)
      lq = log10(m2 / m1)

      !$OMP PARALLEL DO PRIVATE(j, ar) SCHEDULE(dynamic,2)
      do j = 1, s% nz  ! for every cell, compute fp, ft from an interpolating table
         ar = r(j) / r_roche
         ! set values
         fp(j) = eval_fp(lq, ar, ierr)
         ft(j) = eval_ft(lq, ar, ierr)
         ! set log derivatives
         fp(j)% d1Array(i_lnR_00) = (eval_fp(lq, ar + nudge, ierr) - fp(j)% val) / nudge * ar
         ft(j)% d1Array(i_lnR_00) = (eval_ft(lq, ar + nudge, ierr) - ft(j)% val) / nudge * ar
!            if (fp(j) == 0d0 .or. fp(j) == 0d0) write(*, *) j, ar, fp(j), ft(j), ierr  ! debug
         ! fix these to the current radius, they're only used for some wind mass loss enhancement
         r_polar(j) = r(j)
         r_equatorial(j) = r(j)
      end do
      !$OMP END PARALLEL DO
   end subroutine my_fp_ft

   subroutine my_irot(id, k, r00, w_div_w_crit_roche, i_rot)
      integer, intent(in) :: id, k
      real(dp), intent(in) :: r00, w_div_w_crit_roche
      type (auto_diff_real_star_order1), intent(out) :: i_rot
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: ierr, j, this_star=0, other_star=0
      real(dp) :: r_roche, a, m1, m2, ar, lq

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr/=0) return
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr/=0) return
      call assign_stars(id, this_star, other_star, ierr)
      if (ierr/=0) return

      m1 = b% m(this_star)
      m2 = b% m(other_star)
      a = b% separation
      r_roche = b% rl(this_star)
      if (m1 <= 0) m1 = b% m1 * Msun
      if (m2 <= 0) m2 = b% m2 * Msun
      if (a <= 0) a = pow(standard_cgrav * (m1 + m2) * &
         pow((b% initial_period_in_days) * 86400, 2) / (4 * pi2), one_third)
      if (r_roche <= 0) r_roche = binary_eval_rlobe(m1, m2, a)
!
      lq = log10(m2 / m1)
      i_rot = 0d0

      if (inter_ok) then
         ar = r00 / r_roche
         ! set value
         i_rot = eval_irot(lq, ar, ierr) * a * a
         ! set radius derivative
         i_rot% d1Array(i_lnR_00) = (eval_irot(lq, ar + nudge, ierr) * a * a - i_rot% val) / nudge * ar
!         write(*, *) r00, r_roche, i_rot, w_div_w_crit_roche  ! debug
      end if
   end subroutine my_irot

   ! energy transfer
   subroutine energy_transfer(id, ierr)
      integer, intent(in) :: id
      integer, intent(out) :: ierr
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: k, j, i, this_star, other_star
      real(dp) :: total_m, eps_RL, eps_surf

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return

      do k=1,s% nz
         s% extra_heat(k) = 0d0
      end do

      call binary_ptr(s% binary_id, b, ierr)

      if (b% lxtra(ilx_pre_ms) .or. s% model_number < 2) return
      if (s% xtra(ix_ET_RL) == 0d0 .or. s% doing_relax) return

      call assign_stars(id, this_star, other_star, ierr)

      ! once this triggers there has to be non zero energy input!
      if (.not. b% lxtra(ilx_had_ET)) b% lxtra(ilx_had_ET) = .true.
      call linear  ! or another profile of your design

      contains

         subroutine linear
            ! find depth
            i = 1
            do while(1 - s% m(i) / s% mstar < s% xtra(ix_depth))
               i = i + 1
            end do
            k = 1 ! measure mass chuck to dump energy in.
            do while(s% r(k) > (1d0 + et_scale)*s% r(i) .and. k < s% nz)
               k = k + 1
            end do
            j = k
            do while(s% r(j) > (1d0 - et_scale)*s% r(i))
               j = j + 1
            end do

            if (j==k) j = j+1  ! have at least one cell to put energy in

            total_m = s% m(k) - s% m(j)
            eps_RL = s% xtra(ix_ET_RL) / total_m

            do i=k, j-1
               s% extra_heat(i) = eps_RL
            end do

            ! surface component
            total_m = s% m(1) - s% m(j)
            eps_surf = s% xtra(ix_ET_surf) / total_m

            do i=1, j-1
               s% extra_heat(i) = s% extra_heat(i) + eps_surf
            end do
         end subroutine linear

   end subroutine energy_transfer

   ! surface BC
   subroutine my_atm(id, skip_partials,&
               lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
               lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)
      logical, intent(in) :: skip_partials
      integer, intent(in) :: id
      real(dp), intent(out) :: &
         lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
         lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
      integer, intent(out) :: ierr

      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      type (auto_diff_real_4var_order1) :: lnM1, L_surf, R_surf, lnR_surf, &
            lnT, lnP, lnk_surf, k_surf, m1, &
            A_surf, lq, r_roche, ar, fp_surf, ft_surf
      real(dp) :: m2, a
      integer :: this_star=0, other_star=0
      include 'formats'

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr/=0)  return
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr/=0) return

      call assign_stars(id, this_star, other_star, ierr)

      ! base autodiff variables
      L_surf = s% L(1)
      L_surf% d1val1 = 1d0

      R_surf = s% r(1)
      lnR_surf = log(R_surf)
      lnR_surf% d1val2 = 1d0

      m1 = b% m(this_star)
      if (m1% val <= 0) m1 = b% m1 * Msun
      lnM1 = log(m1)
      lnM1% d1val3 = 1d0

      k_surf = s% opacity(1)
      lnk_surf = log(k_surf)
      lnk_surf% d1val4 = 1d0

      ! other variables
      m2 = b% m(other_star)
      if (m2 <= 0) m2 = b% m2 * Msun
      a = b% separation
      r_roche = b% rl(this_star)
      if (r_roche% val <= 0) r_roche = binary_eval_rlobe(m1% val, m2, a)
      if (a <= 0) a = pow(standard_cgrav * (m1% val + m2) * &
         pow(b% initial_period_in_days * 86400, 2) / (4 * pi2), one_third)

      ! other autodiff variables
      lq = log10(m2 / m1)
      ar = R_surf / r_roche
      fp_surf = s% fp_rot(1)% val
      ft_surf = s% ft_rot(1)% val

      fp_surf% d1val2 = s% fp_rot(1)% d1Array(i_lnR_00)
      ft_surf% d1val2 = s% ft_rot(1)% d1Array(i_lnR_00)

      A_surf = eval_area(lq% val, ar% val, ierr) * a * a
      A_surf% d1val2 = (eval_area(lq% val, ar% val + nudge, ierr) * a * a - A_surf% val) / nudge * ar% val
      A_surf% d1val3 = (eval_area(lq% val + nudge, ar% val, ierr) * a * a - A_surf% val) / nudge

      ! calculate pressure and temperature
      lnP = log(two_thirds*(pi4*standard_cgrav*m1*fp_surf/(k_surf*A_surf*ft_surf)&
                  + L_surf/(clight*A_surf)))
      lnT = log(pow(L_surf / (boltz_sigma * A_surf), 0.25d0))

      lnT_surf = lnT% val
      dlnT_dL = lnT% d1val1
      dlnT_dlnR = lnT% d1val2
      dlnT_dlnM = lnT% d1val3
      dlnT_dlnkap = lnT% d1val4

      lnP_surf = lnP% val
      dlnP_dL = lnP% d1val1
      dlnP_dlnR = lnP% d1val2
      dlnP_dlnM = lnP% d1val3
      dlnP_dlnkap = lnP% d1val4

      if (isnan(lnT_surf)) then
         write(*, *) "negative temperature"
         ierr = -1
         return
      end if

   end subroutine my_atm

   ! wind
   subroutine brott_wind(id, Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z, w, ierr)
      use star_def
      integer, intent(in) :: id
      real(dp), intent(in) :: Lsurf, Msurf, Rsurf, Tsurf, X, Y, Z ! surface values (cgs)
      ! NOTE: surface is outermost cell. not necessarily at photosphere.
      ! NOTE: don't assume that vars are set at this point.
      ! so if you want values other than those given as args,
      ! you should use values from s% xh(:,:) and s% xa(:,:) only.
      ! rather than things like s% Teff or s% lnT(:) which have not been set yet.
      real(dp), intent(out) :: w ! wind in units of Msun/year (value is >= 0)
      integer, intent(out) :: ierr
      real(dp) :: Z_div_Z_solar, Teff_jump, alfa, L1, M1, R1, T1, &
            vink_wind, nieu_wind, hamann_wind, lowT_w, highT_w, Twindow
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr /= 0) return

      w = 0d0
      if (b% lxtra(ilx_pre_ms)) return

      L1 = Lsurf
      M1 = Msurf
      R1 = Rsurf
      T1 = Tsurf

      Z_div_Z_solar = s% kap_rq% Zbase / 0.0142d0
      ! use Vink et al 2001, eqns 14 and 15 to set "jump" temperature
      Teff_jump = 1d3 * (61.2d0 + 2.59d0 * (-13.636d0 + 0.889d0 * log10(Z_div_Z_solar)))

      vink_wind = 0d0
      nieu_wind = 0d0
      hamann_wind = 0d0

      call eval_Vink_wind(vink_wind)
      call eval_Nieuwenhuijzen_wind(nieu_wind)
      call eval_Hamann_wind(hamann_wind)

      ! use 1/10 hamann
      hamann_wind = hamann_wind / 10d0

      lowT_w = max(vink_wind, nieu_wind)

      alfa = 0d0
      if (X >= 0.7d0) then
         alfa = 1d0
      else if (X > 0.4d0 .and. X < 0.7d0) then
         alfa = (X - 0.4d0) / 0.3d0
      end if
      highT_w = alfa * vink_wind + (1d0 - alfa) * hamann_wind

      ! have a 10% Teff_jump window to switch from the lowT to the highT wind
      Twindow = Teff_jump * 0.10d0
      alfa = 0d0
      if (T1 <= Teff_jump - Twindow / 2d0) then
         alfa = 1d0
      else if (T1 > Teff_jump - Twindow / 2d0 .and. T1 < Teff_jump + Twindow / 2d0) then
         alfa = ((Teff_jump + Twindow / 2d0) - T1) / Twindow
      end if
      w = alfa * lowT_w + (1d0 - alfa) * highT_w

      ! further soften change in wind to avoid things going bad
      if (s% xtra(ix_wind) /= 0) then
         if(abs(w) > abs(s% xtra(ix_wind)) * 1.05) then
            w = s% xtra(ix_wind) * 1.05
         else if(abs(w) < abs(s% xtra(ix_wind)) * 0.95) then
            w = s% xtra(ix_wind) * 0.95
         end if
      end if

      s% xtra(ix_wind) = w

      ierr = 0

   contains

      subroutine eval_Vink_wind(w)
         real(dp), intent(inout) :: w
         real(dp) :: alfa, w1, w2, logMdot, dT, vinf_div_vesc

         ! alfa = 1 for hot side, = 0 for cool side
         if (T1 > 27500d0) then
            alfa = 1
         else if (T1 < 22500d0) then
            alfa = 0
         else
            dT = 100d0
            if (T1 > Teff_jump + dT) then
               alfa = 1
            else if (T1 < Teff_jump - dT) then
               alfa = 0
            else
               alfa = (T1 - (Teff_jump - dT)) / (2 * dT)
            end if
         end if

         if (alfa > 0) then ! eval hot side wind (eqn 24)
            vinf_div_vesc = 2.6d0 ! this is the hot side galactic value
            vinf_div_vesc = vinf_div_vesc * pow(Z_div_Z_solar, 0.13d0) ! corrected for Z
            logMdot = &
                  - 6.697d0 &
                        + 2.194d0 * log10(L1 / Lsun / 1d5) &
                        - 1.313d0 * log10(M1 / Msun / 30) &
                        - 1.226d0 * log10(vinf_div_vesc / 2d0) &
                        + 0.933d0 * log10(T1 / 4d4) &
                        - 10.92d0 * pow2(log10(T1 / 4d4)) &
                        + 0.85d0 * log10(Z_div_Z_solar)
            w1 = 10**(logMdot)
         else
            w1 = 0
         end if

         if (alfa < 1) then ! eval cool side wind (eqn 25)
            vinf_div_vesc = 1.3d0 ! this is the cool side galactic value
            vinf_div_vesc = vinf_div_vesc * pow(Z_div_Z_solar, 0.13d0) ! corrected for Z
            logMdot = &
                  - 6.688d0 &
                        + 2.210d0 * log10(L1 / Lsun / 1d5) &
                        - 1.339d0 * log10(M1 / Msun / 30) &
                        - 1.601d0 * log10(vinf_div_vesc / 2d0) &
                        + 1.07d0 * log10(T1 / 2d4) &
                        + 0.85d0 * log10(Z_div_Z_solar)
            w2 = 10**(logMdot)
         else
            w2 = 0
         end if

         w = alfa * w1 + (1 - alfa) * w2

      end subroutine eval_Vink_wind

      subroutine eval_Nieuwenhuijzen_wind(w)
         ! Nieuwenhuijzen, H.; de Jager, C. 1990, A&A, 231, 134 (eqn 2)
         real(dp), intent(out) :: w
         real(dp) :: log10w
         include 'formats'
         log10w = -14.02d0 &
               + 1.24d0 * log10(L1 / Lsun) &
               + 0.16d0 * log10(M1 / Msun) &
               + 0.81d0 * log10(R1 / Rsun) &
               + 0.85d0 * log10(Z_div_Z_solar)
         w = exp10(log10w)
      end subroutine eval_Nieuwenhuijzen_wind

      subroutine eval_Hamann_wind(w)
         ! Hamann, W.-R.; Koesterke, L.; Wessolowski, U. 1995, A&A, 299, 151
         real(dp), intent(out) :: w
         real(dp) :: log10w
         include 'formats'
         log10w = -11.95d0 &
               + 1.5d0 * log10(L1 / Lsun) &
               - 2.85d0 * X &
               + 0.85d0 * log10(Z_div_Z_solar)
         w = exp10(log10w)
      end subroutine eval_Hamann_wind

   end subroutine brott_wind

   ! helper
   subroutine assign_stars(id, this_star, other_star, ierr)
      ! determine which star is which in the binary,
      ! this module works at the star level after all
      integer, intent(in) :: id
      integer, intent(out) :: this_star, other_star, ierr
      type (binary_info), pointer :: b
      type (star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      call binary_ptr(s% binary_id, b, ierr)

      if (ierr /= 0) return

      if (b% s1% id == s% id) then
         this_star = 1
         other_star = 2
      else if (b% s2% id == s% id) then
         this_star = 2
         other_star = 1
      else
         ierr = 1
         return
      end if

   end subroutine assign_stars