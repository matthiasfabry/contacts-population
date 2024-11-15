   use star_lib
   use star_def
   use const_def
   use chem_def
   use num_lib
   use binary_def
   use binary_lib
   use run_star_extras

   implicit none

   ! logical params
   !!! ilx_pre_ms = 1, ilx_had_ET = 7 defined in run_star_extras !!!
   integer, parameter :: ilx_zams_l2_overflow = 2, ilx_doing_period_relaxation = 3, &
         ilx_had_rlof_at_zams = 4, ilx_had_contact = 5, &
         ilx_had_caseA = 8, ilx_had_max_mdot = 9
   ! integer params
   integer, parameter :: iix_time0 = 1, iix_time1 = 2, iix_clockrate = 3
   ! real params
   integer, parameter :: ix_running_time = 1

   ! inlist settings we wanna remember
   logical :: inlist_L2_setting
   real(dp) :: inlist_fr, inlist_implicit_tolerance, inlist_fr_limit, inlist_fr_hard,&
          inlist_scale_max_correction, inlist_delta_HR, inlist_delta_lgTeff, inlist_delta_lgL

   ! internal variables
   integer :: slow_models = 0, consec_slow = 0
   logical :: too_slow = .false.
   real(dp) :: pre_ms_time, zams_age, deltaJ, deltaT, func_eval
   real(dp) :: et_smoothing = 0.5  ! 0 = no smoothing, -> 1 = extreme smoothing

contains

   ! extras_ functions
   subroutine extras_binary_controls(binary_id, ierr)
      integer :: binary_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      ! Set these function pointers to point to the functions you wish to use in
      ! your run_binary_extras. Any which are not set, default to a null_ version
      ! which does nothing.
      b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
      b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
      b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
      b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

      b% extras_binary_startup => extras_binary_startup
      b% extras_binary_start_step => extras_binary_start_step
      b% extras_binary_check_model => extras_binary_check_model
      b% extras_binary_finish_step => extras_binary_finish_step
      b% extras_binary_after_evolve => extras_binary_after_evolve

      ! period relaxation if overflow L1 at ZAMS
      b% other_extra_jdot => my_extra_jdot
      if (b% use_other_implicit_function_to_solve) then
         b% other_implicit_function_to_solve => my_contact_func_to_solve
      end if

      ! Once you have set the function pointers you want, then uncomment this (or set it in your
      ! star_job inlist)
      ! to disable the printed warning message,
      b% warn_binary_extra = .false.

   end subroutine extras_binary_controls

   integer function extras_binary_startup(binary_id, restart, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      logical, intent(in) :: restart
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if

      inlist_l2_setting = b% terminate_if_L2_overflow
      inlist_fr = b% fr
      inlist_implicit_tolerance = b% implicit_scheme_tolerance
      inlist_fr_limit = b% fr_limit
      inlist_scale_max_correction = b% s1% scale_max_correction
      inlist_delta_HR = b% s1% delta_HR_limit
      inlist_delta_lgTeff = b% s1% delta_lgTeff_limit
      inlist_delta_lgL = b% s1% delta_lgL_limit
      inlist_fr_hard = b% fr_hard
      call system_clock(b% ixtra(iix_time0), b% ixtra(iix_clockrate))
      ! also set old in case retry on a restart
      call system_clock(b% ixtra_old(iix_time0), b% ixtra_old(iix_clockrate))

      if (.not. restart) then
         b% terminate_if_L2_overflow = .false.  ! prevent l2 stoppage on PreMS
         b% lxtra(:) = .false.
         b% lxtra(ilx_pre_ms) = .true.
         b% xtra(:) = 0d0
         b% ixtra(:) = 0d0
      end if

      pre_ms_time = pow(b% m1, -2.3d0) * 1e7  ! rough estimate of the PreMS lifetime

      extras_binary_startup = keep_going
   end function  extras_binary_startup

   integer function extras_binary_start_step(binary_id, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      real(dp) :: new_limit
      call binary_ptr(binary_id, b, ierr)

      extras_binary_start_step = keep_going
      if (ierr /= 0) return

      ! if pre ms, do only this and continue
      if (b% lxtra(ilx_pre_ms)) then
         b% fr = -1d0
         return
      else
         b% fr = inlist_fr
         b% fr_hard = inlist_fr_hard
      end if

      !!! if here, binary must be past pre-ms !!!

!      ! smooth composition of outer 0.01 Msuns of accretor
      b% s_accretor% smooth_outer_xa_big = 0.02
      b% s_accretor% smooth_outer_xa_small = 0.01
      b% s_donor% smooth_outer_xa_big = -1d0
      b% s_donor% smooth_outer_xa_small = -1d0

      ! be more careful when ET is working
      if (b% r(1) > b% rl(1) .and. b% r(2) > b% rl(2) .and. b% s1% use_other_energy) then
         b% fr = 5d-1 * b% fr
      end if

      if (b% mtransfer_rate /= 0d0 .and. &
            max(abs(b% mdot_system_wind(b% d_i)/b% mtransfer_rate), &
                abs(b% mdot_system_wind(b% a_i)/b% mtransfer_rate)) > 5d-1) then
         ! if here, the mdot scheme will struggle to find a solution, either by
         ! looking to keep the donor exactly at the RL, or because the accretor
         ! has a total mdot of around zero. Let's ignore fr so we smooth over
         ! the phase where mdot_trans \approx mdot_wind, we always respect the
         ! contact condition, so let's hope things don't go bad. b% fm will
         ! kick in I think.
         b% fr = -1d0
         b% fr_hard = -1d0
         write(*, *) "high wind/mtransfer ratio, ignoring fr", &
               max(abs(b% mdot_system_wind(b% d_i)/b% mtransfer_rate), &
                abs(b% mdot_system_wind(b% a_i)/b% mtransfer_rate))
      end if

      call set_hr_limits(b% s1)
      call set_hr_limits(b% s2)


   end function extras_binary_start_step

   subroutine set_hr_limits(s)
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: ierr
      logical :: contact

      call binary_ptr(s% binary_id, b, ierr)
      contact = (b% r(1) >= b% rl(1) .and. b% r(2) >= b% rl(2))

      if (contact .and. s% use_other_energy) then  ! be more careful as ET turns on.
         s% delta_HR_limit = 1d-1 * inlist_delta_HR
         s% delta_lgTeff_limit = 1d-1 * inlist_delta_lgTeff
         s% delta_lgL_limit = 1d-1 * inlist_delta_lgL
      else
         s% delta_HR_limit = inlist_delta_HR
         s% delta_lgTeff_limit = inlist_delta_lgTeff
         s% delta_lgL_limit = inlist_delta_lgL
      end if

   end subroutine set_hr_limits

   integer function extras_binary_check_model(binary_id)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer :: ierr
      logical :: use_sum, detached, explicit
      real(dp) :: explicit_mdot, new_limit

      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then ! failure in  binary_ptr
         return
      end if
      extras_binary_check_model = keep_going

      if (b% lxtra(ilx_doing_period_relaxation) .and. &
            b% period < b% initial_period_in_days * secday) then
         extras_binary_check_model = retry  ! retry when we overshoot our period shrinking
         write(*, *) "overshot period shrinkage, retrying"
      end if

   end function extras_binary_check_model

   integer function extras_binary_finish_step(binary_id)
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      type (binary_info), pointer :: b
      type (star_info), pointer :: s1, s2
      integer, intent(in) :: binary_id
      integer :: ierr, id, i
      logical :: contact
      real(dp) :: r2, lq, q, sig, dummy_mdot
      character(10) :: time

      extras_binary_finish_step = keep_going
      call binary_ptr(binary_id, b, ierr)

      if (b% s1% model_number > 150 .and. b% lxtra(ilx_pre_ms)) then
         extras_binary_finish_step = terminate
         write(*, *) "Terminate due to pre-MS evolution taking too long"
         return
      end if

      ! check if models are progressing with decent pace
      if (log(b% s1% dt / secyer) < -5) then
         slow_models = slow_models + 1
         consec_slow = consec_slow + 1
      else
         consec_slow = 0
      end if
      if (slow_models > 500 .or. consec_slow > 100) then  ! let's not bother evolving a snail
         extras_binary_finish_step = terminate
         write(*, *) 'models are evolving too slowly, stopping'
         too_slow = .true.
      end if

      contact = b% rl_relative_gap(1) > 0d0 .and. b% rl_relative_gap(2) > 0d0

      ! Pre main sequence check
      if (b% lxtra(ilx_pre_ms) .and. &
            abs(log10(abs(b% s1% L_nuc_burn_total * Lsun / b% s1% L(1)))) < 0.005 .and. &
            abs(log10(abs(b% s2% L_nuc_burn_total * Lsun / b% s2% L(1)))) < 0.005 .and. &
            b% s1% star_age > pre_ms_time) then
         ! if here, both stars have thermal equilibrium (reached ZAMS), so activate RLOF
         if (b% m(1) > b% m(2)) then  ! make q always > 1 and select smaller star
            q = b% m(1) / b% m(2)
            id = 2
         else
            q = b% m(2) / b% m(1)
            id = 1
         end if
         sig = 62.9237d0 / (15.9839d0 + pow(q, 0.2240d0))
         if (b% rl_relative_gap(id) > 3.3752d0 / &
               (1 + ((log(q) + 1.0105d0) / sig) ** 2) &
               / (9.0087d0 + pow(q, -0.4022d0))) then
            write(*, *) "model is overflowing L2 at ZAMS"
            b% lxtra(ilx_zams_l2_overflow) = .true.
            extras_binary_finish_step = terminate
            return
         else
            write(*, *) "model is not overflowing L2 at ZAMS"
         end if
         b% lxtra(ilx_pre_ms) = .false.
         b% ignore_rlof_flag = .false.
         b% terminate_if_L2_overflow = inlist_L2_setting
         write(*, *) "Engage RLOF!"
         if (b% rl_relative_gap(b% d_i) > 0) then
            write(*, *) "overflowing L1 at ZAMS! Doing period relaxation procedure!"
            b% lxtra(ilx_had_rlof_at_zams) = .true.
            b% lxtra(ilx_doing_period_relaxation) = .true.
            deltaJ = b% angular_momentum_J
            ! boost period to put binary out of contact
            call binary_set_period_eccentricity(binary_id, &
                  b% initial_period_in_days * secday &
                  * (1.3d0 * max(b% rl_relative_gap(1) + 1, b% rl_relative_gap(2) + 1)), 0d0, ierr)
            ! keep stars synchronized
            if (b% point_mass_i /= 1 .and. b% s1% rotation_flag) then
               call star_set_uniform_omega(b% s1% id, 2 * pi / b% period, ierr)
            end if
            if (b% point_mass_i /= 2 .and. b% s2% rotation_flag) then
               call star_set_uniform_omega(b% s2% id, 2 * pi / b% period, ierr)
            end if
            deltaJ = b% angular_momentum_j - deltaJ ! remove added J.
            deltaT = b% s1% kh_timescale * secyer  ! remove it on thermal timescale
         end if
      else if (b% lxtra(ilx_pre_ms) .and. &
            (abs(log10(abs(b% s1% L_nuc_burn_total * Lsun / b% s1% L(1)))) > 0.005 .or. &
                  abs(log10(abs(b% s2% L_nuc_burn_total * Lsun / b% s2% L(1)))) > 0.005 .or. &
                  b% s1% star_age < pre_ms_time)) then
         write(*, *) "still not at ZAMS, keep period fixed"
         call binary_set_period_eccentricity(binary_id, &
               b% initial_period_in_days * secday, 0d0, ierr)
         ! keep stars synchronized
         if (b% point_mass_i /= 1 .and. b% s1% rotation_flag) then
            call star_set_uniform_omega(b% s1% id, 2 * pi / b% period, ierr)
         end if
         if (b% point_mass_i /= 2 .and. b% s2% rotation_flag) then
            call star_set_uniform_omega(b% s2% id, 2 * pi / b% period, ierr)
         end if
      end if

      ! check for completion of period relaxation
      if (b% lxtra(ilx_doing_period_relaxation)) then
         if (b% period <= 1.01 * b% initial_period_in_days * secday) then
            b% lxtra(ilx_doing_period_relaxation) = .false.
            write(*, *) "period relaxation complete! resuming without other_jdot"
         else
            write(*, *) "still doing period relaxation"
            ! keep stars synchronized
            if (b% point_mass_i /= 1 .and. b% s1% rotation_flag) then
               call star_set_uniform_omega(b% s1% id, 2 * pi / b% period, ierr)
            end if
            if (b% point_mass_i /= 2 .and. b% s2% rotation_flag) then
               call star_set_uniform_omega(b% s2% id, 2 * pi / b% period, ierr)
            end if
         end if
      end if

      ! don't do anything below if we're pre-MS
      if (b% lxtra(ilx_pre_ms)) return

      !!! do post-step checks !!!
      ! mark if case A MT started
      if (.not. b% lxtra(ilx_had_caseA) .and. abs(b% mtransfer_rate) >= 1d-10) then
         b% lxtra(ilx_had_caseA) = .true.
      end if

      ! mark if main sequence contact phase occurred
      if (.not. b% lxtra(ilx_had_contact) .and. contact) then
         b% lxtra(ilx_had_contact) = .true.
      end if

      ! terminate if no case A occurred
      if ((b% s1% lxtra(ilx_depleted_H) .or. b% s2% lxtra(ilx_depleted_H)) .and. &
            .not. b% lxtra(ilx_had_caseA)) then
         extras_binary_finish_step = terminate
         write(*, *) "no case A interaction!"
      end if

      ! terminate if max_mdot (shouldn't trigger? let's put it anyway to be sure...)
      if (b% mtransfer_rate <= -b% max_implicit_abs_mdot*Msun/secyer) then
         extras_binary_finish_step = terminate
         b% lxtra(ilx_had_max_mdot) = .true.
         write(*, *) "Terminate because of reaching max mdot_trans"
      end if

      ! set energy transfer stuff for coming step
      if (contact .and. b% s1% use_other_energy) then
         i = 1
         do while (b% s1% r(i) >= b% rl(1))
            i = i + 1
         end do
         b% s1% xtra(ix_depth) = 1 - b% s1% m(i) / b% s1% mstar
         !         write(*, *) " primary i, mstar", i, b% s1% xtra(ix_depth), b% s1% mstar
         i = 1
         do while (b% s2% r(i) >= b% rl(2))
            i = i + 1
         end do
         b% s2% xtra(ix_depth) = 1 - b% s2% m(i) / b% s2% mstar
         !         write(*, *) 'secondary i', i, b% s2% xtra(ix_depth)
         call set_luminosity_transfers(binary_id, ierr)
      else  ! no ET coming step
         b% s1% xtra(ix_ET_RL) = 0d0
         b% s2% xtra(ix_ET_RL) = 0d0
         b% s1% xtra(ix_ET_surf) = 0d0
         b% s2% xtra(ix_ET_surf) = 0d0
         b% s1% xtra(ix_ET_RL_old) = 0d0
         b% s2% xtra(ix_ET_RL_old) = 0d0
         b% s1% xtra(ix_ET_surf_old) = 0d0
         b% s2% xtra(ix_ET_surf_old) = 0d0
         b% s1% xtra(ix_depth) = 0d0
         b% s2% xtra(ix_depth) = 0d0
      end if

   end function extras_binary_finish_step

   subroutine extras_binary_after_evolve(binary_id, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr

      call binary_ptr(binary_id, b, ierr)

      if (b% s1% use_other_energy .and. .not. b% lxtra(ilx_had_ET)) then
         open(99, file = 'had_no_ET.txt')  ! signal with file ET was not triggered
         close(99)
      end if

   end subroutine extras_binary_after_evolve

   ! functions for extra data
   integer function how_many_extra_binary_history_header_items(binary_id)
      use binary_def, only: binary_info
      integer, intent(in) :: binary_id
      how_many_extra_binary_history_header_items = 0
   end function how_many_extra_binary_history_header_items

   subroutine data_for_extra_binary_history_header_items(binary_id, n, names, vals, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id, n
      character (len = maxlen_binary_history_column_name) :: names(n)
      real(dp) :: vals(n)
      integer, intent(out) :: ierr
   end subroutine data_for_extra_binary_history_header_items

   integer function how_many_extra_binary_history_columns(binary_id)
      integer, intent(in) :: binary_id
      how_many_extra_binary_history_columns = 13
   end function how_many_extra_binary_history_columns

   subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
      type (binary_info), pointer :: b
      integer, intent(in) :: binary_id, n
      character (len = maxlen_binary_history_column_name) :: names(n)
      real(dp) :: vals(n), dt
      integer, intent(out) :: ierr
      integer :: count_max

      call binary_ptr(binary_id, b, ierr)
      if(ierr /= 0) return

      ! if max mdot, mark contact manually since termination called by binary_mdot
      if (abs(b% mtransfer_rate) >= b% max_implicit_abs_mdot * Msun / secyer) then
         b% lxtra(ilx_had_max_mdot) = .true.
         b% lxtra(ilx_had_contact) = .true.
      end if

      ! poll the clock for this step
      call system_clock(b% ixtra(iix_time1), b% ixtra(iix_clockrate), count_max)

      dt = dble(b% ixtra(iix_time1) - b% ixtra(iix_time0)) / b% ixtra(iix_clockrate) / 60
      if (dt < 0) then  ! clock folded back to zero
         dt = dt + dble(count_max) / b% ixtra(iix_clockrate) / 60
      end if
      b% xtra(ix_running_time) = b% xtra(ix_running_time) + dt
      b% ixtra(iix_time0) = b% ixtra(iix_time1)

      ! write extra history
      names(1) = 'ET_RL_1'
      vals(1) = b% s1% xtra(ix_ET_RL) / Lsun
      names(2) = 'ET_RL_2'
      vals(2) = b% s2% xtra(ix_ET_RL) / Lsun
      names(3) = 'ET_surf_1'
      vals(3) = b% s1% xtra(ix_ET_surf) / Lsun
      names(4) = 'ET_surf_2'
      vals(4) = b% s2% xtra(ix_ET_surf) / Lsun
      call contact_condition(b, func_eval)
      names(5) = 'contact_condition'
      vals(5) = func_eval
      names(6) = 'stopping_condition'
      if (b% s1% lxtra(ilx_depleted_H) .or. b% s2% lxtra(ilx_depleted_H)) then
         vals(6) = 1  ! 'survival condition'
      else if (b% lxtra(ilx_zams_l2_overflow)) then
         vals(6) = 2  ! l2 at zams
      else if (b% s_donor% termination_code == t_xtra1 .and. &
            termination_code_str(t_xtra1) == "Terminate because of L2 overflow") then
         vals(6) = 3  ! l2 overflow
      else if (.not. b% lxtra(ilx_had_caseA) .and. b% s1% lxtra(ilx_depleted_H)) then
         vals(6) = 4  ! no case A
      else if (b% lxtra(ilx_had_max_mdot)) then
         vals(6) = 5  ! max mdot
      else if (b% s1% termination_code == t_min_timestep_limit .or. &
            b% s2% termination_code == t_min_timestep_limit .or. &
            b% s1% termination_code == t_max_number_retries .or. &
            b% s2% termination_code == t_max_number_retries .or. &
            too_slow) then
         vals(6) = -99  ! problems (no bueno)
      else  ! still running
         vals(6) = -1
      end if
      names(7) = 'mtransfer_timescale'
      if (b% mtransfer_rate /= 0d0) then
         vals(7) = b% m(b% d_i) / b% mtransfer_rate
      else
         vals(7) = 0d0
      end if
      names(8) = 'q'
      vals(8) = b% m(2) / b% m(1)
      names(9) = 'pre_zams_or_relax'
      vals(9) = bool_number(b% lxtra(ilx_pre_ms) .or. b% lxtra(ilx_doing_period_relaxation))
      names(10) = 'runtime_minutes'
      vals(10) = b% xtra(ix_running_time)
      names(11) = 'had_RLOF_at_ZAMS'
      vals(11) = bool_number(b% lxtra(ilx_had_rlof_at_zams))
      names(12) = 'had_ET'
      vals(12) = bool_number(b% lxtra(ilx_had_ET))
      names(13) = 'had_contact'
      vals(13) = bool_number(b% lxtra(ilx_had_contact))

   end subroutine data_for_extra_binary_history_columns

   ! make float from a boolean
   real(dp) function bool_number(bool)
      logical, intent(in) :: bool
      if (bool) then
         bool_number = 1d0
      else
         bool_number = 0d0
      end if
   end function bool_number

   ! custom physics functions
   ! ET
   subroutine set_luminosity_transfers(binary_id, ierr)
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      integer :: roche_cell1, roche_cell2
      type (binary_info), pointer :: b
      real(dp) :: new_l_trans_RL_1, new_l_trans_RL_2, new_l_trans_surf_1, &
            new_l_trans_surf_2, p, scale

      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr/=0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      if(b% rl_relative_gap(1) >= et_scale .and. b% rl_relative_gap(2) >= et_scale) then
         scale = 1d0
      else if(b% rl_relative_gap(1) >= 0.0 .and. b% rl_relative_gap(2) >= 0.0) then
         scale = min(b% rl_relative_gap(1), b% rl_relative_gap(2)) / et_scale
      else
         scale = 0d0
      end if

      ! if timestep is low, smooth the ET a lot
      if (log10(b% s1% dt_next/secyer) < 1) then
         et_smoothing = 0.95
      else if (log10(b% s1% dt_next/secyer) < 3) then
         et_smoothing = 0.95 - (0.95 - 0.5) / 2 * (log10(b% s1% dt_next/secyer) - 1)
      else
         et_smoothing = 0.5
      end if


      p = 1 - et_smoothing
      call calculate_transfer(binary_id, b% rl(1), b% rl(2), &
            new_l_trans_RL_1, new_l_trans_RL_2, ierr)
      b% s1% xtra(ix_ET_RL) = &
            (new_l_trans_RL_1 * p + b% s1% xtra(ix_ET_RL_old) * (1 - p)) * scale
      b% s2% xtra(ix_ET_RL) = &
            (new_l_trans_RL_2 * p + b% s2% xtra(ix_ET_RL_old) * (1 - p)) * scale

      call calculate_transfer(binary_id, b% r(1), b% r(2), &
            new_l_trans_surf_1, new_l_trans_surf_2, ierr)
      b% s1% xtra(ix_ET_surf) = &
            (new_l_trans_surf_1 * p + b% s1% xtra(ix_ET_surf_old) * (1 - p)) * scale
      b% s2% xtra(ix_ET_surf) = &
            (new_l_trans_surf_2 * p + b% s2% xtra(ix_ET_surf_old) * (1 - p)) * scale
      write(*, *) "et smoothing", et_smoothing
      write(*, *) "RL1, new, old", new_l_trans_RL_1 / Lsun, b% s1% xtra(ix_ET_RL_old) / Lsun
      write(*, *) "ET calculations of s1, at RL, at surf", &
            b% s1% xtra(ix_ET_RL) / Lsun, b% s1% xtra(ix_ET_surf) / Lsun

   end subroutine set_luminosity_transfers

   subroutine calculate_transfer(binary_id, loc_1, loc_2, trans_1, trans_2, ierr)
      integer, intent(in) :: binary_id
      real(dp), intent(in) :: loc_1, loc_2
      real(dp), intent(out) :: trans_1, trans_2
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      real(dp) :: lq, l1, l2, s_g1, s_g2, new_l1, new_l2
      integer :: cell_1, cell_2

      ierr = 0
      trans_1 = 0d0
      trans_2 = 0d0

      call binary_ptr(binary_id, b, ierr)
      if (ierr/=0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if
      lq = log10(b% m(2) / b% m(1))

      ! get luminosities at the relevant cells
      cell_1 = 1
      do while (b% s1% r(cell_1) > loc_1 .and. cell_1 < b% s1% nz)
         cell_1 = cell_1 + 1
      end do
      if (cell_1 >= b% s1% nz) return
      l1 = b% s1% L(cell_1)

      cell_2 = 1
      do while (b% s2% r(cell_2) > loc_2 .and. cell_2 < b% s2% nz)
         cell_2 = cell_2 + 1
      end do
      if (cell_2 >= b% s2% nz) return
      l2 = b% s2% L(cell_2)
      ! compute area times average gravity via fp/ft (see Fabry+2022)
      s_g1 = pi4 * standard_cgrav * b% s1% m(cell_1) * &
            eval_fp(lq, b% s1% r(cell_1) / b% rl(1), ierr) / &
            eval_ft(lq, b% s1% r(cell_1) / b% rl(1), ierr)
      s_g2 = pi4 * standard_cgrav * b% s2% m(cell_2) * &
            eval_fp(-lq, b% s2% r(cell_2) / b% rl(2), ierr) / &
            eval_ft(-lq, b% s2% r(cell_2) / b% rl(2), ierr)
      ! redistribute radiative luminosity proportional to area times gravity
      new_l1 = (l1 + l2) * s_g1 / (s_g1 + s_g2)
      new_l2 = (l1 + l2) * s_g2 / (s_g1 + s_g2)
      trans_1 = new_l1 - l1
      trans_2 = new_l2 - l2
   end subroutine calculate_transfer

   ! MT
   subroutine contact_condition(b, func)
      type (binary_info), pointer :: b
      integer :: ierr
      real(dp), intent(out) :: func
      real(dp) :: q, slope_contact

      ierr = 0
      q = b% m(b% a_i) / b% m(b% d_i)
      slope_contact = pow(q, -0.52d0)
      ! If accretor is overflowing its Roche lobe, then the contact scheme needs to be used.
      ! Otherwise, if accretor radius is (within tolerance) below the equipotential
      ! of the donor, or donor is below tolerance for detachment, then use regular roche_lobe scheme.
      if (b% rl_relative_gap(b% a_i) < 0 .and. &
          (b% rl_relative_gap(b% d_i)*slope_contact - b% rl_relative_gap(b% a_i) &
           > b% implicit_scheme_tolerance .or. &
           b% rl_relative_gap(b% d_i) < - b% implicit_scheme_tolerance)) then

         func = (b% rl_relative_gap(b% d_i) &
             + b% implicit_scheme_tolerance/2.0d0) * 2.0d0
      else
         if (q < 1d0) then
            func = (b% rl_relative_gap(b% d_i)*slope_contact - b% rl_relative_gap(b% a_i))
         else
            func = (b% rl_relative_gap(b% d_i) - b% rl_relative_gap(b% a_i)/slope_contact)
         end if
      end if
   end subroutine contact_condition

   subroutine my_contact_func_to_solve(binary_id, function_to_solve, use_sum, &
         detached, explicit, explicit_mdot, ierr)
      integer, intent(in) :: binary_id
      real(dp), intent(out) :: function_to_solve, explicit_mdot
      integer, intent(out) :: ierr
      logical, intent(out) :: use_sum, detached, explicit
      type (binary_info), pointer :: b

      real(dp) :: lq, rmax, rafromd, rdfroma, kolb_function_to_solve

      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      use_sum = .false.
      detached = .false.
      explicit = .false.
      explicit_mdot = 0d0

      if (b% point_mass_i /= 0) then
         ierr = -1
         write(*, *) "WARNING: contact scheme requires evolve_both_stars=.true."
         write(*, *) "Not transfering mass"
         return
      end if
      lq = log10(b% m(b% a_i) / b% m(b% d_i))

      ! check contact condition
      rafromd = eval_other_r(lq, b% r(b% d_i) / b% rl(b% d_i), ierr) - 1
      ! If accretor is overflowing its Roche lobe, then the contact scheme needs to be used.
      ! Otherwise, if accretor is (within tolerance) below the equipotential of the donor,
      ! or donor is below tolerance for detachment, then use regular roche_lobe scheme.
      if (b% rl_relative_gap(b% a_i) < 0 .and. &
            (rafromd - b% rl_relative_gap(b% a_i) > b% implicit_scheme_tolerance .or. &
                  b% rl_relative_gap(b% d_i) < -b% implicit_scheme_tolerance)) then
         function_to_solve = (b% rl_relative_gap(b% d_i) &
               + b% implicit_scheme_tolerance / 2.0d0) * 2.0d0
         if (function_to_solve < 0 .and. abs(b% mtransfer_rate) == 0) then
            detached = .true.
            return
         end if
      else
         call contact_condition(b, function_to_solve)
         use_sum = .true.
      end if

   end subroutine my_contact_func_to_solve

   ! Jdot
   subroutine my_extra_jdot(binary_id, ierr)
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b

      call binary_ptr(binary_id, b, ierr)
      if (.not. b% lxtra(ilx_doing_period_relaxation)) then
         b% extra_jdot = 0
      else
         b% extra_jdot = -deltaJ / deltaT
      end if

   end subroutine my_extra_jdot
