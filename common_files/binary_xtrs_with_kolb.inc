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
         ilx_had_rlof_at_zams = 4, ilx_had_ms_contact = 5, &
         ilx_had_caseA = 8, ilx_had_max_mdot = 9, ilx_has_merger = 10
   ! integer params
   integer, parameter :: iix_time0 = 1, iix_time1 = 2, iix_clockrate = 3
   ! real params
   integer, parameter :: ix_running_time = 1, ix_L2_thin = 3, ix_L2_thick = 4, &
         ix_L2_mdot = 5, ix_L3_thin = 6, ix_L3_thick = 7, ix_L3_mdot = 7

   logical :: inlist_L2_setting
   real(dp) :: inlist_fr, inlist_implicit_tolerance, inlist_fr_limit, inlist_fr_hard,&
          inlist_scale_max_correction, inlist_delta_HR, inlist_delta_lgTeff, inlist_delta_lgL

   integer :: slow_models = 0, consec_slow = 0
   logical :: using_kolb, too_slow = .false.
   logical, parameter :: use_Lout_MT = .false.
   real(dp) :: pre_ms_time, zams_age, deltaJ, deltaT, func_eval
   real(dp), parameter :: et_smoothing = 0.3  ! 0 = no smoothing, -> 1 = extreme smoothing

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

      b% other_adjust_mdots => my_adjust_mdots
      b% other_jdot_ml => my_jdot_ml
      b% other_extra_jdot => my_extra_jdot
      b% other_implicit_function_to_solve => my_contact_func_to_solve  ! uses r2(r1) relation from integrations of Fabry+2022
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
         b% terminate_if_L2_overflow = .false.  ! prevent l2 stoppage
         b% lxtra(ilx_pre_ms) = .true.
         b% lxtra(ilx_zams_l2_overflow) = .false.
         b% lxtra(ilx_doing_period_relaxation) = .false.
         b% lxtra(ilx_had_rlof_at_zams) = .false.
         b% lxtra(ilx_had_ms_contact) = .false.
         b% lxtra(ilx_had_ET) = .false.
         b% lxtra(ilx_had_caseA) = .false.
         b% lxtra(ilx_had_max_mdot) = .false.
         b% lxtra(ilx_has_merger) = .false.
         b% xtra(:) = 0d0
         b% ixtra(iix_time1) = 0d0
      end if

      pre_ms_time = pow(b% m1, -2.3d0) * 1e7  ! rough estimate of the prems lifetime

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
      if (b% lxtra(ilx_pre_ms) .or. b% lxtra(ilx_doing_period_relaxation)) then
         b% fr = -1d0
         return
      else
         b% fr = inlist_fr
         b% fr_hard = inlist_fr_hard
      end if

      !!! if here, binary must be past pre-ms !!!
      if (.not. use_Lout_MT) then
         b% terminate_if_L2_overflow = .true.
      else
         b% terminate_if_L2_overflow = .false.
      end if

!      ! smooth composition of outer 0.01 Msuns of accretor
!      b% s_accretor% smooth_outer_xa_big = 0.02
!      b% s_accretor% smooth_outer_xa_small = 0.01
!      b% s_donor% smooth_outer_xa_big = -1d0
!      b% s_donor% smooth_outer_xa_small = -1d0

      if (b% mtransfer_rate /= 0d0 .and. &
            max(abs(b% mdot_system_wind(b% d_i)/b% mtransfer_rate), &
                abs(b% mdot_system_wind(b% a_i)/b% mtransfer_rate)) > 5d-1) then

         ! if here, the mdot scheme will struggle to find a solution, either by
         ! looking to keep the donor exactly at the RL, or because the accretor
         ! has a total mdot of around zero
         ! let's relax fr so we smooth over the phase where mdot_trans \approx mdot_wind
         write(*, *) "high wind/mtransfer ratio", max(abs(b% mdot_system_wind(b% d_i)/b% mtransfer_rate), &
                abs(b% mdot_system_wind(b% a_i)/b% mtransfer_rate))
         new_limit = min(1d2, 1d1*abs(b% mdot_system_wind(b% d_i)/b% mtransfer_rate)) * inlist_fr
         new_limit = min(new_limit, 1d1*abs(b% mdot_system_wind(b% a_i)/b% mtransfer_rate) * inlist_fr)
         b% fr = min(new_limit, 1d-1)  ! cap it at 0.1
         b% fr_hard = min(5d-1, 2*inlist_fr_hard)
         write(*, *) "relaxing fr to", b% fr
      end if

      call set_hr_limits(b% s1)
      call set_hr_limits(b% s2)

      ! be more careful when we overflow L2 or when ET
      if (b% r(b% d_i) > 0.99 * min(Lout_radius(b), Lfar_radius(b))) then
         b% fr = 2d-1 * b% fr
      end if

      if (b% r(1) > b% rl(1) .and. b% r(2) > b% rl(2) .and. b% s1% use_other_energy) then
         b% fr = 5d-1 * b% fr
      end if

   end function extras_binary_start_step

   subroutine set_hr_limits(s)
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: ierr
      logical :: contact

      call binary_ptr(s% binary_id, b, ierr)
      contact = (b% r(1) >= b% rl(1) .and. b% r(2) >= b% rl(2))

      if (contact .and. s% use_other_energy) then  ! be more careful as ET turns on.
         s% delta_HR_limit = 2d-1 * inlist_delta_HR
         s% delta_lgTeff_limit = 2d-1 * inlist_delta_lgTeff
         s% delta_lgL_limit = 2d-1 * inlist_delta_lgL
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

      if (b% lxtra(ilx_doing_period_relaxation) .and. b% period < b% initial_period_in_days * secday) then
         extras_binary_check_model = retry
      end if

   end function extras_binary_check_model

   integer function extras_binary_finish_step(binary_id)
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      type (binary_info), pointer :: b
      type (star_info), pointer :: s1, s2
      integer, intent(in) :: binary_id
      integer :: ierr, id, i
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
         if (b% rl_relative_gap(1) >= 0d0  .or. b% rl_relative_gap(2) >= 0d0) then
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
         call binary_set_period_eccentricity(binary_id, b% initial_period_in_days * secday, 0d0, ierr)
         ! keep stars synchronized
         if (b% point_mass_i /= 1 .and. b% s1% rotation_flag) then
            call star_set_uniform_omega(b% s1% id, 2 * pi / b% period, ierr)
         end if
         if (b% point_mass_i /= 2 .and. b% s2% rotation_flag) then
            call star_set_uniform_omega(b% s2% id, 2 * pi / b% period, ierr)
         end if
      end if

      ! check for completion of period relaxation
      if (b% lxtra(ilx_doing_period_relaxation) .and. &
            b% period <= 1.01 * b% initial_period_in_days * secday) then
         b% lxtra(ilx_doing_period_relaxation) = .false.
      end if

      ! don't do anything below if we're pre-MS
      if (b% lxtra(ilx_pre_ms) .or. b% lxtra(ilx_doing_period_relaxation)) return

      !!! do post-step checks !!!
      ! check if case A MT started
      if (.not. b% lxtra(ilx_had_caseA) .and. abs(b% mtransfer_rate) >= 1d-10) then
         b% lxtra(ilx_had_caseA) = .true.
      end if

      ! mark if main sequence contact phase occurred
      if (.not. b% lxtra(ilx_had_ms_contact) .and. &
            .not. (b% s1% lxtra(ilx_depleted_H) .or. b% s2% lxtra(ilx_depleted_H)) .and. &
                  b% rl_relative_gap(1) > 0d0 .and. b% rl_relative_gap(2) > 0d0) then
         b% lxtra(ilx_had_ms_contact) = .true.
      end if

      ! terminate when mass loss becomes dynamical in L2 overflow
      if ((b% s1% star_mass/max(1d-99,abs(b% s1% star_mdot)) < 1d-1 * b% s1% kh_timescale .or. &
            b% s2% star_mass/max(1d-99,abs(b% s2% star_mdot)) < 1d-1 * b% s2% kh_timescale) .and. &
            b% r(b% d_i) > min(Lout_radius(b), Lfar_radius(b))) then
         b% lxtra(ilx_has_merger) = .true.
         extras_binary_finish_step = terminate
         write(*, *) "mass loss becoming dynamical (probably merger event starts), stopping MESA"
      end if

      ! terminate if no case A occurred
      if (b% s1% lxtra(ilx_depleted_H) .and. .not. b% lxtra(ilx_had_caseA)) then
         extras_binary_finish_step = terminate
         write(*, *) "no case A interaction!"
      end if

      ! terminate if max_mdot (shouldn't trigger? let's put it anyway to be sure...)
      if (b% mtransfer_rate <= -b% max_implicit_abs_mdot*Msun/secyer) then
         extras_binary_finish_step = terminate
         b% lxtra(ilx_had_max_mdot) = .true.
         write(*, *) "Terminate because of reaching max mdot_trans"
      end if

      ! set explicit mdot through Louts for next step
      ! we do it here because any redos in the next step do not restore MT/ET to these values
      ! because we rely on some hydro_vars and they aren't being reset before binary start step triggers
      call my_kolb_Lout_mdot(b% binary_id, dummy_mdot, ierr)

      ! set energy transfer stuff for coming step
      if (b% r(1) > b% rl(1) .and. b% r(2) > b% rl(2) .and. b% s1% use_other_energy) then
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
      how_many_extra_binary_history_columns = 14
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
         b% lxtra(ilx_had_ms_contact) = .true.
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
      else if (b% lxtra(ilx_has_merger)) then
         vals(6) = 3  ! ms l2 overflow / merger
      else if (.not. b% lxtra(ilx_had_caseA) .and. b% s1% lxtra(ilx_depleted_H)) then
         vals(6) = 4  ! no case A
      else if (b% lxtra(ilx_had_max_mdot)) then
         vals(6) = 5  ! max mdot
      else if (b% s1% termination_code == t_min_timestep_limit .or. &
            b% s2% termination_code == t_min_timestep_limit .or. &
            b% s1% termination_code == t_max_number_retries .or. &
            b% s2% termination_code == t_max_number_retries .or. &
            too_slow) then
         vals(6) = -99  ! problems before TAMS (no bueno)
      else  ! still running
         vals(6) = -1 ! still pre tams
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
      names(13) = 'had_ms_contact'
      vals(13) = bool_number(b% lxtra(ilx_had_ms_contact))
      names(14) = 'overflowing_L2'
      vals(14) = bool_number(b% r(b% d_i) > min(Lfar_radius(b), Lout_radius(b)))

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
      real(dp) :: new_l_trans_RL_1, new_l_trans_RL_2, new_l_trans_surf_1, new_l_trans_surf_2, p, scale

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

      write(*, *) "ET calculations, RL, surf", b% s1% xtra(ix_ET_RL) / Lsun, b% s1% xtra(ix_ET_surf) / Lsun

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

      ! get luminosities at the roche lobe
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
      !      write(*, *) "radii selected to base lum transfer on: ", &
      !         b% s1% r(roche_cell1), b% s2% r(roche_cell2)

      s_g1 = pi4 * standard_cgrav * b% s1% m(cell_1) * &
            eval_fp(lq, b% s1% r(cell_1) / b% rl(1), ierr) / &
            eval_ft(lq, b% s1% r(cell_1) / b% rl(1), ierr)
      s_g2 = pi4 * standard_cgrav * b% s2% m(cell_2) * &
            eval_fp(-lq, b% s2% r(cell_2) / b% rl(2), ierr) / &
            eval_ft(-lq, b% s2% r(cell_2) / b% rl(2), ierr)

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
      real(dp) :: lq, rafromd, rdfroma

      ierr = 0
      lq = log10(b% m(b% a_i) / b% m(b% d_i))

      if (lq < 0d0) then
         rafromd = eval_other_r(lq, b% r(b% d_i) / b% rl(b% d_i), ierr) - 1
         func = (rafromd - b% rl_relative_gap(b% a_i))
      else
         rdfroma = eval_other_r(-lq, b% r(b% a_i) / b% rl(b% a_i), ierr) - 1
         func = (b% rl_relative_gap(b% d_i) - rdfroma)
      end if
   end subroutine contact_condition

   subroutine my_contact_func_to_solve(binary_id, function_to_solve, use_sum, detached, explicit, explicit_mdot, ierr)
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
      ! Otherwise, if accretor radius is (within tolerance) below the equipotential
      ! of the donor, or donor is below tolerance for detachment, then use regular roche_lobe scheme.
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

   subroutine my_kolb_Lout_mdot(binary_id, mdot_extra, ierr)
      integer, intent(in) :: binary_id
      real(dp), intent(out) :: mdot_extra
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      if (.not. use_Lout_MT) then
         b% xtra(ix_L2_mdot) = 0d0
         b% xtra(ix_L3_mdot) = 0d0
      else ! compute MT rates
         call thin_mt_rate(b, 2, b% xtra(ix_L2_thin))
         call thick_mt_rate(b, 2, b% xtra(ix_L2_thick))
         call thin_mt_rate(b, 3, b% xtra(ix_L3_thin))
         call thick_mt_rate(b, 3, b% xtra(ix_L3_thick))

!         write(*, *) b% xtra(ix_L2_thin), b% xtra(ix_L2_thick), b% xtra(ix_L3_thin), b% xtra(ix_L3_thick)
         ! set the first mtransfer try on the previous MT + Lout loss - previous Lout loss
         b% xtra(ix_L2_mdot) = b% xtra(ix_L2_thin) + b% xtra(ix_L2_thick)
         b% xtra(ix_L3_mdot) = b% xtra(ix_L3_thin) + b% xtra(ix_L3_thick)
         write(*, *) "Computed L2,3 mdot, thin, thick", &
            abs(b% xtra(ix_L2_thin) / Msun * secyer), abs(b% xtra(ix_L2_thick) / Msun * secyer), &
            abs(b% xtra(ix_L3_thin) / Msun * secyer), abs(b% xtra(ix_L3_thick) / Msun * secyer)
      end if
   end subroutine my_kolb_Lout_mdot

   logical function near_or_not(lq, point)
      real(dp), intent(in) :: lq
      integer, intent(in) :: point
      if (lq <= 0) then
         if (point == 2) then
            near_or_not = .false.  ! other L point behind accretor
         elseif (point == 3) then
            near_or_not = .true.  ! other L point on this side of binary
         end if
      else
         if (point == 2) then
            near_or_not = .true.
         else
            near_or_not = .false.
         end if
      end if
   end function near_or_not

   subroutine thick_mt_rate(b, point, mdot_thick)
      type (binary_info), pointer :: b
      integer, intent(in) :: point  ! either 2 or 3 for L2 or L3
      real(dp), intent(out) :: mdot_thick
      real(dp) :: radius, rl_radius, rlother_radius
      integer :: indexRL, indexLother
      real(dp) :: lq, a, G1, F, dr_rl, r, S, dS
      integer :: i, ierr
      logical :: near  ! is the L point near (T) or far (F)?

      mdot_thick = 0d0
      if (b% rl_relative_gap(b% d_i) < 0.0d0) then
         return  ! not overflowing L1
      end if

      radius = b% r(b% d_i)
      rl_radius = b% rl(b% d_i)
      lq = log10(b% m(b% a_i) / b% m(b% d_i))
      a = b% separation

      near = near_or_not(lq, point)

      if (near) then
         rlother_radius = Lout_radius(b)
      else
         rlother_radius = Lfar_radius(b)
      end if

      indexRL = 0
      indexLother = 0
      i = 1
      do while (b% s_donor% r(i) > rl_radius)  ! search RL index
         if (b% s_donor% r(i) < rlother_radius .and. indexLother == 0) then
            if (i == 1) then
               indexLother = -1  ! no overflow of other L point
            else
               indexLother = i - 1  ! index just above L other equipotential
            end if
         end if
         i = i + 1
      end do
      if (i == 1) then  ! shouldn't happen (not overflowing L1)
         return
      else
         indexRL = i - 1  ! mark zone just above formal RL
      end if

      !--------------------- Optically thick MT rate ------------------------
      do i = indexRL, 1, -1
         G1 = b% s_donor% gamma1(i)
         F = sqrt(2 * G1 / (G1 + 1)) * pow(2 / (G1 + 1), 1 / (G1 - 1))
         r = b% s_donor% r(i)
         if (r < rlother_radius) cycle  ! integrate only layers from Lout up
         mdot_thick = mdot_thick + Lother_integrand(i)
      end do
      ! mdot_thick is defined negative
      mdot_thick = -mdot_thick
!      write(*, *) "L2 kolb thick", log10(abs(mdot_thick / Msun * secyer))
!            write(*, *) radius / Rsun, rl_radius / Rsun, rl2_radius / Rsun

   contains

      real(dp) function Lother_integrand(i)
         integer, intent(in) :: i
         integer :: j

         if (i == indexLother) then
            dr_rl = (r / rlother_radius - 1)  ! only fraction of the cell overflows other L point
            S = 0d0  ! area at the L point is zero of course
         else
            dr_rl = (r - b% s_donor% r(i + 1)) / rlother_radius
            if (near) then
               S = eval_Lout_split_area(lq, b% s_donor% r(i + 1) / rlother_radius, ierr) * a * a
            else
               S = eval_Lfar_split_area(lq, b% s_donor% r(i + 1) / rlother_radius, ierr) * a * a
            end if
         end if

         if (near) then
            dS = eval_Lout_split_area(lq, r / rlother_radius, ierr) * a * a - S
         else
            dS = eval_Lfar_split_area(lq, r / rlother_radius, ierr) * a * a - S
         end if

         j = 1
         do while (dS <= 0)  ! make sure area is positive
            if (r / rlother_radius + j * dr_rl > 1.5d0) then
               write(*, *) 'warning: radii > 1.5 r_rl needed in kolb thick, setting mdot for this layer to 0'
               dS = 0d0
               exit
            end if
            j = j + 1
            if (near) then
               dS = (eval_Lout_split_area(lq, r / rlother_radius + j * dr_rl, ierr) * a * a - S) / j
            else
               dS = (eval_Lfar_split_area(lq, r / rlother_radius + j * dr_rl, ierr) * a * a - S) / j
            end if
            if (j >= 100) then
               write(*, *) "cannot find proper numerical derivative of Lout split", exp10(lq)
               exit
            end if
         end do

         Lother_integrand = sqrt(b% s_donor% Peos(i) * b% s_donor% rho(i)) * F * dS

      end function Lother_integrand

   end subroutine thick_mt_rate

   subroutine thin_mt_rate(b, point, mdot_thin)
      type(binary_info), pointer :: b
      real(dp), intent(out) :: mdot_thin
      integer, intent(in) :: point  ! 2 or 3 for other L point
      real(dp) :: normalization, rl1, ritter_exponent, psi, psiL, mdot_thin_0, lq, v_th2, r, a, &
            dS_dr, dPsi_dr, rl_frac
      integer :: ierr, j
      logical :: near

      lq = log10(b% m(b% a_i) / b% m(b% d_i))
      rl1 = b% rl(b% d_i)
      r = b% r(b% d_i)
      a = b% separation

      near = near_or_not(lq, point)

      normalization = standard_cgrav * b% m(b% d_i) / a
      ! compute differential area at the lagrangian point
      if (near) then
         rl_frac = Lout_radius(b) / rl1
      else
         rl_frac = Lfar_radius(b) / rl1
      end if
      if (near) then
         dS_dr = eval_Lout_split_area(lq, 1d0 + nudge, ierr) * a * a / nudge
      else
         dS_dr = eval_Lfar_split_area(lq, 1d0 + nudge, ierr) * a * a / nudge
      end if
      dPsi_dr = (eval_psi(lq, rl_frac + nudge, ierr) - eval_psi(lq, rl_frac, ierr)) * normalization / nudge
      j = 1
      do while (dS_dr <= 0)  ! make sure numerical derivative is positive
         j = j + 1
         if (near) then
            dS_dr = eval_Lout_split_area(lq, 1d0 + nudge, ierr) * a * a / nudge
         else
            dS_dr = eval_Lfar_split_area(lq, 1d0 + nudge, ierr) * a * a / nudge
         end if
         if (j >= 100) then
            write(*, *) "cannot find proper numerical derivative of Lout split", exp10(lq)
            if (using_kolb) then
               stop
            else
               exit
            end if
         end if
      end do

      v_th2 = kerg * b% s_donor% T(1) / (mp * b% s_donor% mu(1))  ! thermal sound speed
      mdot_thin_0 = b% s_donor% rho(1) * pow(v_th2, 1.5d0) * pow(exp(1d0), -0.5d0) * dS_dr / dPsi_dr

      psiL = eval_psi(lq, rl_frac, ierr) * normalization
      psi = eval_psi(lq, r / rl1, ierr) * normalization
      ritter_exponent = - max(psiL - psi, 0d0) / v_th2  ! psiL > psi should be true for thin Mdot, so ritter_exp <= 0

      mdot_thin = - mdot_thin_0 * exp(ritter_exponent)
!      write(*, *) "L2 kolb thin", mdot_thin / Msun * secyer, ritter_exponent, log10(abs(mdot_thin_0 / Msun * secyer))
!      write(*, *) dS_dr, eval_L1_split_area(lq, 1d0 + nudge, ierr) * a * a, sqrt(v_th2) / b% s1% csound(1), b% s_donor% rho(1)

   end subroutine thin_mt_rate

   subroutine my_adjust_mdots(binary_id, ierr)
      integer, intent(in) :: binary_id
      real(dp) :: actual_mtransfer_rate
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      logical :: near
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      actual_mtransfer_rate = 0d0

      b% fixed_xfer_fraction = 1 - b% mass_transfer_alpha - b% mass_transfer_beta - &
            b% mass_transfer_delta

      if (.not. b% use_other_mdot_edd) then
         call binary_eval_mdot_edd(b% binary_id, b% mdot_edd, b% mdot_edd_eta, ierr)
      else
         call b% other_mdot_edd(b% binary_id, b% mdot_edd, b% mdot_edd_eta, ierr)
      end if

      ! solve wind mass transfer
      ! b% mdot_wind_transfer(b% d_i) is a negative number that gives the
      ! amount of mass transferred by unit time from the donor to the
      ! accretor.
      call binary_eval_wind_xfer_fractions(b% binary_id, ierr)
      if (ierr/=0) then
         write(*, *) "Error in eval_wind_xfer_fractions"
         return
      end if
      b% mdot_wind_transfer(b% d_i) = b% s_donor% mstar_dot * &
            b% wind_xfer_fraction(b% d_i)
      if (b% point_mass_i == 0) then
         b% mdot_wind_transfer(b% a_i) = b% s_accretor% mstar_dot * &
               b% wind_xfer_fraction(b% a_i)
      else
         b% mdot_wind_transfer(b% a_i) = 0d0
      end if

      near = near_or_not(log10(b% m(b% a_i) / b% m(b% d_i)), 2)

      ! Set mdot for the donor
      if (near) then
         b% s_donor% mstar_dot = b% s_donor% mstar_dot + b% mtransfer_rate - &
               b% mdot_wind_transfer(b% a_i) + b% xtra(ix_L2_mdot)
      else
         b% s_donor% mstar_dot = b% s_donor% mstar_dot + b% mtransfer_rate - &
               b% mdot_wind_transfer(b% a_i) + b% xtra(ix_L3_mdot)
               end if

!      write(*, *) "donor mstar_dot ", b% s_donor% mstar_dot / Msun * secyer, b% xtra(ix_L2_mdot) / Msun * secyer

      ! Set mdot for the accretor
      if (b% point_mass_i == 0) then
         ! do not repeat if using the implicit wind
         if (.not. (b% num_tries > 0 .and. b% s_accretor% was_in_implicit_wind_limit)) then
            b% accretion_mode = 0
            b% acc_am_div_kep_am = 0.0d0
            if (near) then
               b% s_accretor% mstar_dot = b% s_accretor% mstar_dot - &
                     b% mtransfer_rate * b% fixed_xfer_fraction - b% mdot_wind_transfer(b% d_i) + &
                     b% xtra(ix_L3_mdot)
            else
               b% s_accretor% mstar_dot = b% s_accretor% mstar_dot - &
                     b% mtransfer_rate * b% fixed_xfer_fraction - b% mdot_wind_transfer(b% d_i) + &
                     b% xtra(ix_L2_mdot)
            end if
!            write(*, *) "acctr mstar_dot ", b% s_accretor% mstar_dot / Msun * secyer, b% xtra(ix_L3_mdot) / Msun * secyer

            !set angular momentum accretion as described in A.3.3 of de Mink et al. 2013
            if (b% do_j_accretion) then
               if (.not. b% use_other_accreted_material_j) then
                  call binary_eval_accreted_material_j(b% binary_id, ierr)
               else
                  call b% other_accreted_material_j(b% binary_id, ierr)
               end if
               if (ierr /= 0) then
                  write(*, *) 'error in accreted_material_j'
                  return
               end if
            end if
         end if
         b% accretion_luminosity = 0d0 !only set for point mass
      else
         ! accretor is a point mass
         if (.not. b% model_twins_flag) then
            !combine wind and RLOF mass transfer
            actual_mtransfer_rate = b% mtransfer_rate * b% fixed_xfer_fraction + b% mdot_wind_transfer(b% d_i) !defined negative
            b% component_mdot(b% a_i) = -actual_mtransfer_rate
            ! restrict accretion to the Eddington limit
            if (b% limit_retention_by_mdot_edd .and. b% component_mdot(b% a_i) > b% mdot_edd) then
               b% component_mdot(b% a_i) = b% mdot_edd ! remove all accretion above the edd limit
            end if
            b% accretion_luminosity = &
                  b% mdot_edd_eta * b% component_mdot(b% a_i) * clight * clight
            ! remove rest mass radiated away
            if (b% use_radiation_corrected_transfer_rate) then
               b% component_mdot(b% a_i) = (1 - b% mdot_edd_eta) * b% component_mdot(b% a_i)
            end if
         end if
      end if

      ! mdot_system_transfer is mass lost from the vicinity of each star
      ! due to inefficient rlof mass transfer, mdot_system_cct is mass lost
      ! from a circumbinary coplanar toroid.
      if (b% mtransfer_rate + b% mdot_wind_transfer(b% d_i) >= 0 .or. b% CE_flag) then
         b% mdot_system_transfer(b% d_i) = 0d0
         b% mdot_system_transfer(b% a_i) = 0d0
         b% mdot_system_cct = 0d0
      else
         b% mdot_system_transfer(b% d_i) = b% mtransfer_rate * b% mass_transfer_alpha
         b% mdot_system_cct = b% mtransfer_rate * b% mass_transfer_delta
         if (b% point_mass_i == 0 .or. b% model_twins_flag) then  ! no point mass or twins
            b% mdot_system_transfer(b% a_i) = b% mtransfer_rate * b% mass_transfer_beta
         else  ! point mass
            ! do not compute mass lost from the accretor using just mass_transfer_beta, as
            ! mass transfer can be stopped also by going past the Eddington limit
            b% mdot_system_transfer(b% a_i) = (actual_mtransfer_rate + b% component_mdot(b% a_i)) &
                  + b% mtransfer_rate * b% mass_transfer_beta
         end if
      end if

   end subroutine my_adjust_mdots

   ! Jdot
   subroutine my_jdot_ml(binary_id, ierr)
      integer, intent(in) :: binary_id
      integer, intent(out) :: ierr
      type (binary_info), pointer :: b
      real(dp) :: x_L2, x_L3, mu, om_root_e2
      ierr = 0
      call binary_ptr(binary_id, b, ierr)
      if (ierr /= 0) then
         write(*, *) 'failed in binary_ptr'
         return
      end if

      om_root_e2 = 2 * pi / b% period * sqrt(1 - pow2(b% eccentricity))
      mu = b% m(b% a_i) / (b% m(b% a_i) + b% m(b% d_i))  ! dim-less center of mass
      ! mass lost from vicinity of donor
      b% jdot_ml = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i)) * &
            pow2(mu * b% separation) * om_root_e2
      ! mass lost from vicinity of accretor
      b% jdot_ml = b% jdot_ml + (b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i)) * &
            pow2((1 - mu) * b% separation) * om_root_e2
      ! mass lost from circumbinary coplanar toroid
      b% jdot_ml = b% jdot_ml + b% mdot_system_cct * b% mass_transfer_gamma * &
            sqrt(b% s_donor% cgrav(1) * (b% m(1) + b% m(2)) * b% separation)

      if (use_Lout_MT) then
         if (b% m(b% a_i) > b% m(b% d_i)) then
            x_L2 = find_x_Lout(b)
            x_L3 = find_x_Lfar(b)
         else
            x_L2 = find_x_Lfar(b)
            x_L3 = find_x_Lout(b)
         end if
         b% jdot_ml = b% jdot_ml + &
               b% xtra(ix_L2_mdot) * pow2((mu - x_L2) * b% separation) * om_root_e2 +&
               b% xtra(ix_L3_mdot) * pow2((mu - x_L3) * b% separation) * om_root_e2
      end if

   end subroutine my_jdot_ml

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

   ! helper functions for RL/geometry
   real(dp) function find_x_Lout(b) result (x_lout)
      ! finds x coordinate of outer Lagrangian point of the donor (wrt to center of donor)
      ! in units of separation

      type (binary_info), pointer :: b

      integer, pointer :: ipar(:) ! (lipar)  ! not used
      real(dp), pointer :: rpar(:) ! (lrpar)
      integer :: ierr

      allocate(rpar(1))
      rpar(1) = b% m(b% a_i) / b% m(b% d_i)
      x_lout = safe_root_with_guess(f, -1d-1, 1d-3, &  ! function, guess, dx for bracket
            -1d0, 0d0, & ! left, right bracket
            my_f(-1d0, rpar(1)), my_f(0d0, rpar(1)), & ! f(left, right bracket)
            25, 50, 1d-4, 1d-6, & ! i_newton, i_max, x_tol, y_tol
            1, rpar, 0, ipar, & ! func_params
            ierr)

   contains

      real(dp) function my_f(x, q)
         real(dp), intent(in) :: x, q
         my_f = (1 + q) * pow5(x) - (2 + 3 * q) * pow4(x) + (&
               1 + 3 * q) * pow3(x) + pow2(x) - 2 * x + 1
      end function my_f

      real(dp) function f(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
         real(dp), intent(in) :: x
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(out) :: dfdx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         real(dp) :: q

         q = rpar(1)
         f = my_f(x, q)
         dfdx = 5 * (1 + q) * pow4(x) - 4 * (2 + 3 * q) * pow3(x) + 3 * (&
               1 + 3 * q) * pow2(x) + 2 * x - 2
         ierr = 0
      end function f

   end function find_x_Lout

   real(dp) function find_x_Lfar(b) result (x_lfar)
      ! finds x coordinate of far Lagrangian point of the donor (wrt to center of donor)
      ! in units of separation

      type (binary_info), pointer :: b

      integer, pointer :: ipar(:) ! (lipar)  ! not used
      real(dp), pointer :: rpar(:) ! (lrpar)
      integer :: ierr

      allocate(rpar(1))
      rpar(1) = b% m(b% a_i) / b% m(b% d_i)
      x_lfar = safe_root_with_guess(f, -1d-1, 1d-3, &  ! function, guess, dx for bracket
            1d0, 2d0, & ! left, right bracket
            my_f(1d0, rpar(1)), my_f(2d0, rpar(1)), & ! f(left, right bracket)
            25, 50, 1d-4, 1d-6, & ! i_newton, i_max, x_tol, y_tol
            1, rpar, 0, ipar, & ! func_params
            ierr)

   contains

      real(dp) function my_f(x, q)
         real(dp), intent(in) :: x, q
         my_f = (1 + q) * pow5(x) - (2 + 3 * q) * pow4(x) + (&
               1 + 3 * q) * pow3(x) - (1 + 2 * q) * pow2(x) + 2 * x - 1
      end function my_f

      real(dp) function f(x, dfdx, lrpar, rpar, lipar, ipar, ierr)
         real(dp), intent(in) :: x
         integer, intent(in) :: lrpar, lipar
         real(dp), intent(out) :: dfdx
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr
         real(dp) :: q

         q = rpar(1)
         f = my_f(x, q)
         dfdx = 5 * (1 + q) * pow4(x) - 4 * (2 + 3 * q) * pow3(x) + 3 * (&
               1 + 3 * q) * pow2(x) - 2 * (1 + 2 * q) * x + 2
         ierr = 0
      end function f

   end function find_x_Lfar

   real(dp) function Lout_radius(b)
      ! finds the equivalent radius of donor at which it would fill up to its outer L point, in cm
      type (binary_info), pointer :: b
      real(dp) :: sig, gap, q

      q = b% m(b% a_i) / b% m(b% d_i)
      sig = 62.9237d0 / (15.9839d0 + pow(q, 0.2240d0))
      gap = 3.3752d0 / (1 + ((log(q) + 1.0105d0) / sig) ** 2) / (9.0087d0 + pow(q, -0.4022d0))
      Lout_radius = (1 + gap) * b% rl(b% d_i)

   end function Lout_radius

   real(dp) function Lfar_radius(b)
      ! finds the equivalent radius of donor at which it would fill up to its far L point, in cm
      type (binary_info), pointer :: b
      real(dp) :: gap, q

      q = b% m(b% a_i) / b% m(b% d_i)
      gap = 1.284d0 / (2.168d0 * pow(q, -0.8541d0) + log(pow(q, -3.076d0) + 1d0) + 1.169d0)
      Lfar_radius = (1 + gap) * b% rl(b% d_i)
   end function Lfar_radius
