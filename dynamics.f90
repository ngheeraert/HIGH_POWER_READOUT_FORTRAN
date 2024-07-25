MODULE dynamics

	USE systm, only: param, state, traj, allocate_state, allocate_trajectory, initialise_bare_state, &
		parameterchar, norm, normalise, lvl_occupation, energy, photon_nb_in_sig_mode, &
		print_log, ov, ov_scalar, Ic, pi, update_sums, renyi_entropy, initialise_dressed_state, &
		ds_occupation, apadag_sig_mode!, calculate_derivatives
	USE typedefs, only : cx => c_type, rl => r_type
	USE lapackblas, only : diagonalise_hermitian_zhpev, diagonalise_complex_zgeev


	implicit none

CONTAINS

	SUBROUTINE compute_trajectory( sys ) 
		type(param), intent(in out)				::  sys
		type(param)								::  sys2
		type(state)								::  st
		type(state)								::  ost
		type(traj)								::  tr
		integer									::  min_slow_fact_exp, max_slow_fact_exp, slow_fact_exp
		integer									::  adding_min_slow_fact_exp, adding_max_slow_fact_exp
		real(rl), dimension(11)					::  print_times
		real(rl), dimension( sys%median_pts )   ::  error_buffer 
		real(rl)							    ::  t_first_print, t_final_print !-- for time printing on console
		real(rl)							    ::  t_last_print, t_last_error_calc, t_initial, t1, t2, t_final, tref
		real(rl)							    ::  dt_add, dt_add_incerement
		character(len=100)						::  log_message
		type(state) 		    				::  st_before_adding
		integer 								::  info, info_2, adding_i_time, incr_count, j
		real(rl)								::  dt_add_increment, adding_t_last_print, adding_tref, entr, error_val
		logical									::  back_tracking_allowed
		integer									::  i, p, s, n, k
		real(8)									::  cumu_solve_time, cumu_other_time, cumu_sum_time
		real(8)									::  full_cumu_solve_time, full_cumu_other_time,&
														full_cumu_sum_time, full_cumu_err_time, t_err_1, t_err_2
		real(8)									::  old_full_cumu_solve_time, old_full_cumu_other_time, &
														old_full_cumu_sum_time,old_full_cumu_err_time
		complex(cx), allocatable				::  eigen_vals_0(:), eigen_vals_1(:)


		!======================================================== 
		!-- Allocation and initialisation
		!=======================================================

		CALL allocate_state(sys,st)
		CALL print_log( sys, "-- arrays allocated", sys%logfile )
		CALL allocate_trajectory(sys, tr, 2*int(sys%tmax/sys%print_delay)  )
		CALL print_log( sys, "-- trajectory allocated", sys%logfile )

		if (sys%qb_ini < 0) then
			print*,'-- ds occupation gs'
			print*, ds_occupation( sys, st, 1)
			print*,'-- ds occupation es'
			print*, ds_occupation( sys, st, 2)
			CALL initialise_dressed_state(sys,st,0._rl)
			CALL print_log( sys, "-- DRESSED state initialisation complete", sys%logfile )
		else
			CALL initialise_bare_state(sys,st,0._rl)
			CALL print_log( sys, "-- BARE state initialisation complete", sys%logfile )
		end if

		old_full_cumu_solve_time = 0._8
		old_full_cumu_other_time = 0._8
		old_full_cumu_sum_time = 0._8
		old_full_cumu_err_time = 0._8
		full_cumu_solve_time = 0._8
		full_cumu_other_time = 0._8
		full_cumu_sum_time = 0._8
		full_cumu_err_time = 0._8
		cumu_solve_time = 0._8
		cumu_other_time = 0._8
		sys2 = sys
		dt_add_increment = 0.0
		slow_fact_exp = 0
		min_slow_fact_exp = 30
		max_slow_fact_exp = 0
		print_times = 0._rl
		t_first_print =  st%t - 1e-8
		t_final_print =  sys%tmax - 1e-8
		!p2_above_nb = 0

		do i=1, size( print_times,1 )-2
			print_times(i) = t_first_print &
				+ (i-1)*( t_final_print - t_first_print ) / dble( size( print_times,1 ) - 1 )
		end do
		print_times( size(print_times,1)-1 ) = sys%tmax - 1.1*sys%print_delay
		print_times( size(print_times,1) ) = 1e8

		tref = 0._rl
		t_last_print = -1e10
		t_last_error_calc = -1e10
		error_buffer(:) = 0._rl
		!ps_buffer(:,:) = 0._rl
		p=1
		CALL CPU_TIME(t_initial)
		t1 = t_initial

		!======================================
		!== TIME-EVOLUTION LOOP
		!======================================
		write(log_message,'( a6, a4,a8, a4,a3, a4,a3, a4,a3, a4,a3, a4,a3, a4,a3, a4,a3)') ' TIME ',&
			' || ', ' DEL t  ',&
			' || ',  'err ',&
			' || ',  'solve ',&
			' || ', 'other ',&
			' || ', ' sum  ',&
			' || ', 'ncs',&
			' || ', 'SFE',' -- ', 'SFE'
		CALL print_log( sys, log_message, sys%logfile )
		write(log_message,'( a6, a4,a8, a4,a3, a4,a3, a4,a3, a4,a3, a4,a3, a4,a3, a4,a3)') '------',&
			' || ', '--------',&
			' || ', '--------',&
			' || ',  '------',&
			' || ', '------',&
			' || ', '------',&
			' || ', '---',&
			' || ', '---','----', '---'
		CALL print_log( sys, log_message, sys%logfile )

		WHILE_DO: DO 

			IF ( st%t > sys%tmax )  exit

			if ( (sys%time_bar .eqv. .true.) .and. (print_times(p) < st%t) ) then
				p = p+1
				CALL CPU_TIME(t2)
				write(log_message,'(f6.1,a4,I8,a4,I3,a4,I3,a4,I3,a4,I3,a4,I3,a4,I3,a4,I3)') st%t,&
					' || ', int(t2 - t1),&
					' || ', int( 100*(full_cumu_err_time - old_full_cumu_err_time)/(t2 - t1) ),&
					' || ', int( 100*(full_cumu_solve_time - old_full_cumu_solve_time)/(t2 - t1) ),&
					' || ', int( 100*(full_cumu_other_time - old_full_cumu_other_time)/(t2 - t1) ),&
					' || ', int( 100*(full_cumu_sum_time - old_full_cumu_sum_time)/(t2 - t1) ),&
					' || ',st%ncs,&
					' || ', min_slow_fact_exp,' -- ', max_slow_fact_exp
				CALL print_log( sys, log_message, sys%logfile )

				max_slow_fact_exp = 0
				min_slow_fact_exp = 30

				old_full_cumu_solve_time = full_cumu_solve_time
				old_full_cumu_other_time = full_cumu_other_time
				old_full_cumu_sum_time = full_cumu_sum_time
				old_full_cumu_err_time = full_cumu_err_time
				t1 = t2

				!if ( (st%t>0.1*sys%tmax) .and. (st%t<0.9*sys%tmax) .and. (modulo(p,2)==0) )  then
				!	sys2%tmax = st%t
				!	sys2%ncs_max = st%ncs
				!	CALL print_evolution_data(sys2,tr)
				!	CALL print_yks( sys2, st  )
				!	CALL print_ps( sys2, st  )
				!end if

			end if

			!=========================================================================
			!== data printing
			if ( st%t > t_last_print + sys%print_delay ) then

				tr%i_time = tr%i_time + 1
				tr%time_ar( tr%i_time ) = st%t

				if ( (sys%calc_error==1) &
						.and. (st%t > t_last_error_calc + sys%error_calc_delay) ) then
					CALL CPU_TIME(t_err_1)
					error_val = error( sys, st )
					CALL CPU_TIME(t_err_2)
					full_cumu_err_time = full_cumu_err_time + (t_err_2 - t_err_1)
					error_buffer = CSHIFT( error_buffer, -1 )
					tr%error_ar( tr%i_time,1 ) = median( error_buffer )
					t_last_error_calc = st%t
				end if
				error_buffer(1) = real( error_val )

				do s=1, sys%nl
					if (sys%qb_ini < 0) then
						tr%populations_ar( tr%i_time, s, 1 ) = ds_occupation( sys, st, s )
					else
						tr%populations_ar( tr%i_time, s, 2 ) = lvl_occupation( st, s )
					end if
					do n=1, st%ncs
						tr%ps2_ar( tr%i_time, s, n ) = abs( st%p(s,n) )**2
						tr%y02_ar( tr%i_time, s, n ) = abs( sum( st%y(s,n,:)*sys%Omat(1,:) )  )**2
					end do
				end do
				tr%norm_ar( tr%i_time ) = norm(st)
				tr%renyi_entropy( tr%i_time ) = renyi_entropy(sys, st, 1, 8 )
				tr%apadag_ar( tr%i_time ) = apadag_sig_mode( sys, st, 1 )
				do k=1, sys%nmodes
					tr%n_ar( tr%i_time, k ) = photon_nb_in_sig_mode( sys, st, k )
				end do
				t_last_print = st%t

				allocate( eigen_vals_0( st%ncs ) )
				allocate( eigen_vals_1( st%ncs ) )
				!call diagonalise_hermitian_zhpev( st%ovm(1,:,1,:), eigen_vals_0 )
				!call diagonalise_hermitian_zhpev( st%ovm(2,:,2,:), eigen_vals_1 )
				call diagonalise_complex_zgeev( st%ovm(1,:,1,:), eigen_vals_0 )
				call diagonalise_complex_zgeev( st%ovm(2,:,2,:), eigen_vals_1 )
				tr%ovm_eval_0_ar( tr%i_time, :st%ncs ) = real( eigen_vals_0 )
				tr%ovm_eval_1_ar( tr%i_time, :st%ncs ) = real( eigen_vals_1 )
				deallocate( eigen_vals_0 )
				deallocate( eigen_vals_1 )

			end if

			if (st%ncs == sys%ncs_ini) then
				back_tracking_allowed = .true.
			else
				if (st%ncs > st_before_adding%ncs) then
					back_tracking_allowed = .true.
				else
					back_tracking_allowed = .false.
				end if
			end if

			!-- define the minimumm time delay dt_add beofre adding a new coherent state
			dt_add = sys%dtadd_arr(st%ncs) + dt_add_increment

			!-- evolve until next time printing time, or info > 1.
			!CALL time_evolve_RK45( sys, st, tref, dt_add, min_slow_fact_exp, max_slow_fact_exp, &
			!	slow_fact_exp, info, back_tracking_allowed, error_buffer )


			CALL time_evolve_RK45_no_buffer( sys, st, tref, dt_add, min_slow_fact_exp, max_slow_fact_exp, &
				slow_fact_exp, info, back_tracking_allowed, cumu_solve_time, cumu_other_time, cumu_sum_time )

			full_cumu_solve_time = full_cumu_solve_time + cumu_solve_time
			full_cumu_other_time = full_cumu_other_time + cumu_other_time
			full_cumu_sum_time = full_cumu_sum_time + cumu_sum_time


			!-- if the evolution did not succeed, backtrack in time, to state before adding
			if ( (info == 2) .and. (st%ncs>sys%ncs_ini) )   then

				st = st_before_adding
				tr%i_time = adding_i_time
				t_last_print = adding_t_last_print
				tref = adding_tref
				max_slow_fact_exp = adding_max_slow_fact_exp
				min_slow_fact_exp = adding_min_slow_fact_exp

				!-- making sure cs is not directly added
				incr_count=0
				DO_INCREMENT: DO

					IF ( st%t < tref + sys%dtadd_arr(st%ncs) + dt_add_increment - 0.009  ) THEN
						CALL print_log( sys, "-- dt_add incremented.", sys%logfile )
						EXIT DO_INCREMENT
					END IF

					dt_add_increment = dt_add_increment + 0.01
					!print*,'---------------------------------------------'
					!print*,'-- increment added -> tref =', tref 
					!print*,'-- increment added -> dt_incr =', dt_add_increment 
					!print*,'-- increment added -> sys%dtadd_arr(st%ncs)  =', sys%dtadd_arr(st%ncs) 
					!print*,'-- increment added -> dt_add =', sys%dtadd_arr(st%ncs) + dt_add_increment 
					incr_count = incr_count+1

				END DO DO_INCREMENT

			elseif ( (info == 2) .and. (st%ncs==sys%ncs_ini) ) then
				sys%offset_ang = sys%offset_ang + 0.1*2*pi/( sqrt(3.0) )
				tref = 0._rl
				t_last_print = -1e5
				error_buffer(:) = 0._rl
				tr%i_time = 1
				CALL print_log( sys, "-- offset angle incremented.", sys%logfile )
				CALL initialise_bare_state(sys,st,0._rl)
				CALL print_log( sys, "-- state RE-initialisation complete.", sys%logfile )
			end if
			
			!if ( all_p2_above_threshold( sys, st ) .eqv. .true. ) then
			!	p2_above_nb = p2_above_nb + 1
			!else
			!	p2_above_nb = 0
			!end if

			!-- define the waiting time before adding interms of ncs
			dt_add = sys%dtadd_arr(st%ncs) + dt_add_increment

			!-- add a coherent state in case the conditions are met
			if ( ( median(error_buffer) > sys%error_thr ) &
						.and. ( st%t > tref + dt_add ) &
						.and. ( st%ncs < sys%ncs_max ) ) then

				adding_tref = tref
				adding_i_time = tr%i_time
				adding_t_last_print = t_last_print
				st_before_adding = st
				adding_max_slow_fact_exp = max_slow_fact_exp
				adding_min_slow_fact_exp = min_slow_fact_exp

				CALL add_coherent_state( sys, st )
				tref = st%t

				!CALL evolve_new_coherent_state( sys, st, st%t + sys%dt, sys%dt/100, 1e-8_rl )
				!print*, 'new cs evolve complete'
				!stop

			end if

		END DO WHILE_DO

		CALL CPU_TIME(t_final)
		write(log_message,'(a20,I9)') 'Total running time= ', int( t_final - t_initial )
		CALL print_log( sys, log_message, sys%logfile )

		CALL print_evolution_data(sys,tr)
		CALL print_yks( sys, st  )
		CALL print_ps( sys, st  )

		!CALL print_yks_DB( sys, st  )
		!CALL print_omat(sys)

		if (sys%logfile>=1) then
			close(100)
		end if
	END SUBROUTINE

	SUBROUTINE time_evolve_RK45_no_buffer( sys, st, tref, dt_add, min_slow_fact_exp, max_slow_fact_exp, slow_fact_exp,&
			info, bt_allowed, cumu_solve_time, cumu_other_time, cumu_sum_time )
		type(param), intent(in)       		:: 	sys
		type(state), intent(in out)       	:: 	st
		real(8), intent(in)       			:: 	tref, dt_add
		integer, intent(in out)  	  		:: 	min_slow_fact_exp, max_slow_fact_exp, slow_fact_exp
		integer, intent(out)  		  		:: 	info
		real(8), intent(out)  		  		:: 	cumu_solve_time, cumu_other_time, cumu_sum_time
		logical, intent(in)					::  bt_allowed
		type(state)       					:: 	midst
		real(8)      			  			:: 	dt, error_RK5, t_corrected,&
													last_error_calc, t0, t_solve, t_other, t_sum
		logical								::  dt_corrected
		real(8)					  			:: 	a_mat( 5, 5 ), c_mat( 5 ), b_rk4( 6 ), b_rk5( 6 )
		complex(cx), dimension(size(st%p,1),size(st%p,2))   :: 	kp1, kp2, kp3, kp4, kp5, kp6
		complex(cx), dimension(size(st%y,1), size(st%y,2), size(st%y,3))::  kf1, kf2, kf3, kf4, kf5, kf6
		complex(cx), dimension(size(st%p,1),size(st%p,2))             		::  pdot_RK4, pdot_RK5
		complex(cx), dimension(size(st%y,1), size(st%y,2), size(st%y,3))     ::  fdot_RK4, fdot_RK5

		t0 = st%t
		midst = st
		t_corrected = -1.0e9_8
		dt_corrected = .false.
		last_error_calc = -1e3
		info = 1

		a_mat = transpose( reshape( (/ 1._8/4._8 , 0._8 , 0._8 , 0._8 , 0._8,&
			3._8/32._8 , 9._8/32._8 , 0._8 , 0._8 , 0._8,&
			1932._8/2197._8 , -7200._8/2197._8 , 7296._8/2197._8 , 0._8, 0._8,&
			439._8/216._8 , -8._8 , 3680._8/513._8 , -845._8/4104._8 , 0._8,&
			-8._8/27._8 , 2._8 , -3544._8/2565._8 , 1859._8/4104._8 , -11._8/40._8 /),(/5,5/)) )

		c_mat =  (/ 1._8/4._8, 3._8/8._8, 12._8/13._8, 1._8, 1._8/2._8 /)

		b_rk4 = (/25._8/216._8, 0._8, 1408._8/2565._8, 2197._8/4104._8, -1._8/5._8, 0._8/)
		b_rk5 = (/16._8/135._8, 0._8, 6656._8/12825._8, 28561._8/56430._8, -9._8/50._8, 2._8/55._8 /)

		cumu_solve_time = 0._8
		cumu_other_time = 0._8
		cumu_sum_time = 0._8

		WHILE_DO: DO

			IF ( ( st%t - t0 > sys%print_delay ) &
				.or. (st%t > sys%tmax)  &
				.or. (info > 1) ) EXIT WHILE_DO

			dt = sys%dt * ( 1._8 /3._8 )**slow_fact_exp

			if ( abs(st%t-tref)<sys%dt/100 ) then
				dt = dt/1000
			elseif ( abs(st%t-tref)<sys%dt/10 ) then
				dt = dt/100
			elseif ( abs(st%t-tref)<3*sys%dt ) then
				dt = dt/10
			end if

			!===============
			!== RK 45
			!===============

			!== k1 calculation  ==========
			if ( .not. dt_corrected ) then
				CALL calc_derivatives( sys, st, kp1, kf1, t_solve, t_other )
				cumu_solve_time = cumu_solve_time + t_solve
				cumu_other_time = cumu_other_time + t_other
			end if

			!== k2 calculation  ==========
			midst%y = st%y + a_mat(1,1)*dt*kf1
			midst%p = st%p + a_mat(1,1)*dt*kp1
			midst%t = st%t + c_mat(1)*dt
			CALL update_sums( sys, midst, t_sum )
			CALL calc_derivatives( sys, midst, kp2, kf2, t_solve, t_other  )
			cumu_solve_time = cumu_solve_time + t_solve
			cumu_other_time = cumu_other_time + t_other
			cumu_sum_time = cumu_sum_time + t_sum

			!== k3 calculation  ==========
			midst%y = st%y + dt*( a_mat(2,1)*kf1 + a_mat(2,2)*kf2 )
			midst%p = st%p + dt*( a_mat(2,1)*kp1 + a_mat(2,2)*kp2 )
			midst%t = st%t + c_mat(2)*dt
			CALL update_sums( sys, midst, t_sum )
			CALL calc_derivatives( sys, midst, kp3, kf3, t_solve, t_other  )
			cumu_solve_time = cumu_solve_time + t_solve
			cumu_other_time = cumu_other_time + t_other
			cumu_sum_time = cumu_sum_time + t_sum

			!== k4 calculation  ==========
			midst%y = st%y + dt*( a_mat(3,1)*kf1 + a_mat(3,2)*kf2 + a_mat(3,3)*kf3 )
			midst%p = st%p + dt*( a_mat(3,1)*kp1 + a_mat(3,2)*kp2 + a_mat(3,3)*kp3 )
			midst%t = st%t + c_mat(3)*dt
			CALL update_sums( sys, midst, t_sum )
			CALL calc_derivatives( sys, midst, kp4, kf4, t_solve, t_other  )
			cumu_solve_time = cumu_solve_time + t_solve
			cumu_other_time = cumu_other_time + t_other
			cumu_sum_time = cumu_sum_time + t_sum

			!== k5 calculation  ==========
			midst%y = st%y + dt*( a_mat(4,1)*kf1 + a_mat(4,2)*kf2 + a_mat(4,3)*kf3 + a_mat(4,4)*kf4 )
			midst%p = st%p + dt*( a_mat(4,1)*kp1 + a_mat(4,2)*kp2 + a_mat(4,3)*kp3 + a_mat(4,4)*kp4 )
			midst%t = st%t + c_mat(4)*dt
			CALL update_sums( sys, midst, t_sum )
			CALL calc_derivatives( sys, midst, kp5, kf5, t_solve, t_other  )
			cumu_solve_time = cumu_solve_time + t_solve
			cumu_other_time = cumu_other_time + t_other
			cumu_sum_time = cumu_sum_time + t_sum

			!== k6 calculation  ==========
			midst%y = st%y + dt*( a_mat(5,1)*kf1 + a_mat(5,2)*kf2 + a_mat(5,3)*kf3 + a_mat(5,4)*kf4 + a_mat(5,5)*kf5 )
			midst%p = st%p + dt*( a_mat(5,1)*kp1 + a_mat(5,2)*kp2 + a_mat(5,3)*kp3 + a_mat(5,4)*kp4 + a_mat(5,5)*kp5 )
			midst%t = st%t + c_mat(5)*dt
			CALL update_sums( sys, midst, t_sum )
			CALL calc_derivatives( sys, midst, kp6, kf6, t_solve, t_other  )
			cumu_solve_time = cumu_solve_time + t_solve
			cumu_other_time = cumu_other_time + t_other
			cumu_sum_time = cumu_sum_time + t_sum

			pdot_RK4 = ( b_rk4(1)*kp1 + b_rk4(2)*kp2 + b_rk4(3)*kp3 + b_rk4(4)*kp4 + b_rk4(5)*kp5 + b_rk4(6)*kp6 )
			pdot_RK5 = ( b_rk5(1)*kp1 + b_rk5(2)*kp2 + b_rk5(3)*kp3 + b_rk5(4)*kp4 + b_rk5(5)*kp5 + b_rk5(6)*kp6 )
			fdot_RK4 = ( b_rk4(1)*kf1 + b_rk4(2)*kf2 + b_rk4(3)*kf3 + b_rk4(4)*kf4 + b_rk4(5)*kf5 + b_rk4(6)*kf6 )
			fdot_RK5 = ( b_rk5(1)*kf1 + b_rk5(2)*kf2 + b_rk5(3)*kf3 + b_rk5(4)*kf4 + b_rk5(5)*kf5 + b_rk5(6)*kf6 )

			error_rk5 = dt*( sum( abs(fdot_RK5-fdot_RK4) ) + sum( abs(pdot_RK5-pdot_RK4) ) )

			if (  error_rk5 < sys%err_lim  ) then

				st%p = st%p + dt*pdot_RK5
				st%y = st%y + dt*fdot_RK5
				st%t = st%t + dt
				CALL update_sums( sys, st, t_sum )
				cumu_sum_time = cumu_sum_time + t_sum

				if (slow_fact_exp < min_slow_fact_exp) then
					min_slow_fact_exp = slow_fact_exp
				end if
				if (slow_fact_exp > max_slow_fact_exp) then
					max_slow_fact_exp = slow_fact_exp
				end if


				if ( (slow_fact_exp > 0) .and. (st%t-t_corrected > 50*dt) ) then
					if (slow_fact_exp >= 2) then
						slow_fact_exp = slow_fact_exp - 2
					else
						slow_fact_exp = slow_fact_exp - 1
					endif
				end if
				dt_corrected = .false.

			else 

				t_corrected = st%t
				slow_fact_exp = slow_fact_exp + 1

				if (bt_allowed .eqv. .true.) then
					if ( ( ( slow_fact_exp > sys%lim_slow_fact_exp ) .and. (st%ncs<sys%ncs_max) .and. (st%ncs>sys%ncs_ini) ) &
						.or. ( ( slow_fact_exp > sys%lim_slow_fact_exp+2 ) .and. (st%ncs==sys%ncs_max) ) &
						.or. ( ( slow_fact_exp > sys%lim_slow_fact_exp+4 ) .and. (st%ncs==sys%ncs_ini) ) ) then

						CALL print_log( sys, '=========================================', sys%logfile )
						CALL print_log( sys, '== ERROR -1: value of slow_fact_exp higher than sys%lim_slow_fact_exp', sys%logfile)
						CALL print_log( sys, '== BACKTRACKING ==', sys%logfile)
						CALL print_log( sys, '=========================================', sys%logfile)
						
						!-- specifying info=2 means program should exit routine and backtrack
						info = 2

					end if
				else
					if ( slow_fact_exp > 20) then
						CALL print_log( sys, '=========================================', sys%logfile )
						CALL print_log( sys, '== ERROR -1: value of slow_fact_exp higher than 20', sys%logfile)
						CALL print_log( sys, '== ABORT', sys%logfile)
					end if
				end if

				dt_corrected = .true.

			end if

		END DO WHILE_DO

		RETURN
	END SUBROUTINE time_evolve_RK45_no_buffer

	FUNCTION average_of_five_median(array)
		real(8), intent(in)   				::	array(:)
		real(8), dimension(size(array,1))   ::	sorted_array, array_copy
		integer, dimension(size(array,1))   ::	minloc_array
		real(8)								::  average_of_five_median
		integer								::  min_ind, mid_ind, i, L, avrg_L_2

		L = size(array)

		array_copy = array
		sorted_array = 0._8
		mid_ind = int(L/2)
		avrg_L_2 = int((L/20)/2)

		do i=1, L
			minloc_array = minloc(array_copy)
			sorted_array(i) = array_copy( minloc_array( 1 ) )
			array_copy(minloc_array( 1 ))  = 1.e9_8
		end do

		average_of_five_median = sum(sorted_array( mid_ind-avrg_L_2:mid_ind+avrg_L_2 )) / (2*avrg_L_2+1)
	END FUNCTION

	SUBROUTINE calc_derivatives( sys, st, pdot, fdot, solve_time, other_time  )
		type(param), intent(in)       :: 	sys
		type(state), intent(in)       :: 	st
		real(8)					      :: t, tA1, tA2, tB1, tB2
		real(8), intent(out)	      :: solve_time, other_time
		complex(cx), intent(out)   :: pdot(size(st%p,1),size(st%p,2))
		complex(cx), intent(out)   :: fdot(size(st%y,1),size(st%y,2),size(st%y,3))
		integer                   :: nmodes, nq, ncs
		complex(cx), dimension(size(st%y,2))             ::  bigP, dE_dpc_sj_ST!, p, pc,
		complex(cx), dimension( size(st%y,1), size(st%y,2) )            ::  p, pc
		complex(cx), dimension(size(st%y,2), size(st%y,3))     :: f, bigF, inv_rovm_mul_F_m_fnP
		complex(cx), dimension(size(st%y,1), size(st%y,2), size(st%y,3))     :: dE_dyc_sj_ST
		complex(cx), dimension(size(st%y,2),size(st%y,2))         :: inv_rovm, rovm, b, Amat, RHS, d
		complex(cx), dimension(size(st%y,2),size(st%y,2),size(st%y,2))  :: alphaT
		complex(cx), allocatable          :: packed_RHS(:), d_packed(:)
		complex(cx), allocatable  :: mat2D(:,:)
		integer                                :: info,i,j,k,n,m,s,ii,jj
		complex(cx), dimension(size(st%y,2))             ::  summed_d


		allocate( mat2D( size(st%y,2)**2,size(st%y,2)**2 ) )
		allocate( packed_RHS( size(st%y,2)**2 ) )
		allocate( d_packed( size(st%y,2)**2 ) )

		p = st%p
		pc = conjg(st%p)

		nq = size( st%y, 1 )
		ncs = size( st%y, 2 )
		nmodes= size( st%y, 3 )

		t = st%t
		solve_time = 0._8
		call CPU_TIME( tB1 )
		dE_dyc_sj_ST = 0._8
		do j=1, ncs
			do n=1, ncs
				do s=1, nq
					dE_dyc_sj_ST(s,j,:) = dE_dyc_sj_ST(s,j,:) + pc(s,j)*p(s,n)*st%ovm(s,j,s,n)*(  &
						+ sys%wwk(:)*st%y(s,n,:) &
						+ sys%Ad*dcos( sys%wd*st%t )*sys%Omat(1,:)  &
						+ (st%y(s,n,:)-0.5_8*st%y(s,j,:))*( sys%w_qb(s)+st%bigW(s,j,n) + sys%Ad*dcos(sys%wd*st%t)*st%bigU(s,j,n) ) ) &
						- 0.5_8*pc(s,n)*p(s,j)*st%ovm(s,n,s,j)*st%y(s,j,:)*( sys%w_qb(s) &
						+ st%bigW(s,n,j) &
						+ sys%Ad*dcos(sys%wd*st%t)*st%bigU(s,n,j) &
						)
				end do

				!do s=1, nq
				!	do i=1, size( p,1 )
				do ii=1, sys%nb_el_to_keep
					s = sys%ind( ii, 1 )
					i = sys%ind( ii, 2 )

					dE_dyc_sj_ST(s,j,:) = dE_dyc_sj_ST(s,j,:) + pc(s,j)*p(i,n)*st%ovm(s,j,i,n)*( &
						!+ sys%nij(s,i)*sys%g_qc*sys%Omat(1,:)*( 1 + ( conjg(st%y0(s,j)) + st%y0(i,n) )*(sys%expnt-1) ) &
						+ sys%nij(s,i)*sys%gk(:)*( 1 + ( conjg(st%y0(s,j)) + st%y0(i,n) )*(sys%expnt-1) ) &
						+ ( st%y(i,n,:)-0.5_8*st%y(s,j,:) )*st%bigL(s,j,i,n) ) &
						- 0.5_8*pc(i,n)*p(s,j)*st%ovm(i,n,s,j) * st%y(s,j,:)*st%bigL(i,n,s,j)

				end do

			end do
		end do

		call CPU_TIME( tB2 )
		other_time = tB2 - tB1

		do s=1,nq

			mat2D = 0._8
			packed_RHS = 0._8
			d_packed = 0._8

			f = st%y(s,:,:)

			rovm=st%ovm(s,:,s,:)
			inv_rovm=rovm

			CALL invertH(ncs,inv_rovm,info)

			do n=1,ncs
				dE_dpc_sj_ST( n ) = dE_dpc_sj( sys, st, s, n )
			end do

			bigP = -Ic*dE_dpc_sj_ST
			do k=1, nmodes
				bigF(:,k) = -Ic*( dE_dyc_sj_ST(s,:,k)/conjg(p(s,:))&
					+0.5_8*( dE_dpc_sj_ST + conjg(dE_dpc_sj_ST)*p(s,:)/conjg(p(s,:)) )*f(:,k) )
			end do
			inv_rovm_mul_F_m_fnP = 0._8
			do k=1, nmodes
				inv_rovm_mul_F_m_fnP(:,k) = matmul( inv_rovm(:,:),bigF(:,k) )&
					- f(:,k)*matmul( inv_rovm(:,:),bigP(:))
			end do
			Amat = matmul( conjg(f(:,:)), TRANSPOSE(inv_rovm_mul_F_m_fnP) )
			b(:,:) = matmul( conjg(f(:,:)),TRANSPOSE(f(:,:)) )

			!-- build alphaTensor
			do n=1, ncs
				do m=1, ncs
					alphaT(:,n,m)=matmul( inv_rovm(:,:) , rovm(:,n)*(b(:,m)-b(:,n)) )
				end do
			end do

			!-- build the right-hand side
			RHS = matmul(inv_rovm, rovm * Amat)

			!-- build the d matrix
			do i=1, ncs
				do n=1, ncs
					packed_RHS((n-1)*ncs+i)=RHS(i,n)
					mat2D((n-1)*ncs+i,(n-1)*ncs+i) = 1.0_8
					do m=1, ncs
						mat2D((n-1)*ncs+i,(m-1)*ncs+n) = &
							mat2D((n-1)*ncs+i,(m-1)*ncs+n) + alphaT(i,n,m)
					end do
				end do
			end do

			call CPU_TIME( tA1 )
			CALL solveEq_c( mat2D , packed_RHS, size(mat2D,1), d_packed )
			call CPU_TIME( tA2 )
			solve_time = solve_time + tA2 - tA1

			do i=1, ncs
				do n=1, ncs
					d(i,n) = d_packed((n-1)*ncs+i)
				end do
			end do

			summed_d = sum( d, dim=2 )
			do n=1, ncs
				fdot(s,n,:) = ( inv_rovm_mul_F_m_fnP(n,:) - matmul( TRANSPOSE(f), d(n,:) ) + summed_d(n)*f(n,:) ) / p(s,n)
			end do
			pdot(s,:) = matmul( inv_rovm(:,:), bigP(:) ) - summed_d(:)  &
				+ 0.5_8*p(s,:)*( sum( fdot(s,:,:)*conjg(f(:,:)) + conjg(fdot(s,:,:))*f(:,:), dim=2 ) )

		end do

		return
	END SUBROUTINE

	FUNCTION dE_dpc_sj( sys, st, s, j )
		type(param), intent(in)       :: 	sys
		type(state), intent(in)       :: 	st
		integer,intent(in)            ::  s,j     
		complex(cx) 	  				  ::  dE_dpc_sj

		dE_dpc_sj = sum( st%p(s,:)*st%ovm(s,j,s,:)*( sys%w_qb(s) + st%bigW(s,j,:) &
			+ sys%Ad*dcos( sys%wd * st%t ) * st%bigU(s,j,:)  )) &
			+ sum( sum(st%p(:,:)*st%ovm(s,j,:,:)*st%bigL(s,j,:,:), dim=2) )
	END FUNCTION

	SUBROUTINE calc_wigner( p_in, f_in, ovm, xmin, xmax,  wigner, xnum )
		complex(cx), intent(in)                          ::  p_in(:,:)
		complex(cx), intent(in)                          ::  f_in(:,:)
		integer, intent(in)							    ::  xnum
		real(8), intent(in)							    ::  xmin, xmax
		complex(cx), intent(in)                          ::  ovm(:,:,:,:)
		real(8)							                ::  dx, x, p
		complex(cx)					  					::  tmp, zz
		integer							  				::  i,n,m,xi,xj,ncs,nl
		real(8), intent(out)				            ::  wigner(xnum,xnum)

		nl = size( p_in, 1 )
		ncs = size( p_in, 2 )
		dx =   (xmax - xmin) / dble(xnum-1) 

		do xi=1, xnum
			do xj=1, xnum

				tmp = 0._8

				x = xmin + dx*(xi-1)
				p = xmin + dx*(xj-1)

				do i=1, nl
					do n=1, ncs
						do m=1, ncs

							zz =  x + Ic*p 

							tmp = tmp + conjg(p_in(i,n))*p_in(i,m) &
								* exp( -2._8 * ( real(zz)+Ic*aimag(zz) - f_in(i,m) ) &
								* ( real(zz)-Ic*aimag(zz) - conjg( f_in(i,n) ) ) ) &
								* ovm(i,n,i,m)

						end do
					end do
				end do

				wigner( xi, xj ) = (2._8/pi)*tmp

			end do
		end do

		return
	END SUBROUTINE

	!-- routine to coherent states 
	SUBROUTINE add_coherent_state(sys,st)
		type(param), intent(in)				:: sys
		type(state), intent(in out)   		:: st
		type(state)						    :: st_add, ba_st
		integer								:: n, s, n_heaviest
		character(len=52)					:: log_message
		real(8)								:: highest_weight, dummy

		ba_st = st 

		!-- Re-allocate st with more polarons
		CALL allocate_state(sys,st,ba_st%ncs+1)
		CALL allocate_state(sys,st_add,1)

		st_add%y(:,:,:) = 0._cx
		do n=1, ba_st%ncs
			st%y(:,n,:) = ba_st%y(:,n,:)
			st%p(:,n) = ba_st%p(:,n)
		end do

		highest_weight = -1
		do s=1, sys%nl
			do n=1, ba_st%ncs
				if ( abs( ba_st%p(s,n) ) > highest_weight ) then
					highest_weight = abs( ba_st%p(s,n) )
					n_heaviest = n
				end if
			end do
				
			st%y(s,ba_st%ncs+1,:) = sys%adding_ratio * ba_st%y(s,n_heaviest,:)

		end do
		
		st%p(:,ba_st%ncs+1) = sys%p0
		st%t = ba_st%t

		CALL update_sums( sys, st, dummy )
		CALL normalise(st)

		write(log_message,'(a15,I3,a4,I3,a6,f7.2)') "-- CS ADDED: from ",st%ncs-1," to ",st%ncs, " at t=", st%t
		CALL print_log( sys, log_message, sys%logfile )
	END SUBROUTINE
	
	SUBROUTINE print_evolution_data(sys,tr)
		type(param), intent(in)  		::  sys
		type(traj), intent(in)			::  tr
		character(len=300)		  		::  name_lvls, name_ErEN, name_p0, name_p1, name_n, name_y0_0,&
												name_y0_1, name_entr, name_ovm_eval, name_apadag
		real(rl)						::  t
		integer							::  i, n, s, k

		name_lvls="data/PPLT_"//trim(adjustl(parameterchar(sys)))//".d"
		name_ErEN="data/NEr_"//trim(adjustl(parameterchar(sys)))//".d"
		name_n="data/PHOTONS_"//trim(adjustl(parameterchar(sys)))//".d"
		name_ovm_eval="data/ovm_eval_"//trim(adjustl(parameterchar(sys)))//".d"
		name_apadag ="data/APADAG_"//trim(adjustl(parameterchar(sys)))//".d"

		open (unit=10,file=name_lvls,action="write",status="replace")
		open (unit=11,file= name_ErEN,action="write",status="replace")
		open (unit=13,file= name_n,action="write",status="replace")
		open (unit=14,file= name_ovm_eval,action="write",status="replace")
		open (unit=15,file= name_apadag,action="write",status="replace")

		do i=1, tr%i_time

			t = tr%time_ar(i)
			write(10,'(f25.15)',advance='no') t
			do s=1,sys%nl
				write(10,'(2f25.15)',advance='no') tr%populations_ar(i,s,2)
			end do
			do s=1,sys%nl
				write(10,'(2f25.15)',advance='no') tr%populations_ar(i,s,1)
			end do
			write(10,*)

			write(11,'(4f25.15)') t, tr%norm_ar(i),  tr%error_ar(i,1), tr%error_ar(i,2)
			write(13,'(f25.15)',advance='no') t
			do k=1, sys%nmodes
				write(13,'(f25.15)',advance='no') tr%n_ar(i,k)
			end do
			write(13,*)

			write(14,'(f25.15)',advance='no') t
			do n=1, sys%ncs_max
				write(14,'(f25.15)',advance='no') tr%ovm_eval_0_ar(i,n)
			end do
			do n=1, sys%ncs_max
				write(14,'(f25.15)',advance='no') tr%ovm_eval_1_ar(i,n)
			end do
			write(14,*)

			write(15,'(2f25.15)') t, tr%apadag_ar(i)

		end do
		close(10)
		close(11)
		close(13)
		close(14)
		close(15)
	END SUBROUTINE

	SUBROUTINE print_qusimi(sys,st,res)
		type(state), intent(in) 		    			::  st
		type(param), intent(in) 		    			::  sys
		integer, intent(in)								::  res
		real(rl)					  					::  Q
		integer							  				::  i,j,n,m,s,xnum
		real(rl)										::  xmin,xmax,pmin,pmax,x,p
		complex(cx)					  					::  Q_tmp, zz
		complex(cx), dimension( st%ncs )				::  y0
		character(len=300)		  						::  filename

		xmin = -10
		xmax = 10
		pmin = -10
		pmax = 10
		xnum = int( (xmax - xmin) * res )
		y0 = 0._cx


		filename="data/Q_"//trim(adjustl(parameterchar(sys)))//".d"
		open( unit=10, file=filename, action="write", status="replace" )

		do j=0,xnum
			do i=0,xnum
				Q_tmp = 0._cx
				do s=1, size( st%p,1 )

					x = xmin + (xmax-xmin)*i/dble(xnum)
					p = pmin + (pmax-pmin)*j/dble(xnum)

					do n=1, st%ncs
						y0(n) = sum( st%y(s,n,:)*sys%Omat(1,:) )
					end do

					do n=1,st%ncs
						do m=1,st%ncs
							zz =  x + Ic*p
							Q_tmp = Q_tmp + conjg(st%p(s,n))*st%p(s,m)* ov_scalar( zz, y0(m) ) * ov_scalar( y0(n) , zz )
						end do
					end do

				end do
				Q = (2._rl/pi) * real( Q_tmp )

				write(10,'(f25.10)', advance='no') Q

			end do
			write(10,*) 

		end do

		close(10)
	END SUBROUTINE

	!== printing the state vector
	SUBROUTINE print_yks_DB(sys,st)
		type(param), intent(in)								::  sys
		type(state), intent(in)								::  st
		integer												::  i, n, k
		character(len=300)									::  name_yks

		name_yks="data/YKS_DB_"//trim(adjustl(parameterchar(sys)))//".d"
		open (unit=100,file= name_yks,action="write",status="replace")

		do k = 1, sys%nmodes
			write(100,'(f25.15)',advance='no') sys%wk(k)
			do i=1,sys%nl
				do n=1,st%ncs
					write(100,'(f25.15)', advance='no') real( st%y(i,n,k) )
				end do
			end do
			do i=1,sys%nl
				do n=1,st%ncs
					write(100,'(f25.15)', advance='no') aimag( st%y(i,n,k) )
				end do
			end do
			write(100,*)
		end do
		close(100)
	END SUBROUTINE

	SUBROUTINE print_yks(sys,st)
		type(param), intent(in)								::  sys
		type(state), intent(in)								::  st
		integer												::  i, n, k
		complex(cx), dimension(sys%nl,st%ncs,sys%nmodes) 	::  ynk
		character(len=300)									::  name_yks

		name_yks="data/YKS_"//trim(adjustl(parameterchar(sys)))//".d"
		open (unit=100,file= name_yks,action="write",status="replace")

		do i=1, sys%nl
			do n=1, st%ncs
				do k=1, sys%nmodes
					ynk(i,n,k) = sum( st%y(i,n,:)*sys%Omat(k,:) )
				end do
			end do
		end do

		do k = 1, sys%nmodes
			write(100,'(f25.15)',advance='no') sys%wk(k)
			do i=1,sys%nl
				do n=1,st%ncs
					write(100,'(f25.15)', advance='no') real( ynk(i,n,k) )
				end do
			end do
			do i=1,sys%nl
				do n=1,st%ncs
					write(100,'(f25.15)', advance='no') aimag( ynk(i,n,k) )
				end do
			end do
			write(100,*)
		end do
		close(100)
	END SUBROUTINE

	SUBROUTINE print_omat(sys)
		type(param), intent(in)								::  sys
		integer												::  i, j
		character(len=300)									::  name_omat

		name_omat="data/OMAT_"//trim(adjustl(parameterchar(sys)))//".d"
		open (unit=100,file= name_omat,action="write",status="replace")

		do i=1, sys%nmodes
			do j=1, sys%nmodes
				write(100,'(f25.15)', advance='no') sys%omat( i,j )
			end do
			write(100,*) 
		end do
		write(100,*)
		close(100)
	END SUBROUTINE

	SUBROUTINE print_ps(sys, st)
		type(param), intent(in)			::  sys
		type(state), intent(in)			::  st
		integer							::  i, n
		character(len=300)				::  name_ps

		name_ps="data/PS_"//trim(adjustl(parameterchar(sys)))//".d"
		open (unit=100,file= name_ps,action="write",status="replace")

		do i=1,sys%nl
			do n=1,st%ncs
				write(100,'(f25.15)',advance='no') real( st%p(i,n) )
			end do
		end do
		do i=1,sys%nl
			do n=1,st%ncs
				write(100,'(f25.15)',advance='no') aimag( st%p(i,n) )
			end do
		end do
		write(100,*)

		close(100)
	END SUBROUTINE

	FUNCTION error( sys, st )
		type(param), intent(in)   ::     sys
		type(state), intent(in)   ::     st
		type(param)			      ::     sys2
		type(state)			      ::     ost, new_st
		real(8)				      ::	 energy_t
		real(8)  			 	  ::	 error
		complex(cx)				  ::	 tmp1, tmp2, tmp3, tmp4, tmp_sum, sum_diag_g2mat
		real(8)				  	  ::	 at, dummy1, dummy2, dummy3
		integer                   ::	 nmodes, nl, ncs,m,l,n,i,j,k,pp, msfexp, sfexp, s
		complex(cx), dimension(sys%nl,st%ncs)         		::	 p, pdot, pc, pdotc, new_pdot
		complex(cx), dimension(sys%nl,st%ncs,sys%nmodes)     ::	 y, ydot, yc, ydotc, new_ydot
		complex(cx), dimension(st%ncs)         				::	 opdd
		complex(cx), dimension(st%ncs,sys%nmodes)     		::	 oydd
		complex(cx), dimension(st%ncs,st%ncs)         		::	 ovmr
		complex(cx), dimension(sys%nl,st%ncs,sys%nl,st%ncs) 		::	 kap
		integer									::   info

		tmp1 = 0._8
		tmp2 = 0._8
		tmp3 = 0._8
		tmp4 = 0._8

		msfexp = 0
		sfexp = 0

		pdot = 0._rl
		ydot = 0._rl
		new_pdot = 0._rl
		new_ydot = 0._rl

		sys2 = sys
		sys2%dt = sys%dt*1e-4
		ost = st
		CALL calc_derivatives( sys2, ost, pdot, ydot, dummy1, dummy2 )

		new_st = st
		sys2%tmax = new_st%t + 2*sys2%dt
		CALL time_evolve_RK45_no_buffer( sys2, new_st, 1.0e3_8, 1.0e3_8, msfexp, msfexp, sfexp, info, .false., &
				dummy1, dummy2, dummy3 )
		CALL calc_derivatives( sys2, new_st, new_pdot, new_ydot, dummy1, dummy2 )

		energy_t = energy( sys, ost )

		p = ost%p
		y = ost%y

		pc = conjg( p )
		yc = conjg( y )
		pdotc = conjg( pdot )
		ydotc = conjg( ydot )

		at = sys%ad*dcos( sys%wd*ost%t )

		nl = size( p,1 )
		ncs = size( p,2 )
		nmodes = size( y,3 )

		do s=1, nl
			do m=1, ncs
				do l=1, nl
					do n=1, ncs
						kap( s,m,l,n ) = sum( ydot(s,m,:)*yc(s,m,:) &
							+ ydotc(s,m,:)*y(s,m,:) - 2._8*yc(l,n,:)*ydot(s,m,:) )
					end do
				end do
			end do
		end do

		LOOPi: do i=1, nl

			ovmr = ost%ovm(i,:,i,:)

			opdd = ( new_pdot(i,:) - pdot(i,:) )/(new_st%t-ost%t)
			oydd = ( new_ydot(i,:,:) - ydot(i,:,:) )/(new_st%t-ost%t)

			do m=1, ncs
				do n=1, ncs
					!==== tmp1 cajcujation
					tmp1 = tmp1 + ovmr(m,n)*( &
						+ pdotc(i,m)*pdot(i,n) &
						- 0.5_8 * pdotc(i,m)*p(i,n)*kap(i,n,i,m) &
						- 0.5_8 * pc(i,m)*pdot(i,n)*conjg(kap(i,m,i,n)) &
						+ pc(i,m)*p(i,n)*( sum( ydotc(i,m,:)*ydot(i,n,:) )&
						+ 0.25_8*conjg(kap(i,m,i,n))*kap(i,n,i,m)&
						)&
						)

					!==== tmp4 cajcujation
					tmp4 = tmp4 + pc(i,m)*ovmr(m,n)*( &
						+ opdd(n) &
						- pdot(i,n)*kap(i,n,i,m) &
						+ p(i,n)*( sum( yc(i,m,:)*oydd(n,:)&
						- 0.5_8*( y(i,n,:)*conjg(oydd(n,:))&
						+ yc(i,n,:)*oydd(n,:)&
						+ 2._8*ydotc(i,n,:)*ydot(i,n,:) ) )&
						+ 0.25_8*kap(i,n,i,m)**2 )&
						)

					!==== tmp2 cajcujation
					tmp2 = tmp2 + pdotc(i,m)*p(i,n)*ost%ovm(i,m,i,n)*( &
						+ sys%w_qb(i) &
						+ ost%bigW(i,m,n)  &
						+ at*ost%bigU(i,m,n) ) &
						+ pc(i,m)*p(i,n)*ovmr(m,n)*( &
						- 0.5_8*conjg(kap(i,m,i,n)) &
						*( sys%w_qb(i)+ost%bigW(i,m,n)+at*ost%bigU(i,m,n) ) &
						+ sum( sys%wwk(:)*ydotc(i,m,:)*y(i,n,:) &
						+ at*sys%Omat(1,:)*ydotc(i,m,:) ) )

					do j=1, nl
						tmp2 = tmp2 + pdotc(i,m)*p(j,n)*ost%ovm(i,m,j,n)*ost%bigL(i,m,j,n) &
							+ pc(i,m)*p(j,n)*ost%ovm(i,m,j,n)*( &
							- 0.5_8*conjg(kap(i,m,j,n))*ost%bigL(i,m,j,n)  &
							!+ sys%nij(i,j)*sum( sys%g_qc*sys%Omat(1,:)*ydotc(i,m,:) ) )
							+ sys%nij(i,j)*sum( sys%gk(:)*ydotc(i,m,:) ) )
					end do

					!==== tmp3 calculation NEW
					tmp3 = tmp3 + pc(i,m)*p(i,n)*ost%ovm(i,m,i,n)*( (sys%w_qb(i) &
						+ ost%bigW(i,m,n) &
						+ at*ost%bigU(i,m,n) )**2 &
						+ sum( sys%wwk(:)**2*yc(i,m,:)*y(i,n,:)  &
								+ sys%wwk(:)*at*sys%Omat(1,:)*(yc(i,m,:)+y(i,n,:)) ) &
						+ at**2 ) 

					!-- ================================================
					!-- ================================================
					!-- g_qc terms without omat need to be updated for LF modes
					!-- ================================================
					!-- ================================================
					do j=1, nl
						tmp3 = tmp3 + pc(i,m)*p(j,n)*ost%ovm(i,m,j,n)*( &
							!+ sys%g_qc**2*sum( sys%nij(:,j)*sys%nij(:,i) ) &
							!+ 2._8*at*sys%g_qc*sys%nij(i,j) &
							!+ ( ( conjg(ost%y0(i,m)) + ost%y0(j,n) )**2 ) * ( 2._8*at*sys%g_qc*sys%nij(i,j) ) &
							+0& !( ( conjg(ost%y0(i,m)) + ost%y0(j,n) )**2 ) * (  sys%g_qc**2*sum( sys%nij(:,j)*sys%nij(:,i) ) ) &
							) 
					end do

					!!-- updated code fix attempt
					!do j=1, nl
					!	tmp3 = tmp3 + pc(i,m)*p(j,n)*ost%ovm(i,m,j,n)*( &
					!		+ sys%n2ij(i,j)*( 1 & !sys%sum_g2  &
					!		 			+ sys%g_qc**2*( conjg(ost%y0(i,m)) + ost%y0(j,n) )**2 ) &
					!		+ 2._8*at*( sys%nij(i,j)*(1+sys%sum_og*0) & 
					!					+ ost%bigL(i,m,j,n)*( conjg(ost%y0(i,m)) + ost%y0(j,n) ) ) &
					!	)
					!end do

					do j=1, nl
						tmp3 = tmp3 + pc(i,m)*p(j,n)*ost%ovm(i,m,j,n)*( &
							+ ( sys%w_qb(i)+sys%w_qb(j) )*ost%bigL(i,m,j,n) &
							+ 2._8*sum( sys%wwk*yc(i,m,:)*y(j,n,:) )*ost%bigL(i,m,j,n) &
							+ sys%nij(i,j)*sum( sys%gk(:)*sys%wwk(:)*( y(j,n,:)+yc(i,m,:) ) ) &
							)
					end do

				end do
			end do

		end do LOOPi

		error = abs( ( -0.5*real(tmp4) + 0.5*tmp1 - 2*aimag(tmp2) + tmp3 )  )
	END FUNCTION error

	FUNCTION median( arr )
		real(rl), intent(in)	::  arr(:)
		real(rl)				::  arr2( size(arr,1) )
		real(rl)				::  median
		integer					::  i

		arr2 = arr
		do i=1, int( size(arr,1)/2 )
			arr2( minloc(arr2) ) = 1e12
		end do

		median =  minval(arr2)
	END FUNCTION

	SUBROUTINE SolveEq_c(A,B,N,res)
		complex(cx), dimension(N,N), intent(in)   ::  A
		complex(cx), dimension(N), intent(in)	 ::  B
		INTEGER, intent(in)                	     ::  N
		complex(cx), intent(out)					 ::  res(size(B,1))
		INTEGER                	   		    	 ::  INFO,LDA,LDB,NRHS
		INTEGER, dimension(size(A,1))		  	 ::  IPIV   !-- pivot indices

		NRHS = 1  						!-- number of right hand sides
		LDA = size(A,1)				!-- the leading dimension of A multiply_c
		LDB = size(B,1)				!-- the leading dimension of B
		info=0							!-- 0 is successful

		res = B

		CALL ZGESV(N,NRHS,A,LDA,IPIV,res,LDB,INFO)  !-- Solve by performing the Bunch-Kaufman factorisation

		if (info /= 0) then
			print*, "Failure in DGESV - solving real system of equations"
			print*,"INFO = ",info
		end if
	END SUBROUTINE

	SUBROUTINE InvertH(N,A,info)
		integer, intent(in)						::  N
		complex(cx), intent(in out)              ::  A(N,N)
		integer, intent(out)					::  info
		INTEGER                	   		        ::  LDA,LWORK,i,j
		INTEGER, dimension(size(A,1))		    ::  IPIV
		complex(cx), allocatable 				::  WORK(:)

		LDA = N
		LWORK = N
		info=0
		allocate(WORK(LWORK))

		CALL ZHETRF('U',N,A,LDA,IPIV,WORK,LWORK,INFO)  !-- Performs the Bunch-Kaufman factorisation

		if (info==0) then

			CALL ZHETRI('U',N,A,LDA,IPIV,WORK,INFO) !-- CAREFUL: returns only triangular part
			do i=1,N
				do j=1,N
					if (i>j) then
						a(i,j) = conjg(a(j,i))
					end if
				end do
			end do

			if (info /= 0) then
				print*, "Failure in the inversion step, ZHETRI"
				print*,"info=", info
			end if

		else
			print*, "Failure in ZHETRF, Bunch-Kaufman factorisation"
			print*,"info=", info
		end if
	END SUBROUTINE InvertH

	FUNCTION KD(a,b)
		integer,intent(in)   ::  a,b
		integer     ::  KD

		if (a==b) then
			KD = 1
		else
			KD = 0
		end if
	END FUNCTION

	FUNCTION all_p2_above_threshold( sys, st )
		type(param), intent(in)     ::	sys
		type(state), intent(in)     ::	st
		logical						::  all_p2_above_threshold
		integer						::  n
		
		all_p2_above_threshold = .true.
		do n=1, size( st%p,2 )

			if ( abs( st%p(sys%qb_ini+1,n) )**2 < sys%p2_threshold ) then
				all_p2_above_threshold = .false.
			end if

		end do
	END FUNCTION

END MODULE dynamics

!	FUNCTION dE_dyc_sjk( sys, st, s, j )
!		type(param), intent(in)       :: 	sys
!		type(state), intent(in)       :: 	st
!		integer,intent(in)              			 			 ::  s,j
!		complex(cx),dimension( sys%nmodes ) 	               ::  dE_dyc_sjk
!		complex(cx), dimension( sys%nl, st%ncs )            ::  p, pc
!		integer ::  i,n,ii
!
!		p = st%p
!		pc = conjg(st%p)
!		dE_dyc_sjk = 0._8
!
!		!dE_dyc_sjk  = sys%wwk(:)*matmul( transpose(st%y(s,:,:)), pc(s,j)*p(s,:)*st%ovm(s,j,s,:) ) &
!		!	+ sys%Ad*dcos( sys%wd*st%t )*sys%Omat(1,:) * sum( pc(s,j)*p(s,:)*st%ovm(s,j,s,:) ) &
!		!	+ matmul( transpose(st%y(s,:,:)), pc(s,j)*p(s,:)*st%ovm(s,j,s,:)*(sys%w_qb(s)+st%bigW(s,j,:) &
!		!	+ sys%Ad*dcos(sys%wd*st%t)*st%bigU(s,j,:)) ) &
!		!	- 0.5_8*st%y(s,j,:)*sum( pc(s,j)*p(s,:)*st%ovm(s,j,s,:)*(sys%w_qb(s)+st%bigW(s,j,:) &
!		!	+ sys%Ad*dcos(sys%wd*st%t)*st%bigU(s,j,:)) ) &
!		!	- 0.5_8*p(s,j)*st%y(s,j,:)*sum( pc(s,:)*st%ovm(s,:,s,j)*( sys%w_qb(s) &
!		!	+ st%bigW(s,:,j) + sys%Ad*dcos(sys%wd*st%t)*st%bigU(s,:,j) ) )
!
!
!		!do i=1, sys%nl
!		!	dE_dyc_sjk = dE_dyc_sjk + pc(s,j)*( sys%gk(s,i,:)*sum( p(i,:)*st%ovm(s,j,i,:) )  &
!		!		+ matmul( transpose(st%y(i,:,:)), st%bigL(s,j,i,:)*p(i,:)*st%ovm(s,j,i,:) ) &
!		!		- 0.5_8*st%y(s,j,:)*sum( st%bigL(s,j,i,:)*p(i,:)*st%ovm(s,j,i,:) ) ) &
!		!		- 0.5_8*p(s,j)*st%y(s,j,:)*sum( pc(i,:)*st%ovm(i,:,s,j)*st%bigL(i,:,s,j) )
!		!end do
!
!		do n=1, size( p,2 )
!			dE_dyc_sjk = dE_dyc_sjk + pc(s,j)*p(s,n)*st%ovm(s,j,s,n)*(  &
!				+ sys%wwk(:)*st%y(s,n,:) &
!				+ sys%Ad*dcos( sys%wd*st%t )*sys%Omat(1,:)  &
!				+ (st%y(s,n,:)-0.5_8*st%y(s,j,:))*( sys%w_qb(s)+st%bigW(s,j,n) + sys%Ad*dcos(sys%wd*st%t)*st%bigU(s,j,n) ) ) &
!				- 0.5_8*pc(s,n)*p(s,j)*st%ovm(s,n,s,j)*st%y(s,j,:)*( sys%w_qb(s) &
!				+ st%bigW(s,n,j) &
!				+ sys%Ad*dcos(sys%wd*st%t)*st%bigU(s,n,j) &
!			)
!
!			do i=1, size( p,1 )
!			!do ii=1, sys%nb_el_to_keep
!				!if ( sys%ind( ii, 1 ) == s ) then
!					!i=sys%ind( ii, 2 )
!				
!					dE_dyc_sjk = dE_dyc_sjk + pc(s,j)*p(i,n)*st%ovm(s,j,i,n)*( &
!						+ sys%gk(s,i,:) &
!						+ ( st%y(i,n,:)-0.5_8*st%y(s,j,:) )*st%bigL(s,j,i,n) ) &
!						- 0.5_8*pc(i,n)*p(s,j)*st%ovm(i,n,s,j) &
!						* st%y(s,j,:)*st%bigL(i,n,s,j)
!
!				!end if
!			end do
!
!		end do
!	END FUNCTION

!
!	SUBROUTINE evolve_new_coherent_state(sys,st, tf, dt,  tref )
!		type(param), intent(in)				:: sys
!		type(state), intent(in out)   		:: st
!		real(rl), intent(in)				:: tref, tf, dt
!		integer								:: i, n, max_slow_fact_exp, info
!		real(rl)							:: start_time
!		complex(cx)							:: new_p( size(st%p,1), size(st%p,2) )
!		complex(cx)							:: new_y( size(st%y,1), size(st%y,2), size(st%y,3) )
!		real(rl)							:: tmp_t, new_t
!		character(len=300)		  		 	:: name_ps_l1, name_y0s_l1, name_y0s_l1_prevcs
!		character(len=5)		  		 	:: nb_cs_before_char, nb_cs_now_char
!		complex(cx)							:: y0
!
!		write( nb_cs_before_char, '(I3)' ) st%ncs - 1
!		write( nb_cs_now_char, '(3I3)' ) st%ncs
!		
!		name_ps_l1="data/ADD_l1_ps_"//trim(adjustl(nb_cs_before_char))//"_"//trim(adjustl(nb_cs_now_char))//"_" &
!			//trim(adjustl(parameterchar(sys)))//".d"
!		name_y0s_l1_prevcs="data/ADD_l1_y0s_prevcs_"//trim(adjustl(nb_cs_before_char))//"_"//trim(adjustl(nb_cs_now_char))//"_" &
!			//trim(adjustl(parameterchar(sys)))//".d"
!		name_y0s_l1="data/ADD_l1_y0s_"//trim(adjustl(nb_cs_before_char))//"_"//trim(adjustl(nb_cs_now_char))//"_" &
!			//trim(adjustl(parameterchar(sys)))//".d"
!		open (unit=10,file=name_ps_l1,action="write",status="replace")
!		open (unit=11,file=name_y0s_l1_prevcs,action="write",status="replace")
!		open (unit=12,file=name_y0s_l1,action="write",status="replace")
!
!
!		do n = 1, st%ncs-1
!			y0 = sum( st%y(2,n,:)*sys%Omat(1,:) )
!			write( 11, '(2f25.15)' ) real( y0 ), aimag( y0 )
!		end do
!		write( 11, * )
!
!		max_slow_fact_exp = 1
!		tmp_t = st%t
!		info = 0
!
!		print*, '-- new cs evolve start'
!		WHILE_DO: DO 
!
!			IF ( tmp_t > tf ) exit
!
!			write( 10, '(f25.15)', advance='no' ) tmp_t
!			do n = 1, st%ncs
!				write( 10, '(f25.15)', advance='no' ) abs( st%p(2,n) )**2
!			end do
!			write( 10, * ) 
!
!			y0 = sum( st%y(2,st%ncs,:)*sys%Omat(1,:) )
!			write( 12, '(2f25.15)' ) real( y0 ), aimag( y0 )
!
!			CALL time_evolve( sys%err_lim, sys%w_qb, &
!					sys%wd, sys%Ad, st%t, tref, tmp_t+dt/2, dt, sys%wwk, sys%gk, st%ovm, st%bigW, st%bigL, &
!					st%bigU, sys%Omat, st%p, st%y, .true., sys%lim_slow_fact_exp, new_p, new_y, &
!					new_t, max_slow_fact_exp, info )
!
!			!st%p(:,st%ncs) = new_p(:,st%ncs)
!			st%y(:,st%ncs,:) = new_y(:,st%ncs,:)
!			tmp_t = new_t
!
!			CALL update_sums( sys%wwk, sys%gk, sys%Omat, st%y, st%ovm, st%bigL,st%bigW, st%bigU )
!			print*,'-- in', tmp_t
!
!		END DO WHILE_DO
!		print*, '-- new cs evolve end'
!
!		close(10)
!		close(11)
!		close(12)
!	END SUBROUTINE
!
	!SUBROUTINE time_evolve_RK45( sys, st, tref, dt_add, min_slow_fact_exp, max_slow_fact_exp, slow_fact_exp,&
	!		info, bt_allowed, error_buffer )
	!	type(param), intent(in)       		:: 	sys
	!	type(state), intent(in out)       	:: 	st
	!	real(8), intent(in)       			:: 	tref, dt_add
	!	integer, intent(in out)  	  		:: 	min_slow_fact_exp, max_slow_fact_exp, slow_fact_exp
	!	integer, intent(out)  		  		:: 	info
	!	logical, intent(in)					::  bt_allowed
	!	real(8), intent(in out), optional	::  error_buffer(:)
	!	type(state)       					:: 	midst
	!	real(8)      			  			:: 	dt, error_RK5, t_corrected, last_error_calc, t0
	!	logical								::  dt_corrected
	!	real(8)					  			:: 	a_mat( 5, 5 ), c_mat( 5 ), b_rk4( 6 ), b_rk5( 6 )
	!	complex(cx), dimension(size(st%p,1),size(st%p,2))   :: 	kp1, kp2, kp3, kp4, kp5, kp6
	!	complex(cx), dimension(size(st%y,1), size(st%y,2), size(st%y,3))::  kf1, kf2, kf3, kf4, kf5, kf6
	!	complex(cx), dimension(size(st%p,1),size(st%p,2))             		::  pdot_RK4, pdot_RK5
	!	complex(cx), dimension(size(st%y,1), size(st%y,2), size(st%y,3))     ::  fdot_RK4, fdot_RK5

	!	t0 = st%t
	!	midst = st
	!	t_corrected = -1.0e9_8
	!	dt_corrected = .false.
	!	last_error_calc = -1e3
	!	info = 1

	!	a_mat = transpose( reshape( (/ 1._8/4._8 , 0._8 , 0._8 , 0._8 , 0._8,&
	!		3._8/32._8 , 9._8/32._8 , 0._8 , 0._8 , 0._8,&
	!		1932._8/2197._8 , -7200._8/2197._8 , 7296._8/2197._8 , 0._8, 0._8,&
	!		439._8/216._8 , -8._8 , 3680._8/513._8 , -845._8/4104._8 , 0._8,&
	!		-8._8/27._8 , 2._8 , -3544._8/2565._8 , 1859._8/4104._8 , -11._8/40._8 /),(/5,5/)) )

	!	c_mat =  (/ 1._8/4._8, 3._8/8._8, 12._8/13._8, 1._8, 1._8/2._8 /)

	!	b_rk4 = (/25._8/216._8, 0._8, 1408._8/2565._8, 2197._8/4104._8, -1._8/5._8, 0._8/)
	!	b_rk5 = (/16._8/135._8, 0._8, 6656._8/12825._8, 28561._8/56430._8, -9._8/50._8, 2._8/55._8 /)

	!	WHILE_DO: DO

	!		IF ( ( st%t - t0 > sys%print_delay ) &
	!			.or. (st%t > sys%tmax)  &
	!			.or. (info > 1) ) EXIT WHILE_DO

	!		dt = sys%dt * ( 1._8 /3._8 )**slow_fact_exp

	!		if ( abs(st%t-tref)<sys%dt/100 ) then
	!			dt = dt/1000
	!		elseif ( abs(st%t-tref)<sys%dt/10 ) then
	!			dt = dt/100
	!		elseif ( abs(st%t-tref)<3*sys%dt ) then
	!			dt = dt/10
	!		end if

	!		!===============
	!		!== RK 45
	!		!===============

	!		!== k1 calculation  ==========
	!		if ( .not. dt_corrected ) then
	!			CALL calc_derivatives( sys, st, kp1, kf1 )
	!		end if

	!		!== k2 calculation  ==========
	!		midst%y = st%y + a_mat(1,1)*dt*kf1
	!		midst%p = st%p + a_mat(1,1)*dt*kp1
	!		midst%t = st%t + c_mat(1)*dt
	!		CALL update_sums( sys, midst )
	!		CALL calc_derivatives( sys, midst, kp2, kf2 )

	!		!== k3 calculation  ==========
	!		midst%y = st%y + dt*( a_mat(2,1)*kf1 + a_mat(2,2)*kf2 )
	!		midst%p = st%p + dt*( a_mat(2,1)*kp1 + a_mat(2,2)*kp2 )
	!		midst%t = st%t + c_mat(2)*dt
	!		CALL update_sums( sys, midst )
	!		CALL calc_derivatives( sys, midst, kp3, kf3 )

	!		!== k4 calculation  ==========
	!		midst%y = st%y + dt*( a_mat(3,1)*kf1 + a_mat(3,2)*kf2 + a_mat(3,3)*kf3 )
	!		midst%p = st%p + dt*( a_mat(3,1)*kp1 + a_mat(3,2)*kp2 + a_mat(3,3)*kp3 )
	!		midst%t = st%t + c_mat(3)*dt
	!		CALL update_sums( sys, midst )
	!		CALL calc_derivatives( sys, midst, kp4, kf4 )

	!		!== k5 calculation  ==========
	!		midst%y = st%y + dt*( a_mat(4,1)*kf1 + a_mat(4,2)*kf2 + a_mat(4,3)*kf3 + a_mat(4,4)*kf4 )
	!		midst%p = st%p + dt*( a_mat(4,1)*kp1 + a_mat(4,2)*kp2 + a_mat(4,3)*kp3 + a_mat(4,4)*kp4 )
	!		midst%t = st%t + c_mat(4)*dt
	!		CALL update_sums( sys, midst )
	!		CALL calc_derivatives( sys, midst, kp5, kf5 )

	!		!== k6 calculation  ==========
	!		midst%y = st%y + dt*( a_mat(5,1)*kf1 + a_mat(5,2)*kf2 + a_mat(5,3)*kf3 + a_mat(5,4)*kf4 + a_mat(5,5)*kf5 )
	!		midst%p = st%p + dt*( a_mat(5,1)*kp1 + a_mat(5,2)*kp2 + a_mat(5,3)*kp3 + a_mat(5,4)*kp4 + a_mat(5,5)*kp5 )
	!		midst%t = st%t + c_mat(5)*dt
	!		CALL update_sums( sys, midst )
	!		CALL calc_derivatives( sys, midst, kp6, kf6 )

	!		pdot_RK4 = ( b_rk4(1)*kp1 + b_rk4(2)*kp2 + b_rk4(3)*kp3 + b_rk4(4)*kp4 + b_rk4(5)*kp5 + b_rk4(6)*kp6 )
	!		pdot_RK5 = ( b_rk5(1)*kp1 + b_rk5(2)*kp2 + b_rk5(3)*kp3 + b_rk5(4)*kp4 + b_rk5(5)*kp5 + b_rk5(6)*kp6 )
	!		fdot_RK4 = ( b_rk4(1)*kf1 + b_rk4(2)*kf2 + b_rk4(3)*kf3 + b_rk4(4)*kf4 + b_rk4(5)*kf5 + b_rk4(6)*kf6 )
	!		fdot_RK5 = ( b_rk5(1)*kf1 + b_rk5(2)*kf2 + b_rk5(3)*kf3 + b_rk5(4)*kf4 + b_rk5(5)*kf5 + b_rk5(6)*kf6 )

	!		error_rk5 = dt*( sum( abs(fdot_RK5-fdot_RK4) ) + sum( abs(pdot_RK5-pdot_RK4) ) )

	!		if (  error_rk5 < sys%err_lim  ) then

	!			st%p = st%p + dt*pdot_RK5
	!			st%y = st%y + dt*fdot_RK5
	!			st%t = st%t + dt
	!			CALL update_sums( sys, st )

	!			if (slow_fact_exp < min_slow_fact_exp) then
	!				min_slow_fact_exp = slow_fact_exp
	!			end if
	!			if (slow_fact_exp > max_slow_fact_exp) then
	!				max_slow_fact_exp = slow_fact_exp
	!			end if


	!			if ( (sys%calc_error == 1) &
	!				 .and. present(error_buffer) & 
	!				 .and. (st%t-last_error_calc > sys%error_calc_delay) ) then

	!				!--  calculate and store error in array in rolling array
	!				error_buffer = CSHIFT( error_buffer, -1 )
	!				error_buffer(1) = real( error( sys, st ) )

	!				
	!				if ( ( average_of_five_median(error_buffer) > sys%error_thr ) &
	!				!if ( ( median(error_buffer) > sys%error_thr ) &
	!					.and. ( st%t > tref + dt_add ) &
	!					.and. ( st%ncs < sys%ncs_max ) ) then
	!					print*,'-- ERROR > ERROR_THR: t=', st%t
	!					print*,'-- ERROR > ERROR_THR: median(error)=',median(error_buffer)
	!				end if

	!				if ( ( median(error_buffer) > sys%error_thr ) &
	!					.and. ( st%t > tref + dt_add ) &
	!					.and. ( st%ncs < sys%ncs_max ) ) then

	!					!-- specifying info=2 means program should exit routine and add coherent state
	!					info = 3
	!					
	!				end if

	!			end if

	!			if ( (slow_fact_exp > 0) .and. (st%t-t_corrected > 50*dt) ) then
	!				if (slow_fact_exp >= 2) then
	!					slow_fact_exp = slow_fact_exp - 2
	!				else
	!					slow_fact_exp = slow_fact_exp - 1
	!				endif
	!			end if
	!			dt_corrected = .false.

	!		else 

	!			t_corrected = st%t
	!			slow_fact_exp = slow_fact_exp + 1

	!			if (bt_allowed .eqv. .true.) then
	!				if ( ( ( slow_fact_exp > sys%lim_slow_fact_exp ) .and. (st%ncs<sys%ncs_max) .and. (st%ncs>sys%ncs_ini) ) &
	!					.or. ( ( slow_fact_exp > sys%lim_slow_fact_exp+2 ) .and. (st%ncs==sys%ncs_max) ) &
	!					.or. ( ( slow_fact_exp > sys%lim_slow_fact_exp+4 ) .and. (st%ncs==sys%ncs_ini) ) ) then

	!					CALL print_log( sys, '=========================================', sys%logfile )
	!					CALL print_log( sys, '== ERROR -1: value of slow_fact_exp higher than sys%lim_slow_fact_exp', sys%logfile)
	!					CALL print_log( sys, '== BACKTRACKING ==', sys%logfile)
	!					CALL print_log( sys, '=========================================', sys%logfile)
	!					
	!					!-- specifying info=2 means program should exit routine and backtrack
	!					info = 2

	!				end if
	!			else
	!				if ( slow_fact_exp > 20) then
	!					CALL print_log( sys, '=========================================', sys%logfile )
	!					CALL print_log( sys, '== ERROR -1: value of slow_fact_exp higher than 20', sys%logfile)
	!					CALL print_log( sys, '== ABORT', sys%logfile)
	!				end if
	!			end if

	!			dt_corrected = .true.

	!		end if

	!	END DO WHILE_DO

	!	RETURN
	!END SUBROUTINE time_evolve_RK45
