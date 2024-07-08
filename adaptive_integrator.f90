SUBROUTINE time_evolve_RK45_no_buffer( sys, st )
	type(param), intent(in)       		:: 	sys
	type(state), intent(in out)       	:: 	st
	type(state)       					:: 	midst
	real(8)      			  			:: 	dt, error_RK5, t_corrected,&
												last_error_calc, t0, t_solve, t_other, t_sum
	logical								::  dt_corrected
	real(8)					  			:: 	a_mat( 5, 5 ), c_mat( 5 ), b_rk4( 6 ), b_rk5( 6 )
	complex(8), dimension(size(st%p,1),size(st%p,2))   :: 	kp1, kp2, kp3, kp4, kp5, kp6
	complex(8), dimension(size(st%y,1), size(st%y,2), size(st%y,3))::  kf1, kf2, kf3, kf4, kf5, kf6
	complex(8), dimension(size(st%p,1),size(st%p,2))             		::  pdot_RK4, pdot_RK5
	complex(8), dimension(size(st%y,1), size(st%y,2), size(st%y,3))     ::  fdot_RK4, fdot_RK5

	t0 = st%t
	midst = st
	dt_corrected = .false.

	!-- RK45 parameters
	a_mat = transpose( reshape( (/ 1._8/4._8 , 0._8 , 0._8 , 0._8 , 0._8,&
		3._8/32._8 , 9._8/32._8 , 0._8 , 0._8 , 0._8,&
		1932._8/2197._8 , -7200._8/2197._8 , 7296._8/2197._8 , 0._8, 0._8,&
		439._8/216._8 , -8._8 , 3680._8/513._8 , -845._8/4104._8 , 0._8,&
		-8._8/27._8 , 2._8 , -3544._8/2565._8 , 1859._8/4104._8 , -11._8/40._8 /),(/5,5/)) )

	c_mat =  (/ 1._8/4._8, 3._8/8._8, 12._8/13._8, 1._8, 1._8/2._8 /)

	b_rk4 = (/25._8/216._8, 0._8, 1408._8/2565._8, 2197._8/4104._8, -1._8/5._8, 0._8/)
	b_rk5 = (/16._8/135._8, 0._8, 6656._8/12825._8, 28561._8/56430._8, -9._8/50._8, 2._8/55._8 /)

	WHILE_DO: DO

		!-- exit the loop when t crosses the maximum time tmax
		IF ( st%t > sys%tmax ) EXIT WHILE_DO

		!-- slowfactor is the parameter used to increase or decrease dt
		dt = sys%dt * ( 1._8 /3._8 )**slow_fact_exp

		!===============
		!== RK 45
		!===============

		!== k1 calculation  ==========
		if ( .not. dt_corrected ) then
			!-- (kp1 kf1does not need to be recalculated if state did not change)
			CALL calc_derivatives( sys, st, kp1, kf1, t_solve, t_other )
		end if

		!== k2 calculation  ==========
		midst%y = st%y + a_mat(1,1)*dt*kf1
		midst%p = st%p + a_mat(1,1)*dt*kp1
		midst%t = st%t + c_mat(1)*dt
		CALL calc_derivatives( sys, midst, kp2, kf2, t_solve, t_other  )

		!== k3 calculation  ==========
		midst%y = st%y + dt*( a_mat(2,1)*kf1 + a_mat(2,2)*kf2 )
		midst%p = st%p + dt*( a_mat(2,1)*kp1 + a_mat(2,2)*kp2 )
		midst%t = st%t + c_mat(2)*dt
		CALL calc_derivatives( sys, midst, kp3, kf3, t_solve, t_other  )

		!== k4 calculation  ==========
		midst%y = st%y + dt*( a_mat(3,1)*kf1 + a_mat(3,2)*kf2 + a_mat(3,3)*kf3 )
		midst%p = st%p + dt*( a_mat(3,1)*kp1 + a_mat(3,2)*kp2 + a_mat(3,3)*kp3 )
		midst%t = st%t + c_mat(3)*dt
		CALL calc_derivatives( sys, midst, kp4, kf4, t_solve, t_other  )

		!== k5 calculation  ==========
		midst%y = st%y + dt*( a_mat(4,1)*kf1 + a_mat(4,2)*kf2 + a_mat(4,3)*kf3 + a_mat(4,4)*kf4 )
		midst%p = st%p + dt*( a_mat(4,1)*kp1 + a_mat(4,2)*kp2 + a_mat(4,3)*kp3 + a_mat(4,4)*kp4 )
		midst%t = st%t + c_mat(4)*dt
		CALL calc_derivatives( sys, midst, kp5, kf5, t_solve, t_other  )

		!== k6 calculation  ==========
		midst%y = st%y + dt*( a_mat(5,1)*kf1 + a_mat(5,2)*kf2 + a_mat(5,3)*kf3 + a_mat(5,4)*kf4 + a_mat(5,5)*kf5 )
		midst%p = st%p + dt*( a_mat(5,1)*kp1 + a_mat(5,2)*kp2 + a_mat(5,3)*kp3 + a_mat(5,4)*kp4 + a_mat(5,5)*kp5 )
		midst%t = st%t + c_mat(5)*dt
		CALL calc_derivatives( sys, midst, kp6, kf6, t_solve, t_other  )

		!== final derivatives
		pdot_RK4 = ( b_rk4(1)*kp1 + b_rk4(2)*kp2 + b_rk4(3)*kp3 + b_rk4(4)*kp4 + b_rk4(5)*kp5 + b_rk4(6)*kp6 )
		pdot_RK5 = ( b_rk5(1)*kp1 + b_rk5(2)*kp2 + b_rk5(3)*kp3 + b_rk5(4)*kp4 + b_rk5(5)*kp5 + b_rk5(6)*kp6 )
		fdot_RK4 = ( b_rk4(1)*kf1 + b_rk4(2)*kf2 + b_rk4(3)*kf3 + b_rk4(4)*kf4 + b_rk4(5)*kf5 + b_rk4(6)*kf6 )
		fdot_RK5 = ( b_rk5(1)*kf1 + b_rk5(2)*kf2 + b_rk5(3)*kf3 + b_rk5(4)*kf4 + b_rk5(5)*kf5 + b_rk5(6)*kf6 )

		!== error calculated based on RK4 and RK5 results
		error_rk5 = dt*( sum( abs(fdot_RK5-fdot_RK4) ) + sum( abs(pdot_RK5-pdot_RK4) ) )

		if (  error_rk5 < sys%err_lim  ) then
			!== if error below threshold, update the state variables, and continue with next step
			st%p = st%p + dt*pdot_RK5
			st%y = st%y + dt*fdot_RK5
			st%t = st%t + dt
			CALL update_sums( sys, st, t_sum )

			!== incerase dt in case it had been decreased (if 50 steps have elapsed)
			if ( (slow_fact_exp > 0) .and. (st%t-t_corrected > 50*dt) ) then
				if (slow_fact_exp >= 2) then
					slow_fact_exp = slow_fact_exp - 2
				else
					slow_fact_exp = slow_fact_exp - 1
				endif
			end if
			dt_corrected = .false.  !-- record that the dt was NOT corrected at this time step

		else 

			!== if error is too large, disacrd result and decrease dt
			t_corrected = st%t
			slow_fact_exp = slow_fact_exp + 1
			dt_corrected = .true.

			!== if dt was slowed beyond a certain value, assume there is an error and abort the code
			if ( slow_fact_exp > 20) then
				CALL print_log( sys, '=========================================', sys%logfile )
				CALL print_log( sys, '== ERROR -1: value of slow_fact_exp higher than 20', sys%logfile)
				CALL print_log( sys, '== ABORT', sys%logfile)
			end if


		end if

	END DO WHILE_DO

	RETURN
END SUBROUTINE time_evolve_RK45_no_buffer
