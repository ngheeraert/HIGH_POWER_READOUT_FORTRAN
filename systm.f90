MODULE SYSTM 

	USE lapackblas, only: diagonalise_all
	USE typedefs, only : cx => c_type, rl => r_type

	IMPLICIT NONE 

	TYPE state
		real(rl)					::  t 								!-- Time
		integer    					::  ncs
		complex(cx), allocatable  	::  y(:,:,:), y0(:,:)   				!-- value of the fiel displacements and momenta
		complex(cx), allocatable  	::  p(:,:)   					!-- probability amplitudes of the polarons
		complex(cx), allocatable  	::  ovm(:,:,:,:), bigL(:,:,:,:)  	!-- matrices of overlaps
		complex(cx), allocatable  	::  bigU(:,:,:), bigW(:,:,:)
	END TYPE state

	TYPE param
		character(len=5)			::   method
		integer						::   nl
		integer						::   ncs_ini 
		integer						::   ncs_max
		integer						::   median_pts
		integer  					::   nmodes
		integer	    				::   nbathmodes
		integer						::   qb_ini
		integer						::   lim_slow_fact_exp
		integer						::   logfile
		integer						::   job_nb
		integer						::   PID
		integer						::   n_dtadd_change
		integer						::   expnt
		real(rl)					::   print_delay
		real(rl)					::   gcut
		real(rl)					::   error_calc_delay
		real(rl)					::   p2_threshold
		real(rl)					::   dt_print
		real(rl)					::   cav_ini
		real(rl)					::   wc
		real(rl)					::   alpha
		real(rl)					::   w_ge
		real(rl)					::   bw
		character(len=1)			::   bath_type
		real(rl)					::   bc_cav
		real(rl)					::   bc_qb
		real(rl)					::   bc_lf
		real(rl)					::   bc_hf
		real(rl)					::   mode_ratio_qb
		real(rl)					::   mode_ratio_lf
		real(rl)					::   mode_ratio_hf
		real(rl)					::   ad
		real(rl)					::   tmax
		real(rl)					::   dt
		real(rl)					::   err_lim
		real(rl)					::   error_thr
		real(rl)					::   dt_add
		real(rl)					::   dt_add_ini
		real(rl)					::   p0
		real(rl)					::   g_qc
		real(rl), allocatable		::   dtadd_arr(:)
		real(rl), allocatable		::   wk(:)
		real(rl), allocatable		::   wk1(:)
		real(rl), allocatable		::   dk(:)
		real(rl), allocatable		::   wwk(:)
		real(rl), allocatable		::   w_qb(:)
		real(rl), allocatable		::   g_qc_arr(:,:)
		real(rl), allocatable		::   gk(:,:,:)
		real(rl), allocatable		::   gij(:,:)
		real(rl), allocatable		::   ind(:,:)
		integer						::   nb_el_to_keep
		real(rl), allocatable		::   Omat(:,:)
		real(rl)					::   dw
		real(rl)					::   wd
		real(rl)					::   anh
		real(rl)					::   adding_ratio
		real(rl)					::   offset_abs
		real(rl)					::   offset_ang
		complex(cx), allocatable	::   g2mat(:,:,:,:)
		complex(cx), allocatable	::   sum_diag_g2mat(:,:)
		logical						::   time_bar
		integer						::   calc_error
		character(len=5)			::   dvice
		type(state)					::   ground_state
		type(state)					::   excited_state
	END TYPE PARAM

	TYPE TRAJ
		integer						::  i_time
		real(rl), allocatable  	    ::  time_ar(:), populations_ar(:,:,:), norm_ar(:), error_ar(:,:),&
											ps2_ar(:,:,:), y02_ar(:,:,:),n_ar(:,:), renyi_entropy(:),&
											ovm_eval_0_ar(:,:), ovm_eval_1_ar(:,:), apadag_ar(:)
	END TYPE

	real(8),   parameter :: pi=(4._8)*atan(1._8)
	real(8),   parameter :: twopi=2.0_8*(4.0_8)*atan(1.0_8)
	real(8),   parameter :: one=1.0_8
	real(8),   parameter :: zero=0.0_8
	complex(8),parameter :: ic=(0.0_8,1.0_8)

CONTAINS

	!== initialize the parameters and get those in the command line
	SUBROUTINE get_parameters(sys)
		type(param)           			::  sys
		character(len=300)    			::  buffer, charge_op_filename, energy_filename
		character(len=2)    			::  jobnb_char
		integer 				  		::  i, j, k, s, ss,n, nargs, nmodes_qb, nmodes_cav, nmodes_lf, nmodes_hf, mu
		real(rl), allocatable   		::  coupling_op(:,:), h_pp(:,:)
		real(rl)						::  min_w_cav, max_w_cav, min_w_qb, min_w_lf, min_w_hf, dw, mode_ratio_cav
		real(rl)						::  offset_abs, offset_ang
		real(rl), allocatable 			::  gk_abs_sum(:,:)
		integer, allocatable 			::  gk_abs_keep(:,:)
		integer							::  nb_terms_to_keep_on_diag 
		real(8)							::  gk_cutoff, offset, one_photon_f, p1, p2, dummy

		!== these are the parameters provided in the command line 
		nargs = iargc()
		call get_command_argument(1, buffer)

		!-- initialise the log file and overwirte any existing one
		sys%tmax=1
		sys%nl=2
		sys%wc=1.25
		sys%p2_threshold=1e-5
		sys%g_qc=0.016
		sys%alpha=0.01
		sys%mode_ratio_qb=0.0
		sys%mode_ratio_lf=0.0
		sys%mode_ratio_hf=0.0
		sys%nbathmodes = 100
		sys%lim_slow_fact_exp = 10
		sys%bw=0.4
		sys%bath_type='F'
		sys%bc_cav = 1
		sys%bc_qb = 1
		sys%bc_lf = 1
		sys%bc_hf = 1
		sys%wd=1.0
		sys%w_ge=1.0
		sys%Ad=0.0
		sys%tmax=10
		sys%dt=0.1_8
		sys%err_lim=1e-7
		sys%median_pts=25
		sys%adding_ratio=0._8
		sys%error_thr=1e-7
		sys%p0=1e-6
		sys%ncs_ini=1
		sys%ncs_max=1
		sys%dt_add=0.1
		sys%dt_add_ini=sys%dt_add
		sys%n_dtadd_change=4
		sys%qb_ini=1
		sys%cav_ini=0 
		sys%offset_abs=1._rl
		sys%offset_ang=0._rl
		sys%logfile=1
		sys%dvice = 'TRSM1'
		sys%print_delay = sys%dt*50
		sys%error_calc_delay = sys%dt*10
		sys%method='RK45'
		sys%gcut=-1

		sys%time_bar=.true.
		sys%calc_error=1

		if(nargs==0)then
			stop '# Use input parameters parameters'
		else
			do i=1,nargs
				call get_command_argument(i, buffer)
				if(buffer=='-nl')then 
					call get_command_argument(i+1, buffer) 
					read(buffer,*) sys%nl
				else if(buffer=='-p0')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%p0
				else if(buffer=='-tmax')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%tmax
				else if(buffer=='-nbm')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%nbathmodes
				elseif(buffer=='-ncs_ini')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%ncs_ini
				elseif(buffer=='-ncs_max')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%ncs_max
				elseif(buffer=='-cav_ini')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%cav_ini
				elseif(buffer=='-qb_ini')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%qb_ini
				else if(buffer=='-al')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%alpha
					sys%alpha = sys%alpha * twopi
				else if(buffer=='-offset_abs')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%offset_abs
				else if(buffer=='-offset_ang')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%offset_ang
				else if(buffer=='-wd')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%wd
					sys%wd = sys%wd * 2*pi
				else if(buffer=='-w_ge')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%w_ge
					sys%w_ge = sys%w_ge * 2*pi
				else if(buffer=='-wc')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%wc
					sys%wc = sys%wc * 2*pi
				else if(buffer=='-g_qc')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%g_qc
					sys%g_qc = sys%g_qc * 2*pi
				else if(buffer=='-dt')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%dt
				else if(buffer=='-dt_add')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%dt_add
				else if(buffer=='-n_dtadd_change')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%n_dtadd_change
				else if(buffer=='-dt_add_ini')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%dt_add_ini
				else if(buffer=='-ad')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%ad
					sys%Ad = sys%Ad * 2*pi
				else if(buffer=='-bc_cav')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%bc_cav
					sys%bc_cav = sys%bc_cav * 2*pi
				else if(buffer=='-bc_qb')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%bc_qb
					sys%bc_qb = sys%bc_qb * 2*pi
				else if(buffer=='-bc_lf')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%bc_lf
					sys%bc_lf = sys%bc_lf * 2*pi
				else if(buffer=='-bc_hf')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%bc_hf
					sys%bc_hf = sys%bc_hf * 2*pi
				else if(buffer=='-mr_qb')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%mode_ratio_qb
				else if(buffer=='-mr_lf')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%mode_ratio_lf
				else if(buffer=='-mr_hf')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%mode_ratio_hf
				else if(buffer=='-bw')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%bw
					sys%bw = sys%bw * 2*pi
				else if(buffer=='-bath_type')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%bath_type
				else if(buffer=='-lsfe')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%lim_slow_fact_exp
				else if(buffer=='-errlim')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%err_lim
				else if(buffer=='-error_thr')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%error_thr
				else if(buffer=='-median_pts')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%median_pts
				else if(buffer=='-adding_ratio')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%adding_ratio
				else if(buffer=='-dvice')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%dvice
				else if(buffer=='-method')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%method
				else if(buffer=='-logfile')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%logfile
				else if(buffer=='-p2_threshold')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%p2_threshold
				else if(buffer=='-print_delay')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%print_delay
				else if(buffer=='-error_calc_delay')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%error_calc_delay
				else if(buffer=='-calc_error')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%calc_error
				else if(buffer=='-job_nb')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%job_nb
				else if(buffer=='-gcut')then 
					call get_command_argument(i+1, buffer)
					read(buffer,*) sys%gcut
				end if

			end do
		end if
		
		print*,'8888888888888', sys%alpha
		!if (sys%job_nb>0) then
		!	write(jobnb_char, '(I2)') sys%job_nb
		!	open (unit=1,file='PID_'//trim(adjustl(jobnb_char)),action="read",status="old")
		!	read(1,'(I5)',advance='no') sys%PID
		!	close(1)
		!else
		!	sys%PID = 0
		!end if
		allocate( sys%dtadd_arr( sys%ncs_max ) )
		sys%dtadd_arr = 1.0e5
		do n = sys%ncs_ini, sys%ncs_max
			if (n >= sys%n_dtadd_change) then
				sys%dtadd_arr(n) = sys%dt_add
			else
				sys%dtadd_arr(n) = sys%dt_add_ini
			end if
		end do

		allocate( sys%w_qb(sys%nl) )
		allocate( coupling_op(sys%nl,sys%nl) )
		!allocate( sys%g_qc_arr(sys%nl,sys%nl) )
		allocate( sys%gij(sys%nl,sys%nl) )


		IF (sys%dvice == 'TRSM3') then
			open (unit=100,file='qubit_params/T_FOR_E_lvl_Ec0.280_Ej14.000.txt', action="read",status="old")
			open (unit=101,file='qubit_params/T_FOR_Charge_Mat_Ec0.280_Ej14.000.txt', action="read",status="old")
			sys%expnt = 1
		ELSEIF (sys%dvice == 'TRSM2') then
			!-- here the coupling 01 is normalised to 1
			open (unit=100,file='qubit_params/FOR_E_lvl_Ec0.190_Ej14.368.txt', action="read",status="old")
			open (unit=101,file='qubit_params/FOR_Charge_Mat_Ec0.190_Ej14.368.txt', action="read",status="old")
			sys%expnt = 1
		ELSEIF (sys%dvice == 'QUTR2') then
			!-- here the coupling 01 is normalised to 1
			open (unit=100,file='qubit_params/FOR_E_lvl_Ec0.190_Ej14.368.txt', action="read",status="old")
			open (unit=101,file='qubit_params/FOR_Cos_Phi_q_Mat_Ec0.190_Ej14.368.txt', action="read",status="old")
			sys%expnt = 2
		ELSE
			print*,'--ERROR IN DEVICE SPECIFICTION'
		END IF

		do i=1,sys%nl
			read(100,'(f25.15)',advance='yes') sys%w_qb(i)
		end do
		do i=1,sys%nl
			do j=1,sys%nl
				read(101,'(f20.15)',advance='no') coupling_op(i,j)
			end do
			read(101,*)
		end do
		sys%w_qb = sys%w_qb * 2*pi
		sys%w_ge = sys%w_qb(2) - sys%w_qb(1)
		!sys%g_qc_arr = - sys%g_qc*coupling_op
		sys%gij = sys%g_qc*coupling_op( 1:sys%nl, 1:sys%nl )
		!sys%g_qc_arr = sys%g_qc*coupling_op
		!sys%gij = -sys%g_qc*coupling_op

		!-- for the quantromon no minus as it comes from the  cosine


		if (sys%nl>2) then
			sys%anh = sys%w_qb(3) - 2*sys%w_qb(2)
		else
			sys%anh = 1000
		end if

		sys%nmodes = sys%nbathmodes + 1

		if (sys%logfile >= 1) then
			open (unit=100,file='data/LOG_'//trim(adjustl(parameterchar(sys)))//'.d' ,action="write",status="replace")
			write(100 ,*) '========== LOG FILE ============'
			write(100,*) ''
		end if
		CALL print_log( sys, '-- simulation parameter char:', sys%logfile  )
		CALL print_log( sys, parameterchar(sys), sys%logfile  )
		CALL print_log( sys, '', sys%logfile )

		!=======================================  -1.4217379764058791E-004
		!== setting up de mode values array

		allocate( sys%wk(sys%nmodes) )

		!-- first mode is the cavity mode
		sys%wk( 1 ) = sys%wc

		!-- setting the mode range that is around the qubit and cavity frequencies
		min_w_qb = sys%bc_qb - 0.5_rl*sys%mode_ratio_qb*sys%bw
		min_w_lf = sys%bc_lf - 0.5_rl*sys%mode_ratio_lf*sys%bw
		min_w_hf = sys%bc_hf - 0.5_rl*sys%mode_ratio_hf*sys%bw
		mode_ratio_cav = 1.0-sys%mode_ratio_qb-sys%mode_ratio_lf-sys%mode_ratio_hf
		min_w_cav = sys%bc_cav - 0.5_rl*mode_ratio_cav*sys%bw
		max_w_cav = sys%bc_cav + 0.5_rl*mode_ratio_cav*sys%bw
		nmodes_qb = int(sys%mode_ratio_qb*sys%nbathmodes)
		nmodes_lf = int(sys%mode_ratio_lf*sys%nbathmodes)
		nmodes_hf = int(sys%mode_ratio_hf*sys%nbathmodes)
		nmodes_cav = sys%nbathmodes - nmodes_qb - nmodes_lf - nmodes_hf
		dw = ( max_w_cav - min_w_cav )/nmodes_cav
		!-- lf modes
		do i=1, nmodes_lf
			sys%wk( 1+i ) = min_w_lf + dw*(i-1)
		end do
		!-- qubit
		do i=1, nmodes_qb
			sys%wk( 1+nmodes_lf+i ) = min_w_qb + dw*(i-1)
		end do
		!-- cavity
		do i=1, nmodes_cav
			sys%wk( 1+nmodes_qb+nmodes_lf+i ) = min_w_cav + dw*(i-1)
		end do
		!-- hf modes
		do i=1, nmodes_hf
			sys%wk( 1+nmodes_qb+nmodes_lf+nmodes_cav+i ) = min_w_hf + dw*(i-1)
		end do
		sys%dw = sys%wk(3) - sys%wk(2)

		!=======================================
		!== diagonalisaiton the linear part of the Hamiltonian

		allocate( h_pp( sys%nmodes,sys%nmodes ) )
		allocate( sys%wwk( sys%nmodes ) )
		!allocate( sys%gk( sys%nl,sys%nl,sys%nmodes ) )
		allocate( sys%Omat( sys%nmodes, sys%nmodes ) )

		do i=1, sys%nmodes
			h_pp(i,i) = sys%wk(i)
			if (i>1) then
				if ( sys%bath_type=='O' ) then
					h_pp(1,i) = dsqrt( 2._8*sys%alpha*sys%dw*sys%wk(i) )
					h_pp(i,1) = dsqrt( 2._8*sys%alpha*sys%dw*sys%wk(i) )
				else if ( sys%bath_type=='F' ) then
					h_pp(1,i) = dsqrt( 2._8*sys%alpha*sys%dw )
					h_pp(i,1) = dsqrt( 2._8*sys%alpha*sys%dw )
				end if
			end if
		end do

		call diagonalise_all( h_pp, sys%Omat, sys%wwk )

		do i=1, sys%nl
			print*, sys%gij( i , : )
		end do

		allocate( gk_abs_keep( sys%nl, sys%nl ) )
		gk_abs_keep = 0
		do i=1, sys%nl
			do j=1, sys%nl
				if ( abs( sys%gij(i,j) ) > sys%gcut*abs( sys%gij(1,2) )  ) then
					gk_abs_keep(i,j) = 1
				else
					sys%gij(i,j) = 0._8
				end if
			end do
		end do

		allocate( sys%ind( sum(gk_abs_keep), 2 ) )
		sys%nb_el_to_keep = 0
		do i=1, sys%nl
			do j=1, sys%nl
				if ( gk_abs_keep(i,j) == 1 ) then
					sys%nb_el_to_keep = sys%nb_el_to_keep + 1
					sys%ind( sys%nb_el_to_keep, 1 ) = i
					sys%ind( sys%nb_el_to_keep, 2 ) = j
				end if
			end do
		end do
		print*,'-- for error calculation, keeping: ', sys%nb_el_to_keep, ' terms'
		print*,'-- ', sys%nl**2 - sys%nb_el_to_keep, ' terms removed'

		do i=1, size( sys%ind,1 )
			print*,sys%gij( sys%ind(i,1) , sys%ind(i,2) )
		end do

		do i=1, sys%nl
			print*, sys%gij( i , : )
		end do


		!do i=1, sys%nl
		!	do j=1, sys%nl
		!		sys%gk(i,j,:) = sys%Omat(1,:) * sys%gij(i,j)
		!	enddo
		!end do

		!allocate( gk_abs_sum( sys%nl, sys%nl ) )
		!allocate( gk_abs_keep( sys%nl, sys%nl ) )
		!gk_abs_sum = sum( abs(sys%gk(:,:,:)), dim=3 )
		!gk_abs_keep = 0
		!gk_cutoff = -0.001
		!do i=1, sys%nl
		!	do j=1, sys%nl
		!		if ( gk_abs_sum(i,j) > gk_cutoff*gk_abs_sum(1,2)  ) then
		!			gk_abs_keep(i,j) = 1
		!		!else
		!	!		sys%gk(i,j,:) = 0._8
		!		end if
		!	end do
		!end do
		!   
		!allocate( sys%ind( sum(gk_abs_keep), 2 ) )
		!sys%nb_el_to_keep = 0
		!do i=1, sys%nl
		!	do j=1, sys%nl
		!		if ( gk_abs_keep(i,j) == 1 ) then
		!			sys%nb_el_to_keep = sys%nb_el_to_keep + 1
		!			sys%ind( sys%nb_el_to_keep, 1 ) = i
		!			sys%ind( sys%nb_el_to_keep, 2 ) = j
		!		end if
		!	end do
		!end do
		!!print*,'-- for error calculation, keeping: ', sys%nb_el_to_keep, ' terms'
		!!print*,'-- ', sys%nl**2 - sys%nb_el_to_keep, ' terms removed'

		!!do i=1, size( sys%ind_D )
		!!	print*,'here'
		!!	print*,gk_abs_sum( sys%ind_D(i) , sys%ind_D(i) ) 
		!!end do
		

		!if (sys%calc_error == 1) then

		!	allocate( sys%g2mat( sys%nl,sys%nl,sys%nmodes,sys%nmodes ) )
		!	sys%g2mat = 0._8
		!	do s=1, sys%nl
		!		do i=1, sys%nl
		!			do j=1, sys%nl
		!				do ss=1, sys%nmodes
		!					sys%g2mat(i,j,ss,:) = sys%g2mat(i,j,ss,:) + sys%gk(s,i,ss)*sys%gk(s,j,:)
		!				end do
		!			end do
		!		enddo
		!	enddo
		!	allocate( sys%sum_diag_g2mat( sys%nl,sys%nl ) )
		!	sys%sum_diag_g2mat = 0._8
		!	do i=1, sys%nl
		!		do j=1, sys%nl
		!			do k=1, sys%nmodes
		!				sys%sum_diag_g2mat(i,j) = sys%sum_diag_g2mat(i,j) + sys%g2mat(i,j,k,k)
		!			end do
		!		end do
		!	end do

		!endif

		!================================================
		!-- defining dressed ground and excited states

		CALL allocate_state(sys,sys%ground_state,sys%ncs_ini)
		CALL allocate_state(sys,sys%excited_state,sys%ncs_ini)

		offset = sys%offset_abs !* exp( Ic*sys%offset_ang )
		one_photon_f = 0.01

		do n=1, sys%ncs_ini
			sys%ground_state%y(:,n,:) = 0.01*n
			sys%excited_state%y(:,n,:) = 0.01*n
		end do
		sys%ground_state%p = sys%p0
		sys%excited_state%p = sys%p0

		sys%ground_state%p(1,1) = 1.0
		sys%ground_state%p(1,2) = 1.0
		sys%ground_state%p(2,1) = 1/(2*one_photon_f)
		sys%ground_state%p(2,2) = - 1/(2*one_photon_f)
		do mu=1, sys%nmodes
			sys%ground_state%y(2,1,mu) = one_photon_f*( sys%Omat(1,mu) )
			sys%ground_state%y(2,2,mu) = - one_photon_f*( sys%Omat(1,mu) )
		end do
		CALL update_sums( sys, sys%ground_state, dummy )
		CALL normalise( sys%ground_state )

		sys%excited_state%p(2,:) = sys%ground_state%p(1,:)
		sys%excited_state%p(1,:) = sys%ground_state%p(2,:)
		sys%excited_state%y(2,:,:) = sys%ground_state%y(1,:,:)
		sys%excited_state%y(1,:,:) = sys%ground_state%y(2,:,:)
		CALL update_sums( sys, sys%excited_state, dummy )

		p1 = lvl_occupation( sys%ground_state, 1 )
		p2 = lvl_occupation( sys%ground_state, 2 )
		sys%ground_state%p(1,:) = -0.999774217*sys%ground_state%p(1,:)/sqrt(p1)
		sys%ground_state%p(2,:) = -0.0212418917*sys%ground_state%p(2,:)/sqrt(p2)
		CALL normalise( sys%ground_state )

		p1 = lvl_occupation( sys%excited_state, 1 )
		p2 = lvl_occupation( sys%excited_state, 2 )
		sys%excited_state%p(1,:) = -0.12160415*sys%excited_state%p(1,:)/sqrt(p1)
		sys%excited_state%p(2,:) = -0.9925738*sys%excited_state%p(2,:)/sqrt(p2)
		CALL normalise( sys%excited_state )

		!print*, '-- PROBABILITIES'
		!print*, 'g1',lvl_occupation( sys%ground_state, 1 )
		!print*, 'g2',lvl_occupation( sys%ground_state, 2 )
		!print*, 'e1',lvl_occupation( sys%excited_state, 1 )
		!print*, 'e2',lvl_occupation( sys%excited_state, 2 )
		!print*, '-- MEAN_PHOTON_IN_CAVITY'
		!print*, 'g',photon_nb_in_sig_mode( sys, sys%ground_state, 1 )
		!print*, 'e',photon_nb_in_sig_mode( sys, sys%excited_state, 1 )

  		CALL print_log( sys, "-- parameters initialised", sys%logfile )
	END SUBROUTINE  

	SUBROUTINE allocate_state(sys,st,ncs_value) 
		type(param), intent(in)      	::  sys
		type(state), intent(out)  		::  st
		type(state)							::  tmpst !-- if st has to be reallocated
		integer, intent(in),optional		::  ncs_value
		integer									::  ncs, step_numb

		if (present(ncs_value)) then
			tmpst%ncs = ncs_value
		else
			tmpst%ncs = sys%ncs_ini
		end if
		ncs = tmpst%ncs

		allocate(tmpst%y(sys%nl,ncs,sys%nmodes))
		allocate(tmpst%y0(sys%nl,ncs))
		allocate(tmpst%p(sys%nl,ncs))
		allocate(tmpst%ovm(sys%nl,ncs,sys%nl,ncs)) 
		allocate(tmpst%bigU(sys%nl,ncs,ncs)) 
		allocate(tmpst%bigW(sys%nl,ncs,ncs)) 
		allocate(tmpst%bigL(sys%nl,ncs, sys%nl,ncs)) 

		tmpst%t = 0._rl
		tmpst%y(:,:,:) = 0._cx
		tmpst%y0(:,:) = 0._cx
		tmpst%p(:,:) = 0._cx
		tmpst%ovm(:,:,:,:) = 0._cx
		tmpst%bigL(:,:,:,:) = 0._cx
		tmpst%bigW(:,:,:) = 0._cx
		tmpst%bigU(:,:,:) = 0._cx
		st = tmpst
	END SUBROUTINE

	SUBROUTINE allocate_state_old(sys,st,ncs_value) 
		type(param), intent(in)      	::  sys
		type(state), intent(out)  		::  st
		type(state)							::  tmpst !-- if st has to be reallocated
		integer, intent(in),optional		::  ncs_value
		integer									::  ncs, step_numb

		if (present(ncs_value)) then
			tmpst%ncs = ncs_value
		else
			tmpst%ncs = sys%ncs_ini
		end if
		ncs = tmpst%ncs

		allocate(tmpst%y(sys%nl,ncs,sys%nmodes))
		allocate(tmpst%p(sys%nl,ncs))
		allocate(tmpst%ovm(sys%nl,ncs,sys%nl,ncs)) 
		allocate(tmpst%bigU(sys%nl,ncs,ncs)) 
		allocate(tmpst%bigW(sys%nl,ncs,ncs)) 
		allocate(tmpst%bigL(sys%nl,ncs, sys%nl,ncs)) 

		tmpst%t = 0._rl
		tmpst%y(:,:,:) = 0._cx
		tmpst%p(:,:) = 0._cx
		tmpst%ovm(:,:,:,:) = 0._cx
		tmpst%bigL(:,:,:,:) = 0._cx
		tmpst%bigW(:,:,:) = 0._cx
		tmpst%bigU(:,:,:) = 0._cx
		st = tmpst
	END SUBROUTINE

	SUBROUTINE allocate_trajectory(sys,tr, step_numb)
		type(param), intent(in)      	::  sys
		type(traj), intent(out)      	::  tr
		integer, intent(in)				::  step_numb

		allocate( tr%time_ar(step_numb) )
		allocate( tr%error_ar(step_numb,2) )
		allocate( tr%populations_ar(step_numb, sys%nl, 2) )
		allocate( tr%norm_ar(step_numb) )
		allocate( tr%ps2_ar(step_numb,sys%nl,sys%ncs_max) )
		allocate( tr%y02_ar(step_numb,sys%nl,sys%ncs_max) )
		allocate( tr%n_ar(step_numb,sys%nmodes) )
		allocate( tr%renyi_entropy(step_numb) )
		allocate( tr%ovm_eval_0_ar(step_numb, sys%ncs_max) )
		allocate( tr%ovm_eval_1_ar(step_numb, sys%ncs_max) )
		allocate( tr%apadag_ar(step_numb) )

		tr%i_time=0
		tr%time_ar = 0._rl
		tr%populations_ar = 0._rl
		tr%norm_ar = 0._rl
		tr%error_ar = 0._rl
		tr%ps2_ar = 1.0e-8_rl
		tr%y02_ar = 0._rl
		tr%n_ar = 0._rl
		tr%renyi_entropy = 0._rl
		tr%ovm_eval_0_ar = 1e6
		tr%ovm_eval_1_ar = 1e6
		tr%apadag_ar = 0._rl
	END SUBROUTINE

	SUBROUTINE initialise_dressed_state(sys,st,time_in)
		type(param), intent(in)      		::  sys
		type(state), intent(in out) 	 	::  st
		real(rl), intent(in)				::  time_in

		if ( sys%qb_ini == -1 )  then
			st = sys%ground_state
			print*,'-- qubit initialised in dressed GROUND state.'
		else if ( sys%qb_ini == -2 )  then
			st = sys%excited_state
			print*,'-- qubit initialised in dressed EXCITED state.'
		end if
	END SUBROUTINE

	SUBROUTINE initialise_bare_state(sys,st,time_in)
		type(param), intent(in)      		::  sys
		type(state), intent(in out) 	 	::  st
		real(rl), intent(in)				::  time_in
		integer								::  l,mu,n
		real(8)								::  one_photon_f, dummy
		complex(cx)							::  offset

		st%t = time_in
		st%y = 0._rl
		st%p = sys%p0
		one_photon_f = 0.1
		offset = sys%offset_abs * exp( Ic*sys%offset_ang )
		
		if (sys%cav_ini > -1e-7) then

			do l=1,sys%nl

				st%y(l,1,:) = sys%cav_ini * sys%Omat(1,:) 

				if (st%ncs > 1) then
					do n=2, st%ncs
						st%y(l,n,:) = ( sys%cav_ini + offset ) * sys%Omat(1,:)
					end do
				end if

			end do

			st%p( sys%qb_ini, 1 ) = 1.0

		else 

			do n=1, st%ncs
				st%y(:,n,:) = offset*n
			end do

			st%y(sys%qb_ini,:,:) = 0.0
			st%p(sys%qb_ini,1) = 1/(2*one_photon_f)
			st%p(sys%qb_ini,2) = - 1/(2*one_photon_f)
			do mu=1, sys%nmodes
				st%y(sys%qb_ini,1,mu) = one_photon_f*( sys%Omat(1,mu) )
				st%y(sys%qb_ini,2,mu) = - one_photon_f*( sys%Omat(1,mu) )
			end do

		end if

		!-- updating the sums over k
		call update_sums( sys, st, dummy )
		CALL normalise(st)
	END SUBROUTINE

	SUBROUTINE update_sums( sys, st, sum_time )
		type( param ), intent(in)		   :: sys
		type( state ), intent(in out)	   :: st
		real(8), intent(out)		   :: sum_time
		real(8)		   			   :: t1, t2, tmp
		integer ::  m,n,i,j,ii,jj

		call cpu_time(t1)
		do i=1, size(st%y,1)
			st%y0(i,:) = matmul( st%y(i,:,:), sys%Omat(1,:) )
		end do

		!st%ovm = 1._8
		!do i=1, size(st%y,1)
		!	do n=1, size(st%y,2)
		!		tmp = exp( sum( -0.5_8*conjg(st%y(i,n,:))*st%y(i,n,:)) )
		!		st%ovm(i,n,:,:) = st%ovm(i,n,:,:) * tmp
		!		st%ovm(:,:,i,n) = st%ovm(:,:,i,n) * tmp
		!	end do
		!end do
		st%bigU = 0
		do n=1, size(st%y,2)
			st%bigU(:,:,n) = st%bigU(:,:,n) + conjg( st%y0(:,:) )
			st%bigU(:,n,:) = st%bigU(:,n,:) + st%y0(:,:)
		end do

		do m=1, size(st%y,2)

			!-- m == n
			st%bigW(:,m,m) = matmul( conjg(st%y(:,m,:)) * st%y(:,m,:), sys%wwk(:) )
			do i=1, size(st%y,1)
				st%ovm(i,m,i,m) = 1
				st%bigL(i,m,i,m) = sys%gij(i,i) * ( conjg(st%y0(i,m)) + st%y0(i,m) )**sys%expnt
			end do
			do i=1, size(st%y,1)
				do j=i+1, size(st%y,1)
					st%ovm(i,m,j,m) = overlap( st%y(i,m,:), st%y(j,m,:) )
					st%ovm(j,m,i,m) = conjg( st%ovm(i,m,j,m) ) 
					!st%ovm(i,m,j,m) = st%ovm(i,m,j,m)*exp( sum( conjg(st%y(i,m,:))*st%y(j,m,:)) )
					st%bigL(i,m,j,m) = sys%gij(i,j) * (conjg(st%y0(i,m)) + st%y0(j,m))**sys%expnt
					st%bigL(j,m,i,m) = conjg( st%bigL(i,m,j,m) )
				end do
			end do

			!-- m != n
			do n=m+1, size(st%y,2)
				st%bigW(:,m,n) = matmul( conjg(st%y(:,m,:)) * st%y(:,n,:),  sys%wwk(:) )
				st%bigW(:,n,m) = conjg(st%bigW(:,m,n))
				do i=1, size(st%y,1)
					st%ovm(i,m,i,n) = overlap( st%y(i,m,:), st%y(i,n,:) )
					st%ovm(i,n,i,m) = conjg(st%ovm(i,m,i,n))
					!st%ovm(i,m,i,n) = st%ovm(i,m,i,n)*exp( sum( conjg(st%y(i,m,:))*st%y(i,n,:)) )
					!st%ovm(i,n,i,m) = conjg(st%ovm(i,m,i,n))
					st%bigL(i,m,i,n) = sys%gij(i,i) * ( conjg(st%y0(i,m)) + st%y0(i,n) )**sys%expnt
					st%bigL(i,n,i,m) = conjg(st%bigL(i,m,i,n))
				end do
				do i=1, size(st%y,1)
					do j=i+1, size(st%y,1)
						st%ovm(i,m,j,n) = overlap( st%y(i,m,:), st%y(j,n,:) )
						st%ovm(i,n,j,m) = overlap( st%y(i,n,:), st%y(j,m,:) )
						st%ovm(j,n,i,m) = conjg( st%ovm(i,m,j,n) )
						st%ovm(j,m,i,n) = conjg( st%ovm(i,n,j,m) )
						!st%ovm(i,m,j,n) = st%ovm(i,m,j,n)*exp( sum( conjg(st%y(i,m,:))*st%y(j,n,:)) )
						!st%ovm(i,n,j,m) = st%ovm(i,n,j,m)*exp( sum( conjg(st%y(i,n,:))*st%y(j,m,:)) )
						st%bigL(i,m,j,n) = sys%gij(i,j) * ( conjg(st%y0(i,m)) + st%y0(j,n) )**sys%expnt
						st%bigL(i,n,j,m) = sys%gij(i,j) * ( conjg(st%y0(i,n)) + st%y0(j,m) )**sys%expnt
						st%bigL(j,n,i,m) = conjg( st%bigL(i,m,j,n) )
						st%bigL(j,m,i,n) = conjg( st%bigL(i,n,j,m) )
					end do
				end do
			end do

		end do 

		call cpu_time(t2)
		sum_time = t2-t1

		return
	END SUBROUTINE

	FUNCTION overlap( f1, f2 )
		complex(8), intent(in)   ::  f1( : ), f2( : )
		complex(8)				 ::   overlap

		overlap = exp( sum( -0.5_8*conjg(f1)*f1 &
									- 0.5_8*conjg(f2)*f2 &
									+ conjg(f1)*f2 ) )
	END FUNCTION

	FUNCTION parameterchar(sys)
		type(param), intent(in)		::   sys
		character(len=20)      		:: ratio_char,nl_char,wc_char,ncsini_char,ncsmax_char,&
			errorthr_char,dtadd_char,addingratio_char, wq_char,gqc_char, wge_char,&
			alpha_char,nmodes_char,p0_char, qbini_char, cavini_char,  &
			wd_char, ad_char, offset_abs_char, offset_ang_char, anh_char, dt_char,&
			errlim_char, device_char,tmax_char, moderatio_char, p2threshold_char, pid_char, &
			dtadd_ini_char, n_dtadd_change_char, lsfe_char,&
			bw_char, bc_cav_char, bc_qb_char, bc_lf_char, bc_hf_char,&
			moderatio_qb_char, moderatio_lf_char, moderatio_hf_char


		character(len=300)				:: parameterchar, addchar

		write( nl_char, '(I3)' ) sys%nl
		write( wc_char, '(f6.3)' ) sys%wc/twopi
		write( ncsini_char, '(I3)' ) sys%ncs_ini
		write( ncsmax_char, '(I3)' ) sys%ncs_max
		write( errorthr_char, '(e10.2)' ) sys%error_thr
		write( dtadd_char, '(f6.2)' ) sys%dt_add
		write( dtadd_ini_char, '(f6.2)' ) sys%dt_add_ini
		write( n_dtadd_change_char, '(I2)' ) sys%n_dtadd_change
		write( addingratio_char, '(f6.1)' ) sys%adding_ratio
		write( wge_char, '(f8.4)' ) sys%w_ge/twopi
		write( gqc_char, '(f8.3)' ) sys%g_qc/twopi
		write( alpha_char, '(f10.6)' ) sys%alpha/twopi
		write( nmodes_char, '(I4)' ) sys%nbathmodes
		write( p0_char, '(e10.1)' ) sys%p0
		write( bw_char, '(f6.3)' ) sys%bw/twopi
		write( bc_cav_char, '(f8.3)' ) sys%bc_cav/twopi
		write( bc_qb_char, '(f8.3)' ) sys%bc_qb/twopi
		write( bc_lf_char, '(f8.3)' ) sys%bc_lf/twopi
		write( bc_hf_char, '(f8.3)' ) sys%bc_hf/twopi
		write( wd_char, '(f8.3)' ) sys%wd/twopi
		write( ad_char, '(f7.4)' ) sys%ad/twopi
		write( offset_abs_char, '(f6.2)' ) sys%offset_abs
		write( offset_ang_char, '(f6.2)' ) sys%offset_ang
		write( anh_char, '(f9.3)' ) sys%anh/twopi
		write( dt_char, '(f9.4)' ) sys%dt
		write( errlim_char, '(e10.1)' ) sys%err_lim
		write( qbini_char, '(I2)' ) sys%qb_ini
		write( cavini_char, '(f6.1)' ) sys%cav_ini
		write( tmax_char, '(I5)' ) int( sys%tmax )
		write( moderatio_qb_char, '(f4.2)' ) sys%mode_ratio_qb
		write( moderatio_lf_char, '(f4.2)' ) sys%mode_ratio_lf
		write( moderatio_hf_char, '(f4.2)' ) sys%mode_ratio_hf
		write( p2threshold_char, '(e10.1)' ) sys%p2_threshold
		write( pid_char, '(I10)' ) sys%pid
		write( lsfe_char, '(I3)' ) sys%lim_slow_fact_exp

		parameterchar="nl"//trim(adjustl(nl_char))//&
				"_n"//trim(adjustl(ncsini_char))//"_"//trim(adjustl(ncsmax_char))//&
			"_E"//trim(adjustl(errorthr_char))//&
			!"_p2th"//trim(adjustl(p2threshold_char))//&
			"_dtadd"//trim(adjustl(dtadd_char))//&
			"_"//trim(adjustl(dtadd_ini_char))//&
			"_"//trim(adjustl(n_dtadd_change_char))//&
			"_ar"//trim(adjustl(addingratio_char))//&
			"_g"//trim(adjustl(gqc_char))//&
			"_al"//trim(adjustl(alpha_char))//&
			"_nm"//trim(adjustl(nmodes_char))//&
			"_p"//trim(adjustl(p0_char))//&
			"_wq"//trim(adjustl(wge_char))//&
			"_wc"//trim(adjustl(wc_char))//&
			"_bc"//trim(adjustl(bc_lf_char))//&
			"_"//trim(adjustl(bc_qb_char))//&
			"_"//trim(adjustl(bc_cav_char))//&
			"_"//trim(adjustl(bc_hf_char))//&
			"_bw"//trim(adjustl(bw_char))//sys%bath_type//&
			"_mr"//trim(adjustl(moderatio_lf_char))//&
			"_"//trim(adjustl(moderatio_qb_char))//&
			"_"//trim(adjustl(moderatio_hf_char))//&
			"_wd"//trim(adjustl(wd_char))//&
			"_ad"//trim(adjustl(ad_char))//&
			!"_OFT"//trim(adjustl(offset_abs_char))//'_'//trim(adjustl(offset_ang_char))//&
			"_anh"//trim(adjustl(anh_char))//&
			"_dt"//trim(adjustl(dt_char))//&
			"_errl"//trim(adjustl(errlim_char))//&
			"_lsfe"//trim(adjustl(lsfe_char))//&
			"_qb"//trim(adjustl(qbini_char))//&
			"_cv"//trim(adjustl(cavini_char))//&
			"_tmax"//trim(adjustl(tmax_char))//&
			"_"//sys%dvice
			!"_"//trim(adjustl(sys%method))//&
	END FUNCTION

	FUNCTION norm(st)
		real(rl) 			       		::  norm
		complex(cx)						::  tmp
		type(state),intent(in)  		::  st
		integer						 	::  i, m,n

		tmp = 0._cx
		do i=1,size(st%p,1)
			do m=1,size(st%p,2)
				do n=1,size(st%p,2)
					tmp = tmp + conjg(st%p(i,m))*st%p(i,n)*st%ovm(i,m,i,n)
				end do
			end do
		end do

		norm = sqrt( real( tmp ) )
	END FUNCTION

	SUBROUTINE normalise(st)
		type(state), intent(in out)  	::  st
		st%p = st%p/norm(st)

		RETURN
	END SUBROUTINE

	FUNCTION energy( sys, st )
		type(param), intent(in)			::   sys
		type(state), intent(in)			::   st
		real(rl)						::   energy
		complex(cx)						::   tmp
		integer							::   i, m, n, l

		tmp = 0._rl
		do i=1, sys%nl
			do m=1, st%ncs
				do n=1, st%ncs

					tmp = tmp + conjg(st%p(i,m))*st%p(i,n)*st%ovm(i,m,i,n) &
							*( sys%w_qb(i) + st%bigW(i,m,n) + sys%ad*dcos(sys%wd*st%t)*st%bigU(i,m,n) )

					do l=1, sys%nl
						tmp = tmp + conjg(st%p(l,m))*st%p(i,n)*st%ovm(l,m,i,n)*st%bigL(l,m,i,n)
					end do
				end do
			end do
		end do

		energy = real( tmp )
	END FUNCTION

	FUNCTION lvl_occupation( st, s )
		type(state), intent(in)		::   st
		integer, intent(in)			::   s
		real(rl)					::   lvl_occupation
		complex(cx)					::   tmp
		integer						::   m, n

		tmp = 0.0_rl
		do m=1, st%ncs
			do n=1, st%ncs
				tmp = tmp + conjg(st%p(s,m))*st%p(s,n)*st%ovm( s,m,s,n )
			end do
		end do

		lvl_occupation = real( tmp )
	END FUNCTION

	FUNCTION ds_occupation( sys, st, s )
		type(param), intent(in)		::   sys
		type(state), intent(in)		::   st
		type(state)					::   ds
		integer, intent(in)			::   s
		real(rl)					::   ds_occupation
		complex(cx)					::   tmp
		integer						::   m, n, i ,j

		if (s==1) then
			ds = sys%ground_state
		else if (s==2) then
			ds = sys%excited_state
		end if

		tmp = 0.0_rl
		do i=1, sys%nl
			do m=1, ds%ncs
				do n=1, st%ncs
					tmp = tmp + conjg(ds%p(i,m))*st%p(i,n)*ov( ds%y(i,m,:), st%y(i,n,:) )
				end do
			end do
		end do

		ds_occupation  = real( conjg( tmp )*tmp )
	END FUNCTION

	FUNCTION photon_nb_in_sig_mode( sys, st, sig )
		type(param), intent(in)			::   sys
		type(state), intent(in)			::   st
		integer, intent(in)				::   sig
		real(rl)						::   photon_nb_in_sig_mode
		complex(cx)						::   tmp
		complex(cx), dimension(st%ncs)	::   y_OB
		integer							::   i, m, n

		y_OB = 0._rl
		tmp = 0._rl

		do i=1, sys%nl

			y_OB = 0._rl
			do n=1, st%ncs
				y_OB(n) = sum( sys%Omat(sig,:) * st%y(i,n,:) )
			end do

			do m=1, st%ncs
				do n=1, st%ncs
					tmp = tmp + conjg(st%p(i,m))*st%p(i,n)*conjg(y_OB(m))*y_OB(n)*st%ovm( i,m,i,n )
				end do
			end do

		end do

		photon_nb_in_sig_mode = real( tmp )
	END FUNCTION

	FUNCTION apadag_sig_mode( sys, st, sig )
		type(param), intent(in)			::   sys
		type(state), intent(in)			::   st
		integer, intent(in)				::   sig
		real(rl)						::   apadag_sig_mode
		complex(cx)						::   tmp
		complex(cx), dimension(st%ncs)	::   y_OB
		integer							::   i, m, n

		y_OB = 0._rl
		tmp = 0._rl

		do i=1, sys%nl

			y_OB = 0._rl
			do n=1, st%ncs
				y_OB(n) = sum( sys%Omat(sig,:) * st%y(i,n,:) )
			end do

			do m=1, st%ncs
				do n=1, st%ncs
					tmp = tmp + conjg(st%p(i,m))*st%p(i,n)*( conjg(y_OB(m)) + y_OB(n) )*st%ovm( i,m,i,n )
				end do
			end do

		end do

		apadag_sig_mode = real( tmp )
	END FUNCTION

	FUNCTION renyi_entropy( sys, st, s, n_max )
		type(param), intent(in)					::   sys
		type(state), intent(in)					::   st
		integer, intent(in)						::   s, n_max
		real(rl)								::   renyi_entropy
		complex(cx)								::   tmp, tmp_nm, num
		integer									::   i, j, l1, l2, n, m
		complex(cx), dimension(sys%nl,st%ncs)	::   y_OB

		y_OB = 0._rl
		do i=1, sys%nl
			do n=1, st%ncs
				y_OB(i,n) = sum( sys%Omat(s,:) * st%y(i,n,:) )
			end do
		end do

		renyi_entropy  = 1.0

		do i=1, sys%nl
			do j=1, sys%nl
				do l1=0, n_max
					do l2=0, n_max

						tmp_nm = 0.0
						do n=1, st%ncs
							do m=1, st%ncs

								num = conjg(st%p(i,n))*st%p(j,m)*st%ovm( i,n,j,m ) &
									* exp( -conjg(y_OB(i,n))*y_OB(j,m) ) &
									* conjg(y_OB(i,n))**l1 * y_OB(j,m)**l2

								tmp_nm = tmp_nm + num /  dsqrt( dble(factorial(l1)*factorial(l2)) )

							end do
						end do
	
						renyi_entropy = renyi_entropy - conjg( tmp_nm ) * tmp_nm

					end do
				end do
			end do
		end do

		renyi_entropy = real( renyi_entropy )
	END FUNCTION

	!-- Calculate the overlap between two coherent states
	FUNCTION ov_scalar(f1,f2)
		complex(cx), intent(in) :: f1, f2
		complex(cx)             :: ov_scalar
		complex(cx)				  :: tmp1, tmp2, tmp3

		tmp1 = conjg(f1)*f1
		tmp2 = conjg(f2)*f2
		tmp3 = conjg(f1)*f2

		ov_scalar = exp( -0.5_rl*tmp1 - 0.5_rl*tmp2 + tmp3 ) 
	END FUNCTION ov_scalar

	FUNCTION ov(f1,f2)
		complex(cx), intent(in) :: f1( : ), f2( : )
		complex(cx)             :: ov
		!= internal variables
		complex(cx)    :: tmp1, tmp2, tmp3

		!= initialize
		tmp1 = 0._rl
		tmp2 = 0._rl
		tmp3 = 0._rl
		if (size(f1,1) .ne. 1) then
			tmp1 = dot_product(f1, f1)
			tmp2 = dot_product(f2, f2)
			tmp3 = dot_product(f1, f2)
		else if (size(f1,1) == 1) then
			tmp1 = conjg(f1(1))*f1(1)
			tmp2 = conjg(f2(1))*f2(1)
			tmp3 = conjg(f1(1))*f2(1)
		end if

		ov = exp( -0.5_rl*tmp1 - 0.5_rl*tmp2 + tmp3 ) 
	END FUNCTION ov

	!======================================================
	!== MATH FUNCTIONS
	!======================================================

	FUNCTION factorial(n)
		integer, intent(in)		::  n
		integer						::  factorial
		integer						::  tmp, i 

		factorial=1
		if (n > 1) then
			do i=1,n
				factorial = factorial*i
			end do
		end if
	END FUNCTION

	SUBROUTINE print_log(sys, message, print_to, formt)
		type(param), intent(in)					::     sys
		character(len=*), intent(in)			::     message
		integer, intent(in)						::     print_to
		character(len=*), intent(in), optional	::     formt

		if (print_to == 0) then
			!-- write to console
			write(*,*) message
		elseif (print_to  == 1) then
			!-- write to file
			open(unit=100,file='data/LOG_'//trim(adjustl(parameterchar(sys)))//'.d' ,&
				action="write",position='append',status="old")
			write(100,*) message
			close(100)
		elseif (print_to == 2) then
			!-- write to both file and console
			write(*,*) message
			open(unit=100,file='data/LOG_'//trim(adjustl(parameterchar(sys)))//'.d' ,&
				action="write",position='append',status="old")
			write(100,*) message
			close(100)
		end if
	END SUBROUTINE

END MODULE SYSTM

!	SUBROUTINE calculate_derivatives( sys, st )
!		type(param),intent(in)  		::  sys
!		type(state),intent(in out)  	::  st
!
!		CALL calc_derivatives( sys%w_qb, sys%wd, sys%Ad, st%t, sys%wwk,&
!			sys%gk, st%ovm, st%bigW, st%bigL, st%bigU, sys%Omat, st%p, st%y, st%pdot, st%ydot )
!		
!		RETURN
!	END SUBROUTINE

!	SUBROUTINE initialise_from_file(sys,st)
!
!		type(param), intent(in)			  	::  sys
!		type(state), intent(in out)	  		::  st
!		integer  									::  i,j,m,k
!		character(len=200)						::  fks_file,ps_file
!		real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a
!
!		print*, "Initialising from: ", parameterchar(sys)
!		fks_file="data/fks_fst_"//trim(adjustl(parameterchar(sys)))//".d"
!		ps_file="data/ps_fst_"//trim(adjustl(parameterchar(sys)))//".d"
!		open (unit=101,file=ps_file,action="read",status="old")
!		open (unit=100,file=fks_file,action="read",status="old")
!
!		do  k=1,sys%nmode
!			read(100,'(f25.15)',advance='no') a
!			do i=1,st%np
!				read(100,'(2f25.15)',advance='no') f_r, h_r
!				st%f(i,k) = f_r
!				st%h(i,k) = h_r
!			end do
!			do i=1,st%np
!				read(100,'(2f25.15)',advance='no') f_i, h_i
!				st%f(i,k) = st%f(i,k) + Ic*f_i
!				st%h(i,k) = st%h(i,k) + Ic*h_i
!			end do
!			read(100,*)
!		end do
!
!		do i=1,st%np
!			read(101,'(2f25.15)',advance='no') p_r, q_r
!			st%p(i) = p_r
!			st%q(i) = q_r
!		end do
!		do i=1,st%np
!			read(101,'(2f25.15)',advance='no') p_i, q_i
!			st%p(i) = st%p(i) + Ic*p_i
!			st%q(i) = st%q(i) + Ic*q_i
!		end do
!
!		st%t = sys%tmax
!		close(100)
!		close(101)
!
!		!-- updating the sums over k
!		CALL update_sums(sys,st)
!		CALL normalise(st)
!
!	END SUBROUTINE
