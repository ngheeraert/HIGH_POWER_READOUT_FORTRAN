MODULE fortran_module

	implicit none

	real(8),   parameter :: pi=(4._8)*atan(1._8)
	real(8),   parameter :: twopi=2.0_8*(4.0_8)*atan(1.0_8)
	real(8),   parameter :: one=1.0_8
	real(8),   parameter :: zero=0.0_8
	complex(8),parameter :: ic=(0.0_8,1.0_8)


CONTAINS

	SUBROUTINE calc_wigner( p_in, f_in, ovm, xmin, xmax,  wigner, xnum )
		complex(8), intent(in)                          ::  p_in(:,:)
		complex(8), intent(in)                          ::  f_in(:,:)
		integer, intent(in)							    ::  xnum
		real(8), intent(in)							    ::  xmin, xmax
		complex(8), intent(in)                          ::  ovm(:,:,:,:)
		real(8)							                ::  dx, x, p
		complex(8)					  					::  tmp, zz
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

				wigner( xj, xi ) = (2._8/pi)*tmp

			end do
		end do

		return
	END SUBROUTINE

	SUBROUTINE calc_split_wigner( i,  p_in, f_in, ovm, xmin, xmax, xnum, wigner )
		complex(8), intent(in)                          ::  p_in(:,:)
		complex(8), intent(in)                          ::  f_in(:,:)
		integer, intent(in)							    ::  xnum
		real(8), intent(in)							    ::  xmin, xmax
		complex(8), intent(in)                          ::  ovm(:,:,:,:)
		integer, intent(in)							  	::  i
		real(8)							                ::  dx, x, p
		complex(8)					  					::  tmp, zz
		integer							  				::  n,m,xi,xj,ncs,nl
		real(8), intent(out)				            ::  wigner(xnum,xnum)

		nl = size( p_in, 1 )
		ncs = size( p_in, 2 )
		dx =   (xmax - xmin) / dble(xnum-1) 
		wigner = 0._8

		do xi=1, xnum
			do xj=1, xnum

				tmp = 0._8

				x = xmin + dx*(xi-1)
				p = xmin + dx*(xj-1)

				do n=1, ncs
					do m=1, ncs

						zz =  x + Ic*p 

						tmp = tmp + conjg(p_in(i,n))*p_in(i,m) &
							* exp( -2._8 * ( real(zz)+Ic*aimag(zz) - f_in(i,m) ) &
							* ( real(zz)-Ic*aimag(zz) - conjg( f_in(i,n) ) ) ) &
							* ovm(i,n,i,m)

					end do
				end do

				wigner( xj, xi ) = (2._8/pi)*tmp

			end do
		end do

		return
	END SUBROUTINE

END MODULE


