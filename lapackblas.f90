MODULE lapackblas

	USE typedefs, only : cx => c_type, rl => r_type

	IMPLICIT NONE

CONTAINS

	SUBROUTINE diagonalise_all(inmat, evec, eval)
		!= calculate the eigen values and eigen vectors
		!= of a real symmetric matrix and right eigenvectros.
		!= returns also eigen values and eigenvectors.

		real(rl), intent(in)           			::	 inmat(:,:)		
		real(rl), intent(out) 					::	 evec(size(inmat,1),size(inmat,1))
		real(rl), intent(out) 					::	 eval(size(inmat,1))
		real(rl), allocatable             		::	 work(:)
		integer                               	::	 lwork, info, i

		allocate(work(size(inmat,1)))
		!= first query for optimal work space
		lwork = -1
		call dsyev( 'V', 'U', size(inmat,1), inmat, size(inmat,1),&
			eval, work, lwork, info )		

		lwork =  int(work(1))

		deallocate(work)
		allocate(work(lwork))

		!= now call the routine again for the final job
		call dsyev( 'V', 'U', size(inmat,1), inmat, size(inmat,1),&
			eval, work, lwork, info )		

		evec = inmat

		deallocate(work)
		return
	END SUBROUTINE diagonalise_all

	SUBROUTINE diagonalise_hermitian_zhpev(inmat, eval)
		!= calculate the eigen values and eigen vectors
		!= of a real symmetric matrix and right eigenvectros.
		!= returns also eigen values and eigenvectors.

		complex(cx), intent(in)           			::	 inmat(:,:)
		complex(cx)          						::	 inmat_duplicate(size(inmat,1),size(inmat,1))
		complex(cx), intent(out) 					::	 eval(size(inmat,1))
		complex(cx) 								::	 evec(size(inmat,1),size(inmat,1))
		real(rl), allocatable             			::	 rwork(:)
		complex(cx), allocatable             		::	 work(:)
		integer                               		::	 info, i, LDA

		inmat_duplicate = inmat

		LDA = size(inmat,1)
		allocate(work(LDA))
		allocate(rwork(LDA))

		!= now call the routine again for the final job
		call ZHPEV( 'N', 'U', size(inmat,1), inmat_duplicate, eval, evec, size(inmat,1), work, rwork, info )

		deallocate(work)
		return
	END SUBROUTINE diagonalise_hermitian_zhpev

	SUBROUTINE diagonalise_complex_zgeev(inmat, eval)
		complex(cx), intent(in)           			::	 inmat(:,:)
		complex(cx)          						::	 inmat_duplicate(size(inmat,1),size(inmat,1))
		complex(cx)          						::	 inmat_duplicate_2(size(inmat,1),size(inmat,1))
		complex(cx), intent(out) 					::	 eval(size(inmat,1))
		complex(cx)				 					::	 evec_left(size(inmat,1),size(inmat,1))
		complex(cx)				 					::	 evec_right(size(inmat,1),size(inmat,1))
		real(rl), allocatable             			::	 rwork(:)
		complex(cx), allocatable             		::	 work(:)
		integer                               		::	 lwork, info, i, LDA

		inmat_duplicate = inmat
		inmat_duplicate_2 = inmat

		LDA = size(inmat,1)
		allocate(rwork(LDA))

		allocate(work(size(inmat,1)))
		!= first query for optimal work space
		lwork = -1
		call zgeev( 'N', 'N', LDA, inmat_duplicate, LDA,&
			eval, evec_left, LDA, evec_right, LDA,&
			work, lwork, rwork, info )

		lwork =  int(work(1))

		deallocate(work)
		allocate(work(lwork))

		!= now call the routine again for the final job
		call zgeev( 'N', 'N', LDA, inmat_duplicate_2, LDA,&
			eval, evec_left, LDA, evec_right, LDA,&
			work, lwork, rwork, info )

		deallocate(work)
		return
	END SUBROUTINE diagonalise_complex_zgeev


END MODULE lapackblas
