PROGRAM MAIN

	use systm, only: param, get_parameters, parameterchar
	use dynamics, only: compute_trajectory

	IMPLICIT NONE

	type(param)			::    sys

  	CALL get_parameters( sys )
  	CALL compute_trajectory( sys )

END PROGRAM

