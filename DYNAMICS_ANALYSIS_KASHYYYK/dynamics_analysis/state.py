import numpy as np
import matplotlib.pyplot as plt
from .fortran_module import fortran_module

def KD( a, b ):

	if a==b:
		return 1
	else:
		return 0

def ov( f1, f2 ):

	overlap = np.exp( -0.5*np.vdot(f1,f1) - 0.5*np.vdot(f2,f2) + np.vdot(f1,f2) )
	return overlap

class state(object):

	def __init__( self, *args ):
		
		if len( args ) == 2:
			ps = args[0]
			yks = args[1]
			self.p = ps
			self.y = yks
			self.ncs = len( self.p[0,:] )
			self.nl = len( self.p[:,0] )
			self.nmodes = len( self.y[0,0,:] )

			self.ovm = np.zeros( ( self.nl, self.ncs,  self.nl, self.ncs), dtype='complex128' )
			self.update_sums()

		elif len( args ) == 4:
			self.nl = args[0]
			self.nmodes = args[1] 
			path = args[2]
			self.paramchar = args[3]

			ps_data_re_im = np.loadtxt( path+"PS_"+self.paramchar+'.d' )
			yks_data_re_im = np.loadtxt( path+"YKS_"+self.paramchar+'.d' )
			self.ncs = ( len( ps_data_re_im[:] )//2 ) // self.nl  

			ps_data_complex = ps_data_re_im[ :self.ncs*self.nl ] + 1j* ps_data_re_im[ self.ncs*self.nl: ]
			yks_data_complex = yks_data_re_im[:, 1:self.ncs*self.nl+1 ] \
					+ 1j* yks_data_re_im[ : , self.ncs*self.nl+1: ]

			self.p = np.zeros( (self.nl,self.ncs), dtype='complex128' )
			self.y = np.zeros( (self.nl, self.ncs, self.nmodes), dtype='complex128' )
			for i in range( self.nl ):
				for n in range( self.ncs ):
					self.p[i,n] = ps_data_complex[ i*self.ncs + n ]
					for k in range( self.nmodes ):
						self.y[i,n,k] = yks_data_complex[ k, i*self.ncs + n ]

		self.ovm = np.zeros( ( self.nl, self.ncs,  self.nl, self.ncs), dtype='complex128' )
		self.update_sums()

	def update_sums( self ):

		for m in range( self.ncs ):
			for n in range( self.ncs ):
				for i in range( self.nl ):
					for j in range( self.nl ):
						self.ovm[i,m,j,n] = ov( self.y[i,m,:], self.y[j,n,:] )

	def norm( self ):

		tmp = 0
		for i in range( len( self.p[:,0] ) ):
			for m in range( self.ncs ):
				for n in range( self.ncs ):
					tmp += np.conj(self.p[i,m])*self.p[i,n]*self.ovm[ i,m,i,n ]

		return np.sqrt(np.real(tmp))

	def norm2( self ):

		tmp = 0
		for i in range( len( self.p[:,0] ) ):
			for m in range( self.ncs ):
				for n in range( self.ncs ):
					tmp += np.conj(self.p[i,m])*self.p[i,n]*self.ovm[ i,m,i,n ]

		return np.real(tmp)

	def normalize( self ):

		self.p /= self.norm()

	def energy( self, prms ):

		p = self.p
		pc = np.conj(p)
		y = self.y
		yc = np.conj(y)
		ovm = self.ovm

		tmp = 0
		for i in range( prms.nl ):
			for m in range( self.ncs ):
				for n in range( self.ncs ):

					tmp +=  pc[i,m]*p[i,n]*ovm[i,m,i,n]\
							*( prms.w_qb[i] + self.bigW[i,m,n] \
							+ prms.A_d*np.cos(prms.wd*self.t)*self.bigU[i,m,n] )

					for l in range( prms.nl ):

						tmp += pc[l,m]*p[i,n]*ovm[l,m,i,n]*self.bigL[l,m,i,n]

		return tmp 

	def renyi_entropy_without_decimal( self, s ):

		from decimal import Decimal
		from math import factorial
		lmax = 13

		tmp = 1.0
		
		for i in range( self.nl ):
			for j in range( self.nl ):
				for l1 in range( lmax ):
					for l2 in range( lmax ):

						tmp_nm = 0.0
						for n in range( self.ncs ):
							for m in range( self.ncs ):

								num = np.conj(self.p[i,n])*self.p[j,m]*self.ovm[ i,n,j,m ]*np.exp( -np.conj(self.y[i,n,s])*self.y[j,m,s] ) * np.conj(self.y[i,n,s])**l1 * self.y[j,m,s]**l2

								tmp_nm += num / np.sqrt( factorial(l1)*factorial(l2) ) 
	
						tmp -= np.abs( tmp_nm )**2

		return tmp

	def renyi_entropy( self, s, n_max ):

		from decimal import Decimal
		from math import factorial

		tmp = 1.0
		
		for i in range( self.nl ):
			for j in range( self.nl ):
				for l1 in range( n_max+1 ):
					for l2 in range( n_max+1 ):

						tmp_nm = 0.0
						for n in range( self.ncs ):
							for m in range( self.ncs ):

								num = np.conj(self.p[i,n])*self.p[j,m]*self.ovm[ i,n,j,m ]*np.exp( -np.conj(self.y[i,n,s])*self.y[j,m,s] ) * np.conj(self.y[i,n,s])**l1 * self.y[j,m,s]**l2

								tmp_nm += num / float( Decimal( factorial(l1)*factorial(l2) ).sqrt() ) 
	
						tmp -= np.abs( tmp_nm )**2

		return tmp

	def renyi_entropy_from_db( self, s, n_max, omat ):

		from decimal import Decimal
		from math import factorial

		y_OB = np.zeros( (self.nl,self.ncs), dtype='complex128' )
		for i in range(self.nl):
			for n in range(self.ncs):
				y_OB[i,n] = np.sum( omat[s,:]*self.y[i,n,:] )

		entropy = 1.0
		for i in range( self.nl ):
			for j in range( self.nl ):
				for l1 in range( n_max+1 ):
					for l2 in range( n_max+1 ):

						tmp_nm = 0.0
						for n in range( self.ncs ):
							for m in range( self.ncs ):

								num = np.conj(self.p[i,n])*self.p[j,m]*self.ovm[ i,n,j,m ] \
										* np.exp( -np.conj(y_OB[i,n])*y_OB[j,m] ) \
										* np.conj(y_OB[i,n])**l1 * y_OB[j,m]**l2

								tmp_nm = tmp_nm + num / np.sqrt( factorial(l1)*factorial(l2) ) 
	
						entropy = entropy - np.conj( tmp_nm )*tmp_nm

		return entropy

	def lvl_occupation( self, s ):

		tmp = 0.0
		for m in range( self.ncs ):
			for n in range( self.ncs ):
				tmp += np.conj(self.p[s,m])*self.p[s,n]*self.ovm[ s,m,s,n ]

		return np.real(tmp)

	def destroy_k( self, k ):

		for i in range( self.nl ):
			for m in range( self.ncs ):
				self.p[i,m] *= self.y[i,m,k] 

		self.normalize()

	def ph_nb_in_k( self, k ):

		p = self.p
		pc = np.conj( p )
		y = self.y
		yc = np.conj( y )

		tmp = 0
		nl = len( p[:,0] )
		for i in range( nl ):
			for m in range( self.ncs ):
				for n in range( self.ncs ):
					tmp += pc[i,m]*p[i,n]*yc[i,m,k]*y[i,n,k]*self.ovm[ i,m,i,n ]

		if np.imag(tmp) > 1e-8:
			print('ERROR in norm: imaginary of Y is non-zero')

		return np.real(tmp)

	def tot_photon_k_nb( self ):

		tot_photon_k_nb  = np.zeros( self.nmodes, dtype='complex128' )
		p = self.p
		pc = np.conj( p )
		y = self.y
		yc = np.conj( y )

		for k in range( self.nmodes ):
			tmp = 0
			for i in range( self.nl ):
				for m in range( self.ncs ):
					for n in range( self.ncs ):
						tot_photon_k_nb[k] += pc[i,m]*p[i,n]*yc[i,m,k]*y[i,n,k]*self.ovm[ i,m,i,n ]

		return np.real( tot_photon_k_nb )

	def one_photon_k_nb( self ):

		one_photon_k_nb = np.zeros( self.nmodes, dtype='float64' )
		zeros_arr = np.zeros( self.nmodes, dtype='complex128' )
		ov_0y = np.zeros( self.ncs, dtype='complex128'  )

		for k in range(self.nmodes):

			tmp =0.0
			for s in range(self.nl):
				for n in range(self.ncs):
					ov_0y[n] = ov( zeros_arr , self.y[s,n,:] )

				tmp +=  np.abs( sum( self.p[s,:] * self.y[s,:,k] * ov_0y[:] ) ) **2
				
			one_photon_k_nb[k] = tmp

		return one_photon_k_nb

	def two_photon_k_nb( self ):
		
		two_photon_alpha_kk  = np.zeros( (self.nmodes,  self.nmodes), dtype='complex128' )
		two_photon_k_nb  = np.zeros( self.nmodes, dtype='float64' )
		zeros_arr = np.zeros( self.nmodes, dtype='complex128' )
		ov_0y = np.zeros( self.ncs, dtype='complex128'  )

		for s in range(self.nl):

			for n in range(self.ncs):
				ov_0y[n] = ov( zeros_arr , self.y[s,n,:] )

			for k1 in range(self.nmodes):
				for k2 in range(k1, self.nmodes):
					two_photon_alpha_kk[k1,k2] = sum( self.p[s,:] * self.y[s,:,k1] * self.y[s,:,k2] * ov_0y[:] ) / 2
					two_photon_alpha_kk[k2,k1] = two_photon_alpha_kk[k1,k2]

			for k1 in range(self.nmodes):
				for k2 in range(self.nmodes):
						two_photon_k_nb[k1] += 4*np.abs(two_photon_alpha_kk[k1,k2])**2

		return two_photon_k_nb

	def three_photon_k_nb( self ):

		three_photon_alpha_kkk  = np.zeros( (self.nmodes,self.nmodes,self.nmodes), dtype='complex128' )
		three_photon_k_nb  = np.zeros( self.nmodes, dtype='float64' )
		zeros_arr = np.zeros( self.nmodes, dtype='complex128' )
		ov_0y = np.zeros( self.ncs, dtype='complex128'  )

		for s in range(self.nl):

			for n in range(self.ncs):
				ov_0y[n] = ov( zeros_arr , self.y[s,n,:] )


			for k1 in range(self.nmodes):
				for k2 in range(k1,self.nmodes):
					for k3 in range(k2,self.nmodes):

						three_photon_alpha_kkk[k1,k2,k3] = sum( self.p[s,:] * self.y[s,:,k1]\
								* self.y[s,:,k2] * self.y[s,:,k3] * ov_0y[:] ) / 6
						three_photon_alpha_kkk[k1,k3,k2] = three_photon_alpha_kkk[k1,k2,k3]
						three_photon_alpha_kkk[k2,k1,k3] = three_photon_alpha_kkk[k1,k2,k3]
						three_photon_alpha_kkk[k2,k3,k1] = three_photon_alpha_kkk[k1,k2,k3]
						three_photon_alpha_kkk[k3,k1,k2] = three_photon_alpha_kkk[k1,k2,k3]
						three_photon_alpha_kkk[k3,k2,k1] = three_photon_alpha_kkk[k1,k2,k3]

			for k1 in range(self.nmodes):
				for k2 in range(self.nmodes):
					for k3 in range(self.nmodes):
						three_photon_k_nb[k1] += np.abs( three_photon_alpha_kkk[k1,k2,k3] )**2

		return 18*three_photon_k_nb

	def tot_ph_nb( self ):

		p = self.p
		pc = np.conj( p )
		y = self.y
		yc = np.conj( y )

		tmp = 0
		nl = len( p[:,0] )
		for i in range( nl ):
			for m in range( self.ncs ):
				for n in range( self.ncs ):
					tmp += pc[i,m]*p[i,n]*np.sum(yc[i,m,:]*y[i,n,:])*self.ovm[ i,m,i,n ]

		if np.imag(tmp) > 1e-8:
			print('ERROR in norm: imaginary of Y is non-zero')

		return np.real(tmp)

	def plot_split_wigners(self,  xmin, log_min=None, add_cs_centers=False ):

		from matplotlib.colors import LogNorm, LinearSegmentedColormap

		cmap = LinearSegmentedColormap.from_list("", ["white",'yellow','orange','red',\
				'black','royalblue','cornflowerblue','lightsteelblue', 'lightblue'])

		xmax=-xmin
		figure, axis = plt.subplots(1,4,figsize=(16,9),gridspec_kw={'width_ratios': [1,0.05, 1,0.05]}) 

		for s in range(2):
			
			axis[s*2].set_title('state=|'+str(s)+'>')
			wigner_ = fortran_module.calc_split_wigner( s+1, self.p[:,:], \
				self.y[:,:,0], self.ovm, xmin, xmax, 100 )
			axis[s*2].axhline( y=0, dashes=[2,2,2,2] )
			axis[s*2].axvline( x=0, dashes=[2,2,2,2] ) 

			if  add_cs_centers:
				for n in range(self.ncs):
					size = 1 + np.abs(self.p[s,n])*50
					axis[s*2].scatter( np.real(self.y[s,n,0]), np.imag(self.y[s,n,0]), s=size, color='C'+str(2+s) )

			if log_min:
				norm_wigner=LogNorm(vmin=log_min,vmax =wigner_.max()  )
				im=axis[s*2].imshow( wigner_, extent=[xmin,xmax,xmin,xmax],norm=norm_wigner,origin ='lower' )
			else:
				abs_max =  np.maximum( wigner_.max(), np.abs(wigner_.min()) ) 
				norm_wigner = plt.Normalize(-abs_max, abs_max)
				im=axis[s*2].imshow( wigner_, extent=[xmin,xmax,xmin,xmax],origin ='lower',\
						norm=norm_wigner, cmap=cmap )
			figure.colorbar(im, cax=axis[2*s+1], orientation='vertical')

		plt.show()

	def plot_wigner(self,  xmin, log_min=None ):#, add_cs_centers=False ):

		from matplotlib.colors import LogNorm, LinearSegmentedColormap

		cmap = LinearSegmentedColormap.from_list("", ["white",'yellow','orange','red',\
				'black','royalblue','cornflowerblue','lightsteelblue', 'lightblue'])

		xmax=-xmin
		figure, ax = plt.subplots(1,1) 

		wigner_ = fortran_module.calc_wigner( self.p[:,:], \
				self.y[:,:,0], self.ovm, xmin, xmax, 100 )
		ax.axhline( y=0, dashes=[2,2,2,2] )
		ax.axvline( x=0, dashes=[2,2,2,2] ) 


		#if  add_cs_centers:
		#	for n in range(self.ncs):
		#		size = 1 + np.abs(self.p[s,n])*50
		#		axis[s*2].scatter( np.real(self.y[s,n,0]), np.imag(self.y[s,n,0]), s=size, color='C'+str(2+s) )

		if log_min:
			norm_wigner=LogNorm(vmin=log_min,vmax =wigner_.max()  )
			im=ax.imshow( wigner_, extent=[xmin,xmax,xmin,xmax],norm=norm_wigner,origin ='lower' )
		else:
			abs_max =  np.maximum( wigner_.max(), np.abs(wigner_.min()) ) 
			norm_wigner = plt.Normalize(-abs_max, abs_max)
			im=ax.imshow( wigner_, extent=[xmin,xmax,xmin,xmax],origin ='lower',\
					norm=norm_wigner, cmap=cmap )
		figure.colorbar(im, cax=ax, orientation='vertical')

		plt.show()


