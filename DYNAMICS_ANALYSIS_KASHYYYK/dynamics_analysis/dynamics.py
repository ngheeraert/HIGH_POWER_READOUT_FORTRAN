import numpy as np
import matplotlib.pyplot as plt
from .state import state
from .fortran_module import fortran_module
from os.path import exists


class dynamics(object):

	def __init__( self, nl, path, paramchar ):

		self.nl = nl
		self.paramchar = paramchar
		self.path = path

		self.norm = np.loadtxt( self.path+"NEr_"+ self.paramchar+'.d' )[:,1]
		self.error = np.loadtxt( self.path+"NEr_"+ self.paramchar+'.d' )[:,2]
		photon_data = np.loadtxt( self.path+"PHOTONS_"+ self.paramchar+'.d' )
		self.pplt = np.loadtxt( self.path+"PPLT_"+ self.paramchar+'.d' )[:,1:]
		if exists (self.path+"ENTROPY_"+ self.paramchar+'.d'):
			self.entropy = np.loadtxt( self.path+"ENTROPY_"+ self.paramchar+'.d' )[:,1]

		self.nk = photon_data[:,1:]
		self.nmodes = len( photon_data[0,1:])

		if exists (self.path+"APADAG_"+ self.paramchar+'.d'):
			apadag_data = np.loadtxt( self.path+"APADAG_"+ self.paramchar+'.d' )
			self.apadag = apadag_data[:,1]
		else:
			print('apadag does not exist')

		file = open( self.path+"LOG_"+ self.paramchar+'.d', "r")
		self.log = file.read()
		file.close()
		#self.omat = np.reshape( np.loadtxt( self.path+"OMAT_"+ self.paramchar+'.d' ), (self.nmodes, self.nmodes) )

		#-- resconstructing the final states from the weigths and the displacements
		ps_data_re_im = np.loadtxt( self.path+"PS_"+self.paramchar+'.d' )
		yks_data_re_im = np.loadtxt( self.path+"YKS_"+self.paramchar+'.d' )
		#yks_db_data_re_im = np.loadtxt( self.path+"YKS_DB_"+self.paramchar+'.d' )

		#-- deducing the number the final number of coherent states
		self.ncs_max = ( len( ps_data_re_im[:] )//2 ) // self.nl  

		ps_data_complex = ps_data_re_im[ :self.ncs_max*self.nl ] + 1j* ps_data_re_im[ self.ncs_max*self.nl: ]
		yks_data_complex = yks_data_re_im[:, 1:self.ncs_max*self.nl+1 ] \
				+ 1j* yks_data_re_im[ : , self.ncs_max*self.nl+1: ]
		#yks_db_data_complex = yks_db_data_re_im[:, 1:self.ncs_max*self.nl+1 ] \
		#		+ 1j* yks_db_data_re_im[ : , self.ncs_max*self.nl+1: ]

		ps = np.zeros( (self.nl,self.ncs_max), dtype='complex128' )
		yks = np.zeros( (self.nl, self.ncs_max, self.nmodes), dtype='complex128' )
		#yks_db = np.zeros( (self.nl, self.ncs_max, self.nmodes), dtype='complex128' )
		for i in range( self.nl ):
			for n in range( self.ncs_max ):
				ps[i,n] = ps_data_complex[ i*self.ncs_max + n ]
				for k in range( self.nmodes ):
					yks[i,n,k] = yks_data_complex[ k, i*self.ncs_max + n ]
					#yks_db[i,n,k] = yks_db_data_complex[ k, i*self.ncs_max + n ]

		self.times = photon_data[:,0]
		self.wc = yks_data_re_im[0,0]
		self.wk = yks_data_re_im[1:,0]
		self.dw = self.wk[1] - self.wk[0]

		self.final_state = state( ps, yks )
		#self.final_state_db = state( ps, yks_db )

		if exists ( self.path+"ovm_eval_"+ self.paramchar+'.d' ):
			ovm_evals_data = np.loadtxt( self.path+"ovm_eval_"+ self.paramchar+'.d' )[:,1:]
			self.ovm_evals_0 = ovm_evals_data[:,:self.ncs_max]
			self.ovm_evals_1 = ovm_evals_data[:,self.ncs_max:]
			for i in range( len(self.ovm_evals_0[:,0]) ):
				for j in range( len(self.ovm_evals_0[i,:]) ):
					if (self.ovm_evals_0[i,j] > 1e5):
						#self.ovm_evals_0[i,j] = -self.ovm_evals_0[i,j]
						self.ovm_evals_0[i,j] = 1
					if (self.ovm_evals_1[i,j] > 1e5):
						#self.ovm_evals_1[i,j] = -self.ovm_evals_1[i,j]
						self.ovm_evals_1[i,j] = 1
				
				indices = np.argsort( np.abs(1-self.ovm_evals_0[i,:]) )
				self.ovm_evals_0[i,:] = np.flip( self.ovm_evals_0[i,indices] )
				indices = np.argsort( np.abs(1-self.ovm_evals_1[i,:]) )
				self.ovm_evals_1[i,:] = np.flip( self.ovm_evals_1[i,indices] )

		if exists ( self.path+"PS2_"+ self.paramchar+'.d' ):
			ps2_data = np.loadtxt( self.path+"PS2_"+ self.paramchar+'.d' )[:,1:]
			self.ps2_0 = ps2_data[:,:self.ncs_max]
			self.ps2_1 = ps2_data[:,self.ncs_max:]
			for i in range( len(self.ps2_0[:,0]) ):
				indices = np.argsort( self.ps2_0[i,:] )
				self.ps2_0[i,:] = np.flip( self.ps2_0[i,indices] )
				indices = np.argsort( self.ps2_1[i,:] )
				self.ps2_1[i,:] = np.flip( self.ps2_1[i,indices] )


	def plot_wigner(self, xmin, log_min=None, add_cs_centers=False, path=None, label=None ):

		from matplotlib.colors import LogNorm
		from matplotlib.colors import Normalize

		xmax=-xmin
		wigner = fortran_module.calc_wigner( self.final_state.p, \
				self.final_state.y[:,:,0], self.final_state.ovm, xmin, xmax, 100 )

		plt.title('full wigner')
		plt.axhline( y=0, dashes=[2,2,2,2] )
		plt.axvline( x=0, dashes=[2,2,2,2] ) 

		if add_cs_centers:
			for s in range(self.final_state.nl):
				for n in range(self.final_state.ncs):
					size = 1 + np.abs(self.final_state.p[s,n])*50
					plt.scatter( np.real(self.final_state.y[s,n,0]), np.imag(self.final_state.y[s,n,0]), s=size, color='C'+str(2+s) )

		if log_min:
			norm_husimi=LogNorm(vmin=log_min)
			plt.imshow( wigner, extent=[xmin,xmax,xmin,xmax],norm=norm_husimi,origin ='lower' )
		else:
			plt.imshow( wigner, extent=[xmin,xmax,xmin,xmax],origin ='lower' )
		plt.colorbar()
		plt.tight_layout()

		if path:
			plt.savefig( path+label,format='pdf' )

		plt.show()

	def plot_split_wigner(self, s, xmin, log_min=None, add_cs_centers=False ):

		from matplotlib.colors import LogNorm
		from matplotlib.colors import Normalize

		xmax=-xmin
		wigner_ = fortran_module.calc_split_wigner( s+1, self.final_state.p[:,:], \
				self.final_state.y[:,:,0], self.final_state.ovm, xmin, xmax, 100 )

		plt.axhline( y=0, dashes=[2,2,2,2] )
		plt.axvline( x=0, dashes=[2,2,2,2] ) 

		for n in range(self.final_state.ncs):
			size = 1 + np.abs(self.final_state.p[s,n])*50
			plt.scatter( np.real(self.final_state.y[s,n,0]), np.imag(self.final_state.y[s,n,0]), s=size, color='C'+str(2+s) )

		if log_min:
			norm_husimi=LogNorm(vmin=log_min,vmax =wigner_.max()  )
			plt.imshow( wigner_, extent=[xmin,xmax,xmin,xmax],norm=norm_husimi,origin ='lower' )
		else:
			#print('here')
			plt.imshow( wigner_, extent=[xmin,xmax,xmin,xmax],origin ='lower', vmin=wigner_.min(), vmax=wigner_.max() )
		plt.colorbar()
		plt.show()

	def plot_split_wigners(self,  xmin, log_min=None, add_cs_centers=False ):

		from matplotlib.colors import LogNorm
		from matplotlib.colors import Normalize

		xmax=-xmin
		figure, axis = plt.subplots(1,4,figsize=(16,9),gridspec_kw={'width_ratios': [1,0.05, 1,0.05]}) 

		for s in range(2):
			axis[s*2].set_title('state=|'+str(s)+'>')
			wigner_ = fortran_module.calc_split_wigner( s+1, self.final_state.p[:,:], \
				self.final_state.y[:,:,0], self.final_state.ovm, xmin, xmax, 100 )
			axis[s*2].axhline( y=0, dashes=[2,2,2,2] )
			axis[s*2].axvline( x=0, dashes=[2,2,2,2] ) 
			for n in range(self.final_state.ncs):
				size = 1 + np.abs(self.final_state.p[s,n])*50
				axis[s*2].scatter( np.real(self.final_state.y[s,n,0]), np.imag(self.final_state.y[s,n,0]), s=size, color='C'+str(2+s) )

			if log_min:
				norm_husimi=LogNorm(vmin=log_min,vmax =wigner_.max()  )
				im=axis[s*2].imshow( wigner_, extent=[xmin,xmax,xmin,xmax],norm=norm_husimi,origin ='lower' )
			else:
				#print('here')
				im=axis[s*2].imshow( wigner_, extent=[xmin,xmax,xmin,xmax],origin ='lower', vmin=wigner_.min(), vmax=wigner_.max() )
			figure.colorbar(im, cax=axis[2*s+1], orientation='vertical')

		plt.show()


