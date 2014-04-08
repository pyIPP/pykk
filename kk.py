import ctypes
import scipy
import numpy
import point

libkk = ctypes.cdll.LoadLibrary('/afs/ipp-garching.mpg.de/aug/ads/lib64/amd64_sles11/libkk.so')


class point:
	R = 0.0
	z = 0.0
	rho = 0.0
	def __init__( self , z=0.0 , R=0.0 , rho=0.0 ):
		self.z = z
		self.R = R
		self.rho = rho

	def __lt__( self , right ):
		if self.z < right.z:
			return True
		return False

	def __le__( self , right ):
		if self.z <= right.z:
			return True
		return False

	def __gt__( self , right ):
		if self.z > right.z:
			return True
		return False

	def __ge__( self , right ):
		if self.z >= right.z:
			return True
		return False



class kk:
	__status = False
	__shotnumber = ctypes.c_int(0)
	__edition = ctypes.c_int(0)
	
	def __init__( self , shotnumber=None,experiment='AUGD',diagnostic='FPP',edition=0):
		self.__status = False
		self.__shotnumber = ctypes.c_int(0)
		self.__edition = ctypes.c_int(0)
		if shotnumber != None:
			self.Open( shotnumber , experiment , diagnostic , edition )
	
	
	def __del__( self ):
		self.Close()	
	
	
	def Open( self , shotnumber , exper='AUGD' , diag='FPP' , edition=0 ):
		self.Close()
		if shotnumber > 0:
			self.__status = True
			self.__shotnumber = ctypes.c_int(shotnumber)
			self.__edition = ctypes.c_int(edition)
			self.__diag = ctypes.c_char_p(diag)
			vars(self)['experiment'] = exper
			vars(self)['diagnostic'] = diag
			self.__exper = ctypes.c_char_p(exper)
		return True
		
	
	
	def Close( self ):
		if self.__status:
			self.__status = False
			self.__shotnumber = ctypes.c_int(0)
			self.__edition = ctypes.c_int(0)
			del self.experiment
			del self.diagnostic
		
	
	
	def get_B( self , time , R , z , exper='AUGD' , diag='FPP' ):
		result = {}
		result['Bt'] = 0.0
		result['Bpol'] = 0.0
		result['Br'] = 0.0
		result['Bz'] = 0.0
		result['time'] = 0.0
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			rin = ctypes.c_float(R)
			_rin = ctypes.byref(rin)
			zin = ctypes.c_float(z)
			_zin = ctypes.byref(zin)
			lin = ctypes.c_int(1)
			_lin = ctypes.byref(lin)
			br = ctypes.c_float(0)
			_br = ctypes.byref(br)
			bz = ctypes.c_float(0)
			_bz = ctypes.byref(bz)
			bt = ctypes.c_float(0)
			_bt = ctypes.byref(bt)
			fpf = ctypes.c_float(0)
			_fpf = ctypes.byref(fpf)
			jpol = ctypes.c_float(0)
			_jpol = ctypes.byref(jpol)
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)
			libkk.kkrzbrzt_( _error , self.__exper, self.__diag , _shotnumber , _edition , _t , _rin , _zin , _lin , _br , _bz , _bt , _fpf , _jpol , lexper , ldiag )
			result['Bt'] = bt.value
			result['Bpol'] = scipy.sqrt( br.value*br.value + bz.value*bz.value )
			result['Br'] = br.value
			result['Bz'] = bz.value
			result['time'] = t.value
		return result
	
	
	def Rz_to_rhopol( self , time , R , z ):
		N = numpy.size(R)
		result = (ctypes.c_float*N)()
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			if N == 1:
				rin = ctypes.c_float(R)
				zin = ctypes.c_float(z)
			else:
				rin = (ctypes.c_float*N)()
				zin = (ctypes.c_float*N)()
				for i in range(N):
					rin[i] = R[i]
					zin[i] = z[i]
			_rin = ctypes.byref(rin)
			_zin = ctypes.byref(zin)
			lin = ctypes.c_int(N)
			_lin = ctypes.byref(lin)
			fpf = (ctypes.c_float*N)()
			_fpf = ctypes.byref( fpf )
			_rhopol = ctypes.byref( result )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)
			libkk.kkrzpfn_(_error, self.__exper,self.__diag,_shotnumber,_edition,_t,_rin,_zin,_lin,_fpf,_rhopol,lexper,ldiag)
		return numpy.frombuffer(result, numpy.float32)
	
	
	def s_to_Rz( self ,  s ):
		N = numpy.size(s)
		if self.__status:
			error = ctypes.c_int( 0 )
			_error = ctypes.byref( error )
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			if N == 1:
				sin = ctypes.c_float( s )
			else:
				sin = (ctypes.c_float*N)()
				for i in range( N ):
					sin[i] = s[i]
			_sin = ctypes.byref( sin )
			Raus = (ctypes.c_float*N)()
			zaus = (ctypes.c_float*N)()
			aaus = (ctypes.c_float*N)()
			_Raus = ctypes.byref( Raus )
			_zaus = ctypes.byref( zaus )
			_aaus = ctypes.byref( aaus )
			length = ctypes.c_int( N )
			_length = ctypes.byref( length )
			lexper = ctypes.c_long( len(self.experiment)  )
			ldiag = ctypes.c_long( len(self.diagnostic) )
			libkk.kkgcsrza_(_error ,self.__exper,self.__diag,_shotnumber,_edition,_length,_sin,_Raus,_zaus,_aaus,lexper,ldiag)
			output = []
			for i in range( N ):
				output.append( point( zaus[i] , Raus[i] , 0.0 ) )
			return output
		return numpy.nan
	
	
	def rhopol_to_Rz( self, time, rhopol, angle, degrees=False ): # angle in degrees...
		N = numpy.size(rhopol)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			
			if N == 1:
				rhopf = ctypes.c_float(rhopol) if isinstance(rhopol, float) else ctypes.c_float(rhopol[0]) 
				Rn = ctypes.c_float(0)
				zn = ctypes.c_float(0)
			else:
				rhopf = (ctypes.c_float*N)()
				Rn = (ctypes.c_float*N)()
				zn = (ctypes.c_float*N)()
				for i in range(N):
					rhopf[i] = rhopol[i]
					Rn[i] = 0
					zn[i] = 0
			_rhopf = ctypes.byref(rhopf)
			_Rn = ctypes.byref(Rn)
			_zn = ctypes.byref(zn)
			lrho = ctypes.c_int(N)
			_lrho = ctypes.byref(lrho)
			iorg = ctypes.c_int(0)
			_iorg = ctypes.byref(iorg)
			ang = ctypes.c_float(angle) if degrees else ctypes.c_float(angle/numpy.pi*180.)
			_angle = ctypes.byref(ang)
			
			lexper = ctypes.c_long(len(self.__exper))
			ldiag = ctypes.c_long(3)
			
			libkk.kkrhorz_(_error, self.__exper,self.__diag,_shotnumber,_edition,_t,\
							_rhopf, _lrho, _iorg, _angle,\
							_Rn, _zn,\
							lexper, ldiag)
			if N == 1:
				return {'error' : numpy.int_(error),\
				        'time'  : numpy.float32(t),\
				        'R'     : numpy.float32(Rn),\
				        'z'     : numpy.float32(zn)}
			else:
				return {'error' : numpy.int_(error),\
				        'time'  : numpy.float32(t),\
				        'R'     : numpy.array(Rn),\
				        'z'     : numpy.array(zn)}
	
	
	def rhopol_to_q(self, time, rhopol):
		N = numpy.size(rhopol)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			
			if N == 1:
				rhopf = ctypes.c_float(rhopol) if isinstance(rhopol, float) else ctypes.c_float(rhopol[0]) 
				q = ctypes.c_float(0)
				pf = ctypes.c_float(0)
			else:
				rhopf = (ctypes.c_float*N)()
				q = (ctypes.c_float*N)()
				pf = (ctypes.c_float*N)()
				for i in range(N):
					rhopf[i] = rhopol[i]
					q[i] = 0
					pf[i] = 0
			_rhopf = ctypes.byref(rhopf)
			_q = ctypes.byref(q)
			_pf = ctypes.byref(pf)
			lrho = ctypes.c_int(N)
			_lrho = ctypes.byref(lrho)
					
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)
			
			libkk.kkrhopfq_(_error, self.__exper,self.__diag,_shotnumber,_edition,_t,\
			_rhopf, _lrho,\
			_q,_pf,\
			lexper, ldiag)
					
			if N == 1:
				return {'error' : numpy.int_(error),\
				'time'  : numpy.float32(t),\
				'q'     : numpy.float32(q),\
				'pf'     : numpy.float32(pf)}
			else:
				return {'error' : numpy.int_(error),\
				'time'  : numpy.float32(t),\
				'q'     : numpy.frombuffer(q, numpy.float32),\
				'pf'     : numpy.frombuffer(pf, numpy.float32)}
		
	def rhopol_to_rhotor(self, time, rhopol):
		N = numpy.size(rhopol)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			
			if N == 1:
				rhopf = ctypes.c_float(rhopol) if isinstance(rhopol, float) else ctypes.c_float(rhopol[0]) 
				rhot = ctypes.c_float(0)
				pf = ctypes.c_float(0)
				fpf = ctypes.c_float(0)
				ftf = ctypes.c_float(0)
			else:
				rhopf = (ctypes.c_float*N)()
				rhot = (ctypes.c_float*N)()
				pf = (ctypes.c_float*N)()
				fpf = (ctypes.c_float*N)()
				ftf = (ctypes.c_float*N)()
				for i in range(N):
					rhopf[i] = rhopol[i]
					rhot[i] = 0
					pf[i] = 0
					fpf[i] = 0
					ftf[i] = 0
			_rhopf = ctypes.byref(rhopf)
			_rhot = ctypes.byref(rhot)
			_pf = ctypes.byref(pf)
			lrho = ctypes.c_int(N)
			_lrho = ctypes.byref(lrho)
			_fpf = ctypes.byref(fpf)
			_ftf = ctypes.byref(ftf)
					
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)
			
			libkk.kkrhopto_(_error, self.__exper,self.__diag,_shotnumber,_edition,_t,\
			_rhopf, _lrho,\
			_rhot,_fpf,_ftf,\
			lexper, ldiag)

			#kkrhopto (&error,exp,diag,&shot,&edition, &time,
			#rhopf, &lrho,
			#rhopf, fpf, ftf , strlen(exp), strlen(diag) );
					
			if N == 1:
				return {'error' : numpy.int_(error),\
				'time'  : numpy.float32(t),\
				'rhotor'     : numpy.float32(rhot),\
				'fpf'     : numpy.float32(fpf),\
				'ftf': numpy.float32(ftf)}
			else:
				return {'error' : numpy.int_(error),\
				'time'  : numpy.float32(t),\
				'rhotor'     : numpy.array(rhot, dtype=float),\
				'fpf'     : numpy.array(fpf, dtype=float),\
				'ftf': numpy.array(ftf, dtype=float)}	
	
	def theta_to_sfla(self, time, q, theta, degrees=False):
		N = numpy.size(theta)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)
			cq = ctypes.c_float(q)
			_q = ctypes.byref(cq)
			swrad = ctypes.c_int(0) if degrees else ctypes.c_int(1)
			_swrad = ctypes.byref(swrad)
			langle = ctypes.c_int(N)
			_langle = ctypes.byref(langle)
			_Rmag = ctypes.byref(ctypes.c_float(0))
			_zmag = ctypes.byref(ctypes.c_float(0))
			_tSHf = ctypes.byref(ctypes.c_float(0)) #  no idea what that is, docs aren't clear
			
			# todo: N*_angle, N*_Rn/_zn, N* 
			if N == 1:
				angle = ctypes.c_float(theta)
				Rn = ctypes.c_float(0)
				zn = ctypes.c_float(0)
				thetsn = ctypes.c_float(0)
				Brn = ctypes.c_float(0)
				Bzn = ctypes.c_float(0)
				Btn = ctypes.c_float(0)
			else:
				angle  = (ctypes.c_float*N)()
				Rn     = (ctypes.c_float*N)()
				zn     = (ctypes.c_float*N)()
				thetsn = (ctypes.c_float*N)()
				Brn    = (ctypes.c_float*N)()
				Bzn    = (ctypes.c_float*N)()
				Btn    = (ctypes.c_float*N)()
				for i in range(N):
					angle  [i] = theta[i]
					Rn     [i] = 0
					zn     [i] = 0
					thetsn [i] = 0
					Brn    [i] = 0
					Bzn    [i] = 0
					Btn    [i] = 0
			
			_angle  = ctypes.byref(angle )
			_Rn     = ctypes.byref(Rn    )
			_zn     = ctypes.byref(zn    )
			_thetsn = ctypes.byref(thetsn)
			_Brn    = ctypes.byref(Brn   )
			_Bzn    = ctypes.byref(Bzn   )
			_Btn    = ctypes.byref(Btn   )
			
			libkk.kkeqqfl_(_error, self.__exper,self.__diag,_shotnumber,_edition,_t,\
			_q, _langle, _angle, _swrad,\
			_Rmag,_zmag,_Rn,_zn,        _tSHf,\
			_thetsn, _Brn, _Bzn, _Btn,\
			lexper, ldiag)
			
			if numpy.int_(error) != 0:
				print 'kkeqqfl_ error ', numpy.int_(error)
				if numpy.int_(error) == 14:
					print 'libkk note: 2*numpy\'s pi is a little too big for libkk.so; lower it a little bit'
				return {'error': numpy.int_(error)}
			
			if N == 1:
				return {'error' : numpy.int_(error),\
				'time' : numpy.float32(t),\
				'sfla' : numpy.float32(thetsn)}
			else:
				return {'error' : numpy.int_(error),\
				'time'  : numpy.float32(t),\
				'sfla'  : numpy.array(thetsn)}
	
	
	def psi_to_rhopol(self, time, psi):
		N = numpy.size(psi)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)
			lrho = ctypes.c_int(N)
			_lrho = ctypes.byref(lrho)
			if N == 1:
				PFi = ctypes.c_float(psi)
				rhoPF = ctypes.c_float(0)
			else:
				PFi = (ctypes.c_float*N)()
				rhoPF = (ctypes.c_float*N)()
				for i in range(N):
					PFi[i] = psi[i]
					rhoPF[i] = 0

			_PFi = ctypes.byref(PFi)
			_rhoPF = ctypes.byref(rhoPF)	

			#kkPFrhoP (iERR ,expnam,dianam,nSHOT,nEDIT,tSHOT,
			# 			> PFi,Lrho,
			# 			< rhoPF)
			libkk.kkpfrhop_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
				_PFi, _lrho,
				_rhoPF,
				lexper, ldiag)
			if N == 1:
				return numpy.float32(rhoPF)
			else:
				return numpy.frombuffer(rhoPF, numpy.float32)
	
	def get_jpol(self, time, N, NSOL=0):
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)

			PFL = (ctypes.c_float*(N+NSOL+1))()
			Jpol = (ctypes.c_float*(N+NSOL+1))()
			Jpolp = (ctypes.c_float*(N+NSOL+1))()
			_PFL   = ctypes.byref(PFL  )
			_Jpol  = ctypes.byref(Jpol )
			_Jpolp = ctypes.byref(Jpolp)


			N = ctypes.c_int(N)
			NSOL = ctypes.c_int(NSOL)
			_LPFp = ctypes.byref(N)
			_LPFe = ctypes.byref(NSOL)

		# kkEQJpolsol (iERR ,expnam,dianam,nSHOT,nEDIT,tSHOT,
		#            < LPFp,LPFe,PFL,Jpol,Jpolp)
			if NSOL.value != 0:
				libkk.kkeqjpolsol_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
					_LPFp, _LPFe, _PFL, _Jpol, _Jpolp,
					lexper, ldiag)
			else:
				libkk.kkeqjpol_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
					_LPFp, _PFL, _Jpol, _Jpolp,
					lexper, ldiag)
				# kkeqjpol returns values from the separatrix to the magnetic axis,
				# but kkeqjpolsol from the axis via the separatrix to the scrape off layer.
			return {'pfl':numpy.frombuffer(PFL, numpy.float32)[:N.value+NSOL.value+1],
				    'Jpol':numpy.frombuffer(Jpol, numpy.float32)[:N.value+NSOL.value+1],
				    'Jpolp':numpy.frombuffer(Jpolp, numpy.float32)[:N.value+NSOL.value+1],
				    'N':N.value,
				    'NSOL':NSOL.value}
		pass	

	def get_p(self, time, psi):
		N = numpy.size(psi)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)

			PFL = (ctypes.c_float*(N+1))()
			pres = (ctypes.c_float*(N+1))()
			presp = (ctypes.c_float*(N+1))()
			_PFL   = ctypes.byref(PFL  )
			_pres  = ctypes.byref(pres )
			_presp = ctypes.byref(presp)
			N = ctypes.c_int(N)
			_LPF = ctypes.byref(N)
			#kkeqpres (&error, exp,diag,&shot,&edition, &time,
			#          &lpf, pfl, pres, presp, strlen(exp), strlen(diag) );)
			libkk.kkeqpres_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
							_LPF, _PFL, _pres, _presp,
							lexper, ldiag)
			return {'pfl':numpy.frombuffer(PFL, numpy.float32)[:N.value+1],
			        'pres':numpy.frombuffer(pres, numpy.float32)[:N.value+1],
				    'presp':numpy.frombuffer(presp, numpy.float32)[:N.value+1],
				    'N':N.value}

	def get_special_points(self, time):
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)
			#kkeqpfx (&error, ,exp,diag,&shot,&edition, &time,
			#         &lpfx, pfxx, rpfx, zpfx,
			#         strlen(exp), strlen(diag) );
			N = 5
			_lpfx = ctypes.byref(ctypes.c_int(N-1))
			pfxx = (ctypes.c_float*N)()
			rpfx = (ctypes.c_float*N)()
			zpfx = (ctypes.c_float*N)()
			_pfxx = ctypes.byref(pfxx)
			_rpfx = ctypes.byref(rpfx)
			_zpfx = ctypes.byref(zpfx)

			#print self.__shotnumber.value

			libkk.kkeqpfx_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
						  _lpfx, _pfxx, _rpfx, _zpfx,
						  lexper, ldiag)

			#0. Magnetic axis
			#1. Primary X point
			#2. Primary limiter point
			#3. Secondary X point
			#4. Secondary limiter point

			return {'pfxx':numpy.frombuffer(pfxx, numpy.float32),
					'rpfx':numpy.frombuffer(rpfx, numpy.float32),
					'zpfx':numpy.frombuffer(zpfx, numpy.float32)}


	def get_ffprime(self, time, psi, typ=11):
		# typ:
		# flag if quantities should be read from shotfile or calculated:
		# 11: read from shotfile if available, otherwise calculate
		#  1: read from shotfile, return error if not available
		#  2: calculate, regardless if available in the shotfile


		N = numpy.size(psi)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)

			PFL = (ctypes.c_float*(N+1))()
			ffp = (ctypes.c_float*(N+1))()
			_PFL   = ctypes.byref(PFL  )
			_ffp  = ctypes.byref(ffp )
			N = ctypes.c_int(N)
			_LPF = ctypes.byref(N)

			_typ = ctypes.byref(ctypes.c_int(typ))		

			#kkeqffp_(&error,exp,diag,&shot,&edition, &time,
			#         &typ, &lpf, pfl, ffp, strlen(exp), strlen(diag) );
			libkk.kkeqffp_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
							_typ, _LPF, _PFL, _ffp,
							lexper, ldiag)
			return {'pfl':numpy.frombuffer(PFL, numpy.float32)[:N.value+1],
				'ffp':numpy.frombuffer(ffp, numpy.float32)[:N.value+1], 'N':N.value}

	def get_Rinv(self, time, psi, typ=11):
		# typ:
		# flag if quantities should be read from shotfile or calculated:
		# 11: read from shotfile if available, otherwise calculate
		#  1: read from shotfile, return error if not available
		#  2: calculate, regardless if available in the shotfile


		N = numpy.size(psi)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)

			PFL = (ctypes.c_float*(N+1))()
			rinv = (ctypes.c_float*(N+1))()
			r2inv = (ctypes.c_float*(N+1))()
			_PFL   = ctypes.byref(PFL  )
			_rinv  = ctypes.byref(rinv )
			_r2inv  = ctypes.byref(r2inv )
			N = ctypes.c_int(N)
			_LPF = ctypes.byref(N)

			_typ = ctypes.byref(ctypes.c_int(typ))		

			#kkeqffp_(&error,exp,diag,&shot,&edition, &time,
			#         &typ, &lpf, pfl, ffp, strlen(exp), strlen(diag) );
			libkk.kkeqrinv_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
							_typ, _LPF, _PFL, _rinv, _r2inv,
							lexper, ldiag)
			return {'pfl':numpy.frombuffer(PFL, numpy.float32)[:N.value+1],
				'rinv':numpy.frombuffer(rinv, numpy.float32)[:N.value+1],
				'r2inv':numpy.frombuffer(r2inv, numpy.float32)[:N.value+1], 'N':N.value}


	def get_jpar(self, time, N, typ=11):
			# void kkeqjpar_( int*, char*, char*, int*, int*, float*, 
			#                int*, int*, 
			#                float*, float*, float*, 
			#                long, long);
			if self.__status:
				error = ctypes.c_int(0)
				_error = ctypes.byref(error)
				_shotnumber = ctypes.byref(self.__shotnumber)
				_edition = ctypes.byref(self.__edition)
				t = ctypes.c_float( time )
				_t = ctypes.byref( t )
				lexper = ctypes.c_long(len(self.experiment))
				ldiag = ctypes.c_long(3)

				PFL = (ctypes.c_float*(N+1))()
				Jpar = (ctypes.c_float*(N+1))()
				Jparp = (ctypes.c_float*(N+1))()
				_PFL   = ctypes.byref(PFL  )
				_Jpar  = ctypes.byref(Jpar )
				_Jparp = ctypes.byref(Jparp)


				N = ctypes.c_int(N)
				_LPFp = ctypes.byref(N)
				_typ = ctypes.byref(ctypes.c_int(typ))		

				libkk.kkeqjpar_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
					_typ, _LPFp,
					_PFL, _Jpar, _Jparp,
					lexper, ldiag)

				return {'pfl':numpy.frombuffer(PFL, numpy.float32)[:N.value+1],
						    'Jpar':numpy.frombuffer(Jpar, numpy.float32)[:N.value+1],
						    'Jparp':numpy.frombuffer(Jparp, numpy.float32)[:N.value+1],
						    'N':N.value}

	def get_pfm(self, time, m=65, n=129):
		o_m = m
		o_n = n
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)

			mdim = ctypes.c_int(m)
			m = ctypes.c_int(m)
			n = ctypes.c_int(n)
			Ri = (ctypes.c_float*(m.value))()
			zj = (ctypes.c_float*(n.value))()
			pfm = ((ctypes.c_float*(m.value))*n.value)()

			_mdim = ctypes.byref(mdim)
			_m = ctypes.byref(m)
			_n = ctypes.byref(n)
			_Ri = ctypes.byref(Ri)
			_zj = ctypes.byref(zj)
			_pfm = ctypes.byref(pfm)

			# kkeqpfm_( int*, char*, char*, int*, int*, float*, 
            #   &mdim, &m, &n, Ri, zj, pfm,
            #   long, long);
			libkk.kkeqpfm_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
							_mdim, _m, _n, _Ri, _zj, _pfm,
							lexper, ldiag)

			return {
				'pfm':numpy.frombuffer(pfm, numpy.float32).reshape(o_n, o_m)[:n.value, :m.value],
				'Ri':numpy.frombuffer(Ri, numpy.float32)[:m.value],
				'zj':numpy.frombuffer(zj, numpy.float32)[:n.value]
				}
	
	def psi_to_v(self, time, psis):
		N = numpy.size(psis)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)

			fpf = (ctypes.c_float*N)()
			for i in range(N):
				fpf[i] = psis[i]
			_fpf = ctypes.byref(fpf)
			lin = ctypes.c_int(N)

			v = (ctypes.c_float*N)()
			work = (ctypes.c_float*N)()
			_v = ctypes.byref(v)
			_work = ctypes.byref(work)

			libkk.kkpfv(_error, self.__exper, self.__diag, self.__shotnumber,_edition,_t,
							fpf, lin,
							_v, _work)#,

			return {'volume': numpy.frombuffer(v, numpy.float32), 'workingarea': numpy.frombuffer(work, numpy.float32)}


			#libkk.kkpfv(byref(error), byref(exp), byref(diag), shot, byref(edition),
		    #        byref(tshot),
		    #        fpf, lin,
		    #        byref(v), byref(work) ) 

	def get_volume(self, time, psi):
		N = numpy.size(psi)
		if self.__status:
			error = ctypes.c_int(0)
			_error = ctypes.byref(error)
			_shotnumber = ctypes.byref(self.__shotnumber)
			_edition = ctypes.byref(self.__edition)
			t = ctypes.c_float( time )
			_t = ctypes.byref( t )
			lexper = ctypes.c_long(len(self.experiment))
			ldiag = ctypes.c_long(3)

			PFL = (ctypes.c_float*(N+1))()
			vol = (ctypes.c_float*(N+1))()
			volp = (ctypes.c_float*(N+1))()
			_PFL   = ctypes.byref(PFL  )
			_vol   = ctypes.byref(vol )
			_volp  = ctypes.byref(volp )
			N = ctypes.c_int(N)
			_LPF = ctypes.byref(N)

			libkk.kkeqvol(_error, self.__exper, self.__diag, self.__shotnumber,_edition,_t,
		                  _LPF, _PFL, _vol, _volp);

			return {'pfl':numpy.frombuffer(PFL, numpy.float32)[:N.value+1],
				'vol':numpy.frombuffer(vol, numpy.float32)[:N.value+1],
				'volp':numpy.frombuffer(volp, numpy.float32)[:N.value+1], 'N':N.value}		


		#libkk.kkeqvol(byref(error), exp,diag, shot, byref(edition), byref(time),
		#            byref(lpf), byref(pfl), byref(vol), byref(volp));

			#kkeqffp_(&error,exp,diag,&shot,&edition, &time,
			#         &typ, &lpf, pfl, ffp, strlen(exp), strlen(diag) );
			#libkk.kkeqffp_(_error, self.__exper, self.__diag, _shotnumber,_edition,_t,
			#				_typ, _LPF, _PFL, _ffp,
			#				lexper, ldiag)
			#return {'pfl':numpy.frombuffer(PFL, numpy.float32)[:N.value+1],
			#	'ffp':numpy.frombuffer(ffp, numpy.float32)[:N.value+1], 'N':N.value}



#void kkeqqfld_( int*, char*, char*, int*, int*, float*, 
#                float*, int*, float*, int*, 
#                float*, float*, float*, float*, float*, 
#                float*, float*, float*, float*, float*, 
#                long, long);

#C                        +---------------------------------+
#C              (3.4)      calc.  straight field line angle
#C                        +---------------------------------+
#C
#C     ___________________ {EQU,FPP} calc. straight f.l.angle for  qval:
#C                kkEQqFL  (iERR  ,expnam,dianam,nSHOT,nEDIT, tSHOT,
#C               >                 qval,Nangl,angle,swrad,
#C               <                 Rmag,zmag,Rn,zn,           tSHf ,
#C               <                 thetsn,     Brn,Bzn,Btn)
#C                kkEQqFLd (iERR  ,expnam,dianam,nSHOT,nEDIT, tSHOT,
#C               >                 qval,Nangl,angle,swrad,
#C               <                 Rmag,zmag,Rn,zn,           tSHf ,
#C               <                 thetsn,dtht,Brn,Bzn,Btn)
#C                -------------------------------------------------------
		

#C     (3.4)      calc.  straight field line angle ..kkEQqFL,d ( )
#C        qval   ...q_value                         (real* 4)      -> <--
#C        Nangl  ...# of elements in angle          (integer)      ->
#C        angle  ...1D array of angle [dgr] | [rad] (real *4)      ->
#C        swrad  ...switch for units of angles:                    ->
#C               := 0 ...angles in [degrees]
#C               := 1 ...angles in [radians]
#C        R,zmag ...position of magn.axis  [m]      (real *4)         <--
#C        Rn,zn  ...1D arrays of surface points     (real *4)         <--
#C        thetsn ...straight_field_line_angle       (real *4)         <--
#C               ...1D array: [dgr] | [rad]
#C        dtht   ...d/dtheta(thetsn)                (real *4)         <--
#C        Br,z,tn...1D arrays of B_components [T]   (real *4)         <--
#




