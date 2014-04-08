#!/usr/bin/env python
# encoding: utf-8
"""
test.py

Created by Alexander Bock on 2011-09-19.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

try:
	reload(kk)
except Exception, e:
	import pykk as kk


eq = kk.kk(); 


#eq.Open(29761, diag='EQI')
eq.Open(20993, 'PJM', 'EQE') #)'ABOCK', 'IDE')

if False:
	t = 3.0
	res = eq.get_pfm(t)
	X, Y = np.meshgrid(res['Ri'], res['zj'])
	plt.contour(X, Y, res['pfm'], 64, colors='k')
	plt.show()
	

if False:
	for t in [1., 2., 3.]:
		res = eq.get_jpar(t, 100)
		m = res['N']+1
		r = eq.psi_to_rhopol(t, res['pfl'])
		#print r


		plt.subplot(231)
		plt.title('j')
		plt.plot(r, res['Jpar'])
		plt.gca().invert_yaxis()
		plt.grid('on')
		plt.subplot(234)
		plt.title('j\'')
		plt.plot(r, res['Jparp'])
		plt.gca().invert_yaxis()
		plt.grid('on')


		p = eq.get_p(t, res['pfl'])
		
		plt.subplot(232)
		plt.title('p')
		plt.plot(r, p['pres'])
		plt.grid('on')
		plt.subplot(235)
		plt.title('p\'')
		plt.plot(r, p['presp'])
		plt.grid('on')


		ffp = eq.get_ffprime(t, res['pfl'])

		plt.subplot(236)
		plt.title('ff\'')
		plt.plot(r, ffp['ffp'])
		plt.grid('on')

	plt.show()

res = eq.get_volume(3.0, np.linspace(0,2, 100))

print res

#rhos = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
#rhos = [0.4, 0.4]
#asd = eq.rhopol_to_rhotor(t, rhos)
#print asd

#plt.plot(rhos, asd['rhotor'])
#plt.show()


#eq.Open(26827, diag='EQI'); 
#rhos = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
#asd = eq.rhopol_to_q(t, rhos);
#
#jkl = eq.theta_to_sfla(t, 2.0, np.linspace(0, 2*np.pi - 0.001))
#
#print jkl
#
#plt.plot(np.linspace(0, 2*np.pi), jkl['sfla'])
#plt.show()

#plt.plot(rhos, np.abs(asd["q"]))
#plt.show()