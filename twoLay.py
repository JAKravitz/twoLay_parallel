#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# SCRIPT USED TO GENERATE BACTERIA/DETRITAL IOPs: (esp phase fn)
# This uses a single junge distribution with slope (-4)

# This is currently the default code for generating file
# "DETRITAL_FLAT_FULL_RES_Qs.mat"

# had flat real RI at 1.04 (Stramski 2001, AO vol 40, no 18 p. 2935)
# size distribs etc. based on Babin et al 2003, Limnol. Oceanogr. 48(2)
# our code can't handle more than 120 micron though, but his fig. 4a
# indicates the relative contributions of sizes to bp and not much beyond
# 120. Small particles nb though. we go downto 0.05 micron, he goes down to
# 0.02 which may explain why our resulting bp(555) (see fig. 10, for junge
# slope 4) is bit lower than his.

@author: jakravit
"""
# from numpy import f2py
# sourcefile = open('/Users/jakravit/pyProjects/EAP/Dmmex_R14B_4.f','rb')
# sourcecode = sourcefile.read()
# f2py.compile(sourcecode, modulename='Dmmex_R14B_4')
import Dmmex_R14B_4
import scipy.io as io
import numpy as np
import pickle
from scipy.interpolate import griddata
from scipy.signal import savgol_filter

def twolay (l, kcore, nreal, j, dmax, rho):

    #l = np.arange(.3, .855, .005) # wavelength range and resolution (changing this changes your interp value when normalising kshell)
    # int_val = 55 # refers to l[55] which is 675 nm. 255 for 1 nm resolution
    
    Vs = 0.2
    D_eff = np.arange(2, 4, 2)# in detrital this is only used for the size of the array to put IOPs in ...
    
    V_eff= 0.6
    Vc=1 - Vs
    FR=(1- Vs) ** (1/ 3)# relative volume of core to shell
    nmedia = 1.334
    wavelength = l/ nmedia
    
    wvno = 2 * np.pi / wavelength # this is a Mie param - combo of size and wavelength
    kshell = kcore;
    
    # hilbert transform
    def analytic_signal(x):
        from scipy.fftpack import fft, ifft
        N = len(x)
        X = fft(x, N)
        h = np.zeros(N)
        h[0] = 1
        h[1:N//2] = 2* np.ones(N// 2-1)
        h[N// 2] = 1
        Z = X * h
        z = ifft(Z, N)
        return z
    
    ncore = nreal + np.imag(analytic_signal(kcore))
    nshell = ncore
    khom = kcore*Vc + kshell*Vs # real refractive index 
    nhom = ncore*Vc + nshell*Vs
    mshell = nshell - kshell*1j
    mcore = ncore - kcore*1j
    mhom = nhom - khom*1j
    
    psd = np.arange(0.1,dmax,0.1)
    #psd = np.arange(0.05,120.05,0.05)
    deltad=1
    deltadm=1
    theta901 = np.arange(0, 90.1, 0.1) # length 901
    
    nang=901
    angles=np.cos(np.deg2rad(theta901)) # length 901
    theta1 = np.deg2rad(theta901) 
    theta2 = np.flipud(np.pi-np.deg2rad(theta901[0:900]))
    theta=np.append(theta1,theta2)
    
    back1=np.where(theta==(np.pi/2)) 
    back2=np.where(theta==(np.pi))
    d1=np.diff(theta)
    dtheta = np.append(d1, d1[0])
    
    # preparing variables to be filled
    VSF = np.zeros((len(D_eff),len(l), 1801))   #dimensions jjj(deff), nii(wavelength), 1801 angles    
    PF_check = np.zeros((len(D_eff),len(l)))
    d_alpha = []  
    PF = np.zeros((len(D_eff), len(wavelength), 1801))
    
    # declare all lists and arrays to fill in the jjj loop (refilled each iteration)
    Qc, Sigma_c, c, Qb, Sigma_b, b, Qa, Sigma_a, a, Qbb, Sigma_bb, bb, bbtilde = (np.zeros([len(D_eff),len(l)]) for i in range(13))
    a_sol, a_solm, Qastar2_dir, Qas21 = (np.zeros([len(D_eff),len(l)]) for i in range(4))
    
    # the loops go as follows:
    # nii - for each wavelength
        # jj - for each size particle in the psd, do the diMilay, get the bbprob for each particle
        # jjj - for each Deff requested, calc the resulting efficiency factors and phase fns
            # ai - for each of 1801 angles, for each Deff
        # then calc the phase fns, package effect and so on. For each Deff, for each wavelength.
    # end.
    for nii in np.arange(0,len(l)): # this is the wavelength loop
        print(nii)
    
        # declare lists to be filled on each iteration of the psd loop
        II, phaseMB, alpha, bbprob, bbprob1, Qbbro, checkMB, Qbro, Qcro, M1 = ([] for i in range (10))
        
    
        for jj in np.arange(0, len(psd)): # this is the psd loop
    
            [QEXT,QSCA,GQSC,QBS,m1,m2,s21,d21] = Dmmex_R14B_4.dmilay((psd[jj]*FR)/2,psd[jj]/2,wvno[nii],mshell[nii],mcore[nii],angles,901,901)
            # the diMilay calculates each particle's optical density, first step in getting the efficiency factors
            # on each iteration of jj, we get a different QEXT and QSCA out. So these must be stored in their own array
            Qcro.insert(jj,QEXT) 
            Qbro.insert(jj,QSCA)    
            
            m1_seta = [num[0] for num in m1]
            m1_setb = [num[1] for num in m1]
            M1 = np.append(m1_seta,m1_setb[900:0:-1])
            M2 = np.append(m2[0:901,0], m2[900:0:-1,1])
            myval = (M1+M2)/2 
            II.insert(jj, myval) 
            
            #Calculate raw phase function [M&B '86, eq. 18] NB this does not produce
    		#phase function that satisfies normalisation conditions
            alpha2=2*np.pi*(psd[jj]/2)/wavelength[nii] 
            alpha.insert(jj, alpha2) 
    
            phaseMB_jj = [II[jj] / (np.pi* Qbro[jj]* (alpha[jj]**2))]
            phaseMB.insert(jj,phaseMB_jj)
    
            checkMB_jj = [2* np.pi* np.sum(phaseMB_jj * np.sin(theta) * dtheta)]
            checkMB.insert(jj,checkMB_jj)
    
            #Calculate backscattering probability from normalised phase function
            section_jj = [item[900:1801] for item in phaseMB_jj]
            bbprob_jj = 2*np.pi* np.sum((section_jj *np.sin(theta[900:1801]) *dtheta[900:1801]))
            bbprob.insert(jj, bbprob_jj) 
            Qbbro_jj = QSCA * bbprob_jj 
            Qbbro.insert(jj,Qbbro_jj) 
        
    	# we are still in the nii loop here! just the jj loop has ended
    
        d_alpha_nii = alpha[1] - alpha[0]
        d_alpha.insert(nii,d_alpha_nii)
    
    	# jjj loop starts here
        for jjj in np.arange(0,len(D_eff)):
            
              #Calculate equivalent Junge PSD
                # psd = np.arange(0.05,120.05,0.05) - defined earlier.
            slope = j; 
            junge = np.power( psd,( - slope));  
            psdm1 = psd / 1e6; # because it was in micron 1:1:100
            psdm2 = junge; # does this need anything doing to it to get into metres?
            deltadm = 1;
            rho_a = rho; # ie. 0.4 e6 g is 400 kg/m^3. Babin et al 2003. will give b of roughly 1. highest relative scattering because least dense.
            distrib_vol = np.pi/ 6 * sum(psdm2 * psdm1 **3 * deltadm); 
            psdm2 = psdm2 * (1./ (distrib_vol * rho_a)); # this serves as per the ci normalisation, so your output IOPs
            # are mass specific!
            rho_vol = np.pi/ 6 * sum(psdm2 * psdm1 **3 * deltadm); # distribution volume after normalisation
       
           # mass_check=rho_vol_check * rho_a; # - must be close to 1 
        
           # #for standard normal distribution 
           # exponent = (-psd/ 2)/ ((D_eff[jjj]/ 2) * V_eff)
           # psd2 = 1.0e20 * np.power((psd/2),((1-3* V_eff)/V_eff)) * np.exp(exponent)
           # # I think the 1e20 here is essentially just random so you can see a nice spread? It gets normalised afterwards.
             
           # psdm1 = psd / 1e6; # because it was in micron 1:1:100
           # psdm2 = psd2 * 1e3; # really not sure why its multiplied by 1000 here.  
           # civol = np.pi/ 6 * sum(psdm2 * psdm1 **3 * deltadm) # the total volume of the distribution - before normalising
           # psdm2 = psdm2 * (1./ (civol * ci)) # normalising to the ci 
           
           # psdvol = np.pi/6 * sum(psdm2 * np.power(psdm1, 3) * deltadm) # the resulting actual volume of the distribution
    		
            
            # calculating the optical efficiencies (proportion of incident light absorbed, scattered on particle cross section) 
            Qc[jjj, nii] = sum(Qcro *psdm2 * np.power(psdm1,2) * deltadm)/ sum(psdm2 * np.power(psdm1,2) *deltadm)
            Sigma_c[jjj,nii] = np.pi/4 * Qc[jjj, nii] * sum(np.power(psdm1, 2) * deltadm)
            c[jjj,nii] = np.pi/4* Qc[jjj, nii]* sum(psdm2* np.power(psdm1,2)* deltadm)
            
            Qb[jjj, nii] = sum(Qbro * psdm2 * np.power(psdm1,2) * deltadm) /sum(psdm2* np.power(psdm1,2)* deltadm) 	            
            Sigma_b[jjj,nii] = np.pi/4 * Qb[jjj,nii]* sum(np.power(psdm1,2)* deltadm)
            b[jjj, nii] = np.pi/4* Qb[jjj, nii]* sum(psdm2* np.power(psdm1,2)* deltadm)
    
            Qbb[jjj, nii] = sum(Qbbro * psdm2 * np.power(psdm1,2) * deltadm) /sum(psdm2* np.power(psdm1,2)* deltadm)
            Sigma_bb[jjj, nii] = np.pi/4 * Qbb[jjj, nii] * sum(np.power(psdm1, 2) * deltadm)
            bb[jjj, nii] =  np.pi/4* Qbb[jjj, nii]* sum(psdm2 * np.power(psdm1, 2) * deltadm)
            
            Qa[jjj, nii] = Qc[jjj, nii] - Qb[jjj, nii]
            Sigma_a[jjj, nii] = np.pi/4 * Qa[jjj, nii]* sum(np.power(psdm1,2)* deltadm)
            a[jjj, nii] = c[jjj, nii] - b[jjj, nii]
            
            # VSF
            betabar, VSF_1 = ([] for i in range(2))
            checkbar = []
            b_check, bb_check = (np.zeros((len(D_eff),len(wavelength))) for i in range(2))
            
            bbtilde[jjj, nii] = bb[jjj, nii] / b[jjj, nii]
    			
    		# this little sub loop is INSIDE the jjj loop		
            for ai in range (0, nang * 2 - 1): # this should go 1801 times - doesn't include the last in range  
        		   # need a variable to get(II(:,ai)):
                varII = [item[ai] for item in II]
                betabar_ai = (1 / np.pi) * (sum(varII * psdm2 * d_alpha[nii]) / sum(Qbro * psdm2 * np.power(alpha, 2) * d_alpha[nii]))
                betabar.insert(ai, betabar_ai)
                VSF_1_ai = betabar[ai] * b[jjj, nii]  
                VSF_1.insert(ai, VSF_1_ai) # this gives VSF_1 of 1801 angles. For the current instance of nii (wavelength) and jjj (Deff)
    
            # checkbar is back outside of the sub loop
            checkbar = (2* np.pi * sum(betabar * np.sin(theta) * dtheta))
            PF_check[jjj,nii] = checkbar
            
            b_check[jjj,nii] = 2 * np.pi * sum((VSF_1) * np.sin(theta) * dtheta)
            
            PF[jjj,nii,:] = betabar
            VSF[jjj,nii,:] = VSF_1 # VSF_1s are put into matrix on each iteration of Deff, then wavelength.
            # We want to get out all wavelengths but backscatter only angles for each Deff:
                
            slice_obj = slice(900, 1801) 
            VSF_b = VSF[jjj, nii, slice_obj] # want to get the backward angles for this instance of Deff and all the wavelengths 
            bb_check[jjj,nii] = 2 * np.pi * sum((VSF_b) * np.sin(theta[900: 1801]) * dtheta[900: 1801])  
            
     	    # both the jjj loop and the nii loop end here.

            result = {'Qc': Qc,
                       'Sigma_c': Sigma_c,
                       'cstar': c,
                       'Qb': Qb,
                       'Sigma_b': Sigma_b,
                       'bstar': b,
                       'Qa': Qa,
                       'Sigma_a': Sigma_a,
                       'astar': a,
                       'Qbb': Qbb,
                       'Sigma_bb': Sigma_bb,
                       'bbstar': savgol_filter(bb, 11, 3),
                       'VSF': VSF,
                       'VSF_theta': theta,
                       'VSF_angles': np.rad2deg(theta)} 
    
    return result


