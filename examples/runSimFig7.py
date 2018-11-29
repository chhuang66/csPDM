#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
from math import log
import os


# set neuron parameters
C_val = 1.0             # membrane capacitance in uF/cm^2
gl_val = 0.05           # leaky conductance in mS/cm^2
Dt_val = 2.0            # sharpness factor in mV
El_val = -65.0          # leaky reversal potential in mV
Vt_val = -50.0          # threshold voltage in mV
Vc_val = -40.0          # cutting voltage in mV
Vr_val = -65.0          # resetting voltage in mV
taum_val = C_val/gl_val
a_val = 0.0             # subthreshold adaptation strength of the adaptive eponential integrate-and-fire model in mS
b_val=0.0               # spike adaptation current of the adaptive eponential integrate-and-fire model in uA
tauw_val = 600.0        # time constant for adaptation current of the adaptive eponential integrate-and-fire model in ms 
taur_val = 3.0          # refractory period in ms

time_step = 0.2         # the size of time step for numerical integration in ms

# set synaptic dynamics
# synaptic time constants
taue_val = 5.0          # for AMPA receptors
taui_val = 10.0         # for GABAa receptors
tausi_val = 100.0       # for GABAb receptors
Ee_val = 0.0            # reversal potential of AMPA in mV
Ei_val = -80.0          # reversal potential of GABAa in mV
Esi_val = -100.0        # reversal potential of GABAb in mV
# sizes of conductance jumps caused by pre-synaptic spikes 
vJ = 1.0                                    # instant voltage jump for excitatory input in mV
jGe_ws = -C_val*np.log(1.-vJ/(Ee_val-El_val))
vJ = 0.25                                   # instant voltage jump for inhibitory input in mV
jGi_ws = -C_val*np.log(1.+vJ/(Ei_val-El_val))
jGsi_ws = -C_val*np.log(1.+vJ/(Esi_val-El_val))

jGe_val = jGe_ws/taue_val                   # conductance jump size for AMPA
jGi_val = jGi_ws/taui_val                   # conductance jump size for GABAa
jGsi_val = jGsi_ws/tausi_val                # conductance jump size for GABAb

# set required parameters for local Galerkin method
N = 600                 # the number of meshes in V-direction
Vlb = -100.             # the lower bound of the V-direction in mV, the upper bound is Vc

# set connetion numbers
extAMPA = 1             # the number of synaptic connections through AMPA receptors
extGABAa = 1            # the number of synaptic connections through GABAa receptors
extGABAb = 1            # the number of synaptic connections through GABAb receptors

# initialize simulations
fname = "simulation.h5"     # file name
param = "--PYparam=a={},b={},tauw={},C={}.Gl={},Vr={},Vl={},Dt={},Vc={},Vt={},Vlb={},jGe={},jGi={},jGsi={},taue={},taui={},tausi={},Ve={},Vi={},Vsi={},t_ref={},N={}".format(a_val,b_val,tauw_val,C_val,gl_val,Vr_val,El_val,Dt_val,Vc_val,Vt_val,Vlb_val,jGe_val,jGi_val,jGsi_val,taue_val,taui_val,tausi_val,Ee_val,Ei_val,Esi_val,taur_val,N)
conn = "-cextAMP={},extGABAa={},extGABAb={}".format(extAMPA,extGABAa,extGABAb)

# run simulation of 3 seconds
os.system('aEIFONETHREE_semiLocal -R+3000 -ptime_step={} {} {} --ext_inFile=input.txt {}'.format(time_step,conn,param,fname))


