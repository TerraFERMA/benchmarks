
from math import pi, cos, cosh, sinh, exp
from numpy import interp

wavelengthfactor = 1.0
dx = 1./80. # nond

eta0 = 1000. # m
rho0 = 4500. # kg/m^3
mu0 = 1.e21  # Pas
g = 10.      # m/s^2
D = 3.e6     # m

def k():
  wavelength = wavelengthfactor*D
  return 2.0*pi/wavelength

def nond_k():
  wavelength = wavelengthfactor
  return 2.0*pi/wavelength

def t0():
  return 2*(D*k() + sinh(D*k())*cosh(D*k()))*k()*mu0/(rho0*g*sinh(D*k())**2)

def nond_factor():
  return rho0*g*D*t0()/mu0

def nond_eta0():
  return eta0/D

def nond_F(x,t):
  F0 = nond_eta0()*cos(nond_k()*x)
  return exp(-t)*F0 # relaxation time is the time scale

def nond_F_amp(t):
  return exp(-t)*nond_eta0()

def numerical_F_amp(detn, t):
  return interp([t], detn["ElapsedTime"]["value"], detn["Stokes"]["ScaledFreeSurface"]["TopLeft"][0]/nond_factor())[0]

def nond_error_amp(detn, t):
  return abs(nond_F_amp(t)-numerical_F_amp(detn, t))


