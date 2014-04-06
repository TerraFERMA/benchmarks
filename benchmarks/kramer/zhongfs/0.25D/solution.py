
from math import pi, cos, cosh, sinh, exp
from numpy import sqrt, interp

wavelengthfactor = 1.0
dx = 1./160. # nond

eta0 = 0. # m
xi0 = 0. # m
depth=750000.0 # m
rho0 = 4500. # kg/m^3
rhou = 2*rho0 # kg/m^3
mu0 = 1.e21  # Pas
g = 10.      # m/s^2
D = 3.e6     # m
deltaT = 2000 #K
alpha = 2.e-5 #/K


def wavenumber():
  wavelength = wavelengthfactor*D
  return 2.0*pi/wavelength

def nond_wavenumber():
  wavelength = wavelengthfactor
  return 2.0*pi/wavelength

def t0_eta():
  delta_rho = (rho0-rhou)*g
  rhog = rho0*g
  mu = mu0
  k = wavenumber()
  return -1.0*(((delta_rho*k**2*mu - k**2*mu*rhog)*D*sinh(D*k)**2 - (delta_rho*k**2*mu - k**2*mu*rhog)*D*cosh(D*k)**2 - (delta_rho*k*mu -\
k*mu*rhog)*sinh(D*k)*cosh(D*k) - sqrt((delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2)*D**2*cosh(D*k)**4 - 2*(delta_rho**2*k -\
2*delta_rho*k*rhog + k*rhog**2)*D*sinh(D*k)**3*cosh(D*k) + 2*(delta_rho**2*k - 2*delta_rho*k*rhog + k*rhog**2)*D*sinh(D*k)*cosh(D*k)**3 +\
((delta_rho**2*k**2 + 2*delta_rho*k**2*rhog + k**2*rhog**2)*D**2 + 4*delta_rho*rhog)*sinh(D*k)**4 - (2*(delta_rho**2*k**2 + k**2*rhog**2)*D**2 -\
delta_rho**2 + 2*delta_rho*rhog - rhog**2)*sinh(D*k)**2*cosh(D*k)**2)*k*mu)/(delta_rho*rhog*sinh(D*k)**2))

def t0_xi():
  delta_rho = (rho0-rhou)*g
  rhog = rho0*g
  mu = mu0
  k = wavenumber()
  return -1.0*(((delta_rho*k**2*mu - k**2*mu*rhog)*D*sinh(D*k)**2 - (delta_rho*k**2*mu - k**2*mu*rhog)*D*cosh(D*k)**2 - (delta_rho*k*mu -\
k*mu*rhog)*sinh(D*k)*cosh(D*k) + sqrt((delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2)*D**2*cosh(D*k)**4 - 2*(delta_rho**2*k -\
2*delta_rho*k*rhog + k*rhog**2)*D*sinh(D*k)**3*cosh(D*k) + 2*(delta_rho**2*k - 2*delta_rho*k*rhog + k*rhog**2)*D*sinh(D*k)*cosh(D*k)**3 +\
((delta_rho**2*k**2 + 2*delta_rho*k**2*rhog + k**2*rhog**2)*D**2 + 4*delta_rho*rhog)*sinh(D*k)**4 - (2*(delta_rho**2*k**2 + k**2*rhog**2)*D**2 -\
delta_rho**2 + 2*delta_rho*rhog - rhog**2)*sinh(D*k)**2*cosh(D*k)**2)*k*mu)/(delta_rho*rhog*sinh(D*k)**2))

def t0():
  return min(t0_eta(), t0_xi())

def nond_factor():
  return rho0*g*D*t0()/mu0

def nond_eta0():
  return eta0/D

def nond_xi0():
  return xi0/D

def nond_F(x,t):
  k = nond_wavenumber()
  F0 = nond_eta0()
  G0 = nond_xi0()
  delta_rho = (rho0-rhou)*nond_factor()/rho0
  zprime = depth/D
  T0 = deltaT
  alphag = nond_factor()*alpha*T0
  rhog = nond_factor() # use this as a proxy for the nondimensional factorisation
  return ((-0.5*(exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k - k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 - k**2*rhog -\
sqrt((delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k -\
2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) - (4*delta_rho*k**2*rhog - delta_rho**2 + 2*delta_rho*rhog -\
rhog**2)*sinh(k)**2)*k))/((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho +\
rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) - 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 -\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))) - exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k - k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 - k**2*rhog +\
sqrt((delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k -\
2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) - (4*delta_rho*k**2*rhog - delta_rho**2 + 2*delta_rho*rhog -\
rhog**2)*sinh(k)**2)*k))/((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho +\
rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) - 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 -\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))))*(G0 + (alphag*k*sinh(k*zprime)*cosh(k) - (alphag*k*zprime*cosh(k*zprime) -\
alphag*sinh(k*zprime))*sinh(k))/(T0*delta_rho*sinh(k)**2)) + (F0 - ((alphag*k*zprime*cosh(k*zprime) -\
alphag*sinh(k*zprime))*sinh(k)*cosh(k) - (alphag*k*zprime*sinh(k*zprime) - alphag*cosh(k*zprime))*sinh(k)**2 -\
alphag*k*sinh(k*zprime))/(T0*rhog*sinh(k)**2))*(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho +\
rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) - 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 -\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k +\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))) + 1)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k - k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 - k**2*rhog\
+ sqrt((delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k\
- 2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) - (4*delta_rho*k**2*rhog - delta_rho**2 + 2*delta_rho*rhog - rhog**2)*sinh(k)**2)*k)) -\
  (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k - k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 - k**2*rhog -\
sqrt((delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k -\
2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) - (4*delta_rho*k**2*rhog - delta_rho**2 + 2*delta_rho*rhog -\
rhog**2)*sinh(k)**2)*k))/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho +\
rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) - 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 -\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k +\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))) + ((alphag*k*zprime*cosh(k*zprime) - alphag*sinh(k*zprime))*sinh(k)*cosh(k) -\
(alphag*k*zprime*sinh(k*zprime) - alphag*cosh(k*zprime))*sinh(k)**2 - alphag*k*sinh(k*zprime))/(T0*rhog*sinh(k)**2)))*cos(k*x)

def nond_G(x,t):
  k = nond_wavenumber()
  F0 = nond_eta0()
  G0 = nond_xi0()
  delta_rho = (rho0-rhou)*nond_factor()/rho0
  zprime = depth/D
  T0 = deltaT
  alphag = nond_factor()*alpha*T0
  rhog = nond_factor() # use this as a proxy for the nondimensional factorisation
  return ((((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k - k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 - k**2*rhog -\
sqrt((delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k -\
2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) - (4*delta_rho*k**2*rhog - delta_rho**2 + 2*delta_rho*rhog -\
rhog**2)*sinh(k)**2)*k))/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho +\
rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) - 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 -\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k +\
k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k - k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 - k**2*rhog +\
sqrt((delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k -\
2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) - (4*delta_rho*k**2*rhog - delta_rho**2 + 2*delta_rho*rhog -\
rhog**2)*sinh(k)**2)*k))/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho +\
rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) - 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 -\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k +\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))))*(G0 + (alphag*k*sinh(k*zprime)*cosh(k) - (alphag*k*zprime*cosh(k*zprime) -\
alphag*sinh(k*zprime))*sinh(k))/(T0*delta_rho*sinh(k)**2)) - 2*(F0 - ((alphag*k*zprime*cosh(k*zprime) -\
alphag*sinh(k*zprime))*sinh(k)*cosh(k) - (alphag*k*zprime*sinh(k*zprime) - alphag*cosh(k*zprime))*sinh(k)**2 -\
alphag*k*sinh(k*zprime))/(T0*rhog*sinh(k)**2))*(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho +\
rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 -\
2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4\
+ k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) +\
2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) - 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 -\
2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 -\
2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3\
- rhog*sinh(k)*cosh(k)**2)/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2\
  + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k +\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))) + 1)*(k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k - k*rhog)*sinh(k)*cosh(k) + delta_rho*k**2 - k**2*rhog +\
sqrt((delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 - 2*delta_rho*k**2*rhog + k**2*rhog**2 + 2*(delta_rho**2*k -\
2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) - (4*delta_rho*k**2*rhog - delta_rho**2 + 2*delta_rho*rhog -\
rhog**2)*sinh(k)**2)*k))/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2 +\
sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4\
- 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 -\
  2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) -\
4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2\
+ 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) -\
k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)**2*exp(-delta_rho*rhog*t*sinh(k)**2/((delta_rho*k - k*rhog)*sinh(k)*cosh(k)\
+ delta_rho*k**2 - k**2*rhog - sqrt((delta_rho**2 + 2*delta_rho*rhog + rhog**2)*sinh(k)**4 + delta_rho**2*k**2 - 2*delta_rho*k**2*rhog +\
k**2*rhog**2 + 2*(delta_rho**2*k - 2*delta_rho*k*rhog + k*rhog**2)*sinh(k)*cosh(k) - (4*delta_rho*k**2*rhog - delta_rho**2 +\
2*delta_rho*rhog - rhog**2)*sinh(k)**2)*k))/(((k*rhog*sinh(k)**2*cosh(k) - k*rhog*cosh(k)**3 + rhog*sinh(k)**3 -\
rhog*sinh(k)*cosh(k)**2)/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k + k*rhog)*cosh(k)**2 -\
sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4\
- 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 -\
  2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k) -\
4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2\
+ 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 + rhog**2*sinh(k)**2*cosh(k)**2)) - (k*rhog*sinh(k)**2*cosh(k) -\
k*rhog*cosh(k)**3 + rhog*sinh(k)**3 - rhog*sinh(k)*cosh(k)**2)/((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 +\
(delta_rho*k + k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 +\
delta_rho**2*k**2*cosh(k)**4 + 2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 -\
2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 + k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 +\
4*delta_rho*k*rhog*sinh(k)**3*cosh(k) - 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) +\
2*k*rhog**2*sinh(k)*cosh(k)**3 + delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))*((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k +\
k*rhog)*cosh(k)**2 - sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2))*((delta_rho + rhog)*sinh(k)*cosh(k) - (delta_rho*k + k*rhog)*sinh(k)**2 + (delta_rho*k +\
k*rhog)*cosh(k)**2 + sqrt(delta_rho**2*k**2*sinh(k)**4 - 2*delta_rho**2*k**2*sinh(k)**2*cosh(k)**2 + delta_rho**2*k**2*cosh(k)**4 +\
2*delta_rho*k**2*rhog*sinh(k)**4 - 2*delta_rho*k**2*rhog*cosh(k)**4 + k**2*rhog**2*sinh(k)**4 - 2*k**2*rhog**2*sinh(k)**2*cosh(k)**2 +\
k**2*rhog**2*cosh(k)**4 - 2*delta_rho**2*k*sinh(k)**3*cosh(k) + 2*delta_rho**2*k*sinh(k)*cosh(k)**3 + 4*delta_rho*k*rhog*sinh(k)**3*cosh(k)\
- 4*delta_rho*k*rhog*sinh(k)*cosh(k)**3 - 2*k*rhog**2*sinh(k)**3*cosh(k) + 2*k*rhog**2*sinh(k)*cosh(k)**3 +\
  delta_rho**2*sinh(k)**2*cosh(k)**2 + 4*delta_rho*rhog*sinh(k)**4 - 2*delta_rho*rhog*sinh(k)**2*cosh(k)**2 +\
rhog**2*sinh(k)**2*cosh(k)**2)))) - (alphag*k*sinh(k*zprime)*cosh(k) - (alphag*k*zprime*cosh(k*zprime) -\
alphag*sinh(k*zprime))*sinh(k))/(T0*delta_rho*sinh(k)**2)))*cos(k*x)


def nond_F_amp(t):
  return nond_F(0.0,t)

def nond_G_amp(t):
  return nond_G(0.0,t)

def numerical_F_amp(detn, t):
  return interp([t], detn["ElapsedTime"]["value"], detn["Stokes"]["ScaledFreeSurface"]["TopLeft"][0]/nond_factor())[0]

def numerical_G_amp(detn, t):
  # minus sign because \Delta\rho is negative
  return interp([t], detn["ElapsedTime"]["value"], -detn["Stokes"]["ScaledFreeSurface"]["BottomLeft"][0]/nond_factor())[0]

def nond_F_error_amp(detn, t):
  return abs(nond_F_amp(t)-numerical_F_amp(detn, t))

def nond_G_error_amp(detn, t):
  return abs(nond_G_amp(t)-numerical_G_amp(detn, t))

def nond_F_ss(x):
  k = nond_wavenumber()
  F0 = nond_eta0()
  G0 = nond_xi0()
  delta_rho = (rho0-rhou)*nond_factor()/rho0
  zprime = depth/D
  T0 = deltaT
  alphag = nond_factor()
  rhog = nond_factor() # use this as a proxy for the nondimensional factorisation
  return (((alphag*k*zprime*cosh(k*zprime) - alphag*sinh(k*zprime))*sinh(k)*cosh(k) - (alphag*k*zprime*sinh(k*zprime) - \
alphag*cosh(k*zprime))*sinh(k)**2 - alphag*k*sinh(k*zprime))/(T0*rhog*sinh(k)**2))*cos(k*x)

 
def nond_G_ss(x):
  k = nond_wavenumber()
  F0 = nond_eta0()
  G0 = nond_xi0()
  delta_rho = (rho0-rhou)*nond_factor()/rho0
  zprime = depth/D
  T0 = deltaT
  alphag = nond_factor()
  rhog = nond_factor() # use this as a proxy for the nondimensional factorisation
  return (-(alphag*k*sinh(k*zprime)*cosh(k) - (alphag*k*zprime*cosh(k*zprime) - alphag*sinh(k*zprime))*sinh(k))/\
(T0*delta_rho*sinh(k)**2))*cos(k*x)
 

