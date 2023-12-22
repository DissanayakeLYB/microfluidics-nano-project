import math

#the values which must be given as inputs
pi = math.pi


#variables
nI = 5.09338*(10**(-3))
kappa = 3.995*(10**(-3)) #nm

epsilonM = 6.903*(10**(-10))

epsilonL_60 = 7.82*(10**(-10))         #for 60nm
epsilonL_30 = 8.73*(10**(-10))         #for 30nm

sigmaL_60 = 0.028         #for 60nm
sigmaL_30 = 0.05         #for 30nm

sigma_m = 5.5*(10**(-6))

zeta_60 = -37.6     #for 60nm for pH = 7.3
zeta_30 = -33.6     #for 30nm for pH = 7.06

miuI = 6.54*(10**(-8)) #miuI = ion mobility

k_b = 1.3806*(10**(-23)) #boltzman constant

e = 1.602*(10**(-19)) #charge of an electron

#r2 = radius of particle + EDL
r1 = 60*(10**(-9))  #radius of the core particle
delta = 3.995*(10**(-9)) #thickness of the double layers
r2 = r1 + delta




#-------------------------------------------------------------------------------




#kappaInverse = thickness of the diffusion layer (Debye Length)
#nI = electrolyte number concentration
#k_b = boltzman constant
#temp = temperature in Kelvin 
#e = charge of an electron
#epsilon_m = permittivity of the medium
def kappaInverse(nI, k_b, temp, e, epsilon_m):
    return ((epsilon_m*k_b*temp)/2*nI*(e**2))**0.5


#sigmaStar = charge density
#phiNote = surface potential (zeta potential)
#epsilon_m = permittivity of the medium
#kappaInverse = thickness of the diffusion layer
#r1 = radius of the core particle
def chargeDensity(r1, phiNote, epsilon_m, kappaInverse):
    return ((phiNote * epsilon_m * (1 + (r1/kappaInverse))) / r1)


#epsilonL = permittivity of the layer
#sigmaStar = charge density
#phiNote = surface potential (zeta potential)
#kappaInverse = function to calculate thickness of the diffusion layer
def epsilonL(sigmaStar,phiNote,kappaInverse):
    return (sigmaStar*kappaInverse)/phiNote


#Ks = the general Surface Conductance
#sigmaStar = charge density
#miuI = ion mobility
def genSurfaceConductance(sigmaStar,miuI):
    return sigmaStar * miuI


#sigmaP = conductivity of the particle 
#sigmaPCore = conductivity of core particle
#Ks = function to calculate the general Surface Conductance
#r2 = radius of particle + EDL
def sigmaP(sigmaPcore, K_s, r2):
    return sigmaPcore + ((2 * K_s) / r2)


#epsilonStar = complex permittivity 
#epsilon = permittivity
#sigma = conductivity
#j = imaginery unit
#w = angular frequency of the electric field
def epsilonStar(epsilon,sigma, w):
    return complex(epsilon,(sigma/w))

#k_w = CM factor 
#epsilonStarP = complex permitivity of particle
#epsilonStarM = conplex permitivty of medium
def k_w(epsilonStarM, epsilonStarP):
    return (epsilonStarP - epsilonStarM)/(epsilonStarP + 2 * epsilonStarM)


#epsilonStarPEff = effective permittivity for the equivalent particle 
#r1 = radius of the core particle
#r2 = radius of particle + EDL
def epsilonStarPEff(r2, r1, epsilonStarL, k_w):
    f1 = (r2/r1)**2
    f2 = k_w
    return epsilonStarL*((f1 + 2*f2)/(f1 - f2))


#slip velocity = velocity difference between fluid and a particle
#miu = mobility
#epsilon_m = permitivity of the medium
#Re_k_w = Real Part of the Clausius Mossotti (CM) factor 
#E_rms = rms value of electric field
def slip_velocity(r2, miu, epsilon_m, Re_k_w, E_rms):
    return (((r2**2) * epsilon_m * Re_k_w * (E_rms**2))/( 3 * miu))**0.5


#assuming F_DEP = -F_drag (neglcting Brownian motion and buoyancy) = dielectrophoretic force acting on them during movement
#vis = viscosity
#r2 = radius of particle + EDL 
#v = velocity
def F_DEP(vis,r2,v):
    return 6*pi*vis*r2*v


#miuDEP = DEP mobility into direction of the electric field 
#r2 = radius of particle + EDL
#epsilon_m = permitivity of the medium
#Re_k_w = Real Part of the Clausius Mossotti (CM) factor 
#vis = viscosity
def miuDEP(r2,epsilon_m,Re_k_w,vis):
    return ((r2**2)*epsilon_m*Re_k_w)/3*vis




#-------------------------------------------------------------------------------
 
#trials for 60nm gold particles
particle_size = 60
temp = 298
epsilon_m = 6.903*(10**(-10)) #relative permitivity
kappa_in = kappaInverse(nI, k_b, temp, e, epsilon_m)

phiNote = -37.6*(10**(-3)) #zeta potential
epsilon_m = 6.903*(10**(-10))
sigmaStar = chargeDensity(r1, phiNote, epsilon_m, kappa_in)

eps_L = epsilonL(sigmaStar, phiNote, kappa_in)

K_s = genSurfaceConductance(sigmaStar,miuI)

#here, estimated conductivity of the core particle to be same as EDL considering the thickness of stern layer is neglegible
sigmaPcore = sigmaL_60
conduct_particle = sigmaP(sigmaPcore, K_s, r2)

epsilonL = epsilonL_60 
f = 10**7
w = 2*pi*f #changeable from 10^2 to 10^10
epsilonStarL = epsilonStar(epsilonL_60,sigmaL_60, w)

epsilonStarP = epsilonStar(epsilonL_60,sigmaL_60, w)
epsilonStarM = epsilonStar(epsilon_m,sigma_m, w)
k_w = k_w(epsilonStarM, epsilonStarP)

#effective permittivity for the equivalent particle
effective_conductivity = epsilonStarPEff(r2, r1, epsilonStarL, k_w)

Re_k_w = k_w.real

vis = 8.9*(10**(-4)) #deionized water

mobility_DEP = miuDEP(r2, epsilon_m, Re_k_w, vis)



x = []
y = []
for f in range(1,10):
    freq = 10**f
    w = 2*pi*freq #changeable from 10^2 to 10^10
    epsilonStarL = epsilonStar(epsilonL_60, sigmaL_60, w)

    epsilonStarP = epsilonStar(epsilonL_60, sigmaL_60, w)
    epsilonStarM = epsilonStar(epsilon_m, sigma_m, w)
    k_w = k_w(epsilonStarM, epsilonStarP)

    #effective permittivity for the equivalent particle
    effective_conductivity = epsilonStarPEff(r2, r1, epsilonStarL, k_w)

    E = 195 #E = applied voltage (vary between a range of 5 to 200) 
    #E_rms = E/(2**0.5)
    Re_k_w = k_w.real
    E_rms = E/(2**0.5)
    v = slip_velocity(r2, vis, epsilon_m, Re_k_w, E_rms)
    DEP_force = F_DEP(vis, r2, v)
    x.append(f)
    y.append(DEP_force)


print(f"Particle Size = {particle_size} nm")
print(f"Temperature = {temp} K")
print(f"Zeta potential = {phiNote} V")
print(f"Frequency = {f} Hz")      
print(f"Viscosity of the medium = {vis} PaS")
print(f"Effective Permittivity = {effective_conductivity} S/m")
print(f"Applied voltage difference = {E} V")
print(f"Velocity difference between the particles and the fluid = {v} m/s")
print(f"Dielectrophoresis force applied on the particle = {DEP_force} N")

#-------------------------------------------------------------------------------



#graphical representation
import matplotlib.pyplot as plt
import numpy as np

mymodel = np.poly1d(np.polyfit(x, y, 3))

myline = np.linspace(40, 160, 100)

plt.scatter(x,y)
plt.xlabel("Frequency (f)")
plt.ylabel("DEP Force (N)")
plt.title("Applied Voltage Vs. DEP Force")

plt.plot(myline, mymodel(myline))

plt.show()






















"""
#timebyF_DEP = time average Dielectrophoresis forces 
#Re_k_w = Real Part of the Clausius Mossotti (CM) factor 
#E_rms = rms value of electric field
#epsilon_m = permittivity of the medium
#r2 = radius of particle + EDL 
def timeByF_DEP(pi, epsilon_m, r2, Re_k_w, E_rms):
    return 2*pi*epsilon_m*(r2**3)*Re_k_w*(E_rms**2)
"""


"""
#used to predict fluid flow patterns (laminar or turbulent)
#N_Re = Reynolds number 
#den = density
#ν = velocity
#dia = diameter of the particle
#vis = viscosity
def Rey_num(den, v, dia, vis):
    return (den*v*dia)/vis
"""


"""
#v_particle = velocity of the particle 
#r2 = radius of particle + EDL
#v_fluid = velocity of the fluid
def v_particle(slip_velocity,v_fluid):
    return slip_velocity + v_fluid
"""