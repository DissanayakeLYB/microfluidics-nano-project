import math

#the values which must be given as inputs
pi = math.pi


#variables
ni = 5.09338*(10**(-3))
kappa = 3.995*(10**(-3)) #nm

epsilonM = 6.903*(10**(-10))

epsilonL_60 = 7.82*(10**(-10))         #for 60nm
epsilonL_30 = 8.73*(10**(-10))         #for 30nm

sigmaL_60 = 0.028         #for 60nm
sigmaL_30 = 0.05         #for 30nm

sigmaM = 5.5*(10**(-6))

zeta_60 = -37.6
zeta_30 = -33.6

miuI = 6.54*(10**(-8))

#r2 = radius of particle + EDL
r1 = 0  #radius of the particle
delta = 0 #thickness of the double layers
r2 = r1 + delta


#kappaInverse = function to calculate thickness of the diffusion layer
#nI = the electrolyte number concentration
#k_b = boltzman constant
#temp = temperature in Kelvin 
#e = charge of an electron
#epsilon = dielectric constant
def kappaInverse(nI, k_b, temp, e, epsilon):
    return ((epsilon*k_b*temp)/2*nI*(e**2))**0.5


#sigmaStar = charge density
#phiNote = surface potential
#epsilon = permitivity of the free space / dielectric constant = zeta potential
#K = CM factor (kappa)
def chargeDensity(r1, phiNote, epsilon, K):
    return ((phiNote * epsilon * (1 + K * r1)) / r1)


#epsilonL = permittivity of the layer
#sigmaStar = charge density
#phiNote = surface potential
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
def epsilonStar(epsilon,sigma, j, w):
    return epsilon-j*(sigma/w)


#epsilonStarPEff = effective permittivity for the equivalent particle 
#r1 = radius of the particle
#r2 = radius of particle + EDL
def epsilonStarPEff(r2, r1, epsilonStarL, k_w):
    f1 = (r2/r1)**2
    f2 = k_w
    return epsilonStarL*((f1 + 2*f2)/(f1 - f2))


#k_w = function to calculate CM factor 
#epsilonStarP = complex permitivity of particle
#epsilonStarM = conplex permitivty of medium
def k_w(epsilonStarM, epsilonStarP):
    return (epsilonStarP - epsilonStarM)/(epsilonStarP + 2 * epsilonStarM)


#timebyF_DEP = function to calculate time average by Dielectrophoresis forces 
#Re_k_w = Real Part of the Clausius Mossotti (CM) factor 
#E_rms = rms value of electric field
#epsilonM = permitivty of medium
#r2 = radius of particle + EDL 
def timeByF_DEP(pi, epsilonM, r2, Re_k_w, E_rms):
    return 2*pi*epsilonM*(r2**3)*Re_k_w*(E_rms**2)


#N_Re = Reynolds number 
#den = density
#ν = velocity
#dia = diameter of the particle
#vis = viscosity
def Rey_num(den, v, dia, vis):
    return (den*v*dia)/vis


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


#v_particle = velocity of the particle 
#r2 = radius of particle + EDL
#epsilon_m = permitivity of the medium
#Re_k_w = Real Part of the Clausius Mossotti (CM) factor 
#E_rms = rms value of electric field
#v_fluid = velocity of the fluid
def v_particle(r2,epsilon_m,Re_k_w,vis,E_rms,v_fluid):
    return (((r2**2)*epsilon_m*Re_k_w*(E_rms**2))/3*vis)+(v_fluid)













































 

"""
#graphical representation
import matplotlib.pyplot as plt
import numpy as np

x = [1,2,3,4,5,6,7,8,9,10]
y = [2,3,2,3,4,3,6,8,10,17]

mymodel = np.poly1d(np.polyfit(x, y, 3))

myline = np.linspace(1, 10, 100)

plt.scatter(x,y)
plt.xlabel("x coordinates")
plt.ylabel("y coordinates")
plt.title("x Vs. y")

plt.plot(myline, mymodel(myline))

plt.show()
"""


















"""#function to calculate thickness of EDL + particle radius
#U_p = particle velocity
#U = fluid velocity
#miu = dynamic velocity of the fluid
#epsilon_m = permitivity of the medium
#Real Part of the Clausius Mossotti (CM) factor (Re_k_w)
#E_rms = electric field
def doubleLayerThickness(U_p, U, miu, epsilonM, Re_k_w, E_rms):
    return (((U_p - U) * 3 * miu)/(epsilonM * Re_k_w * E_rms))**0.5"""
