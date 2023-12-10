import math

#the values which must be given as inputs
pi = math.pi

#function to calculate thickness of the diffusion layer (kappa inverse)
#nI = the electrolyte number concentration
#k_b = boltzman constant
#temp = temperature in Kelvin 
#e = charge of an electron
#epsilon = dielectric constant
def kappa(nI, k_b, temp, e, epsilon):
    return ((2*nI*(e**2))/(epsilon*k_b*temp))**0.5


#function to calculate CM factor k_w
#epsilonStarP = complex permitivity of particle
#epsilonStarM = conplex permitivty of medium
def k_w(epsilonStarM, epsilonStarP):
    return (epsilonStarP - epsilonStarM)/(epsilonStarP + 2 * epsilonStarM)


#function to calculate time average Dielectrophoresis forces (F_DEP)
#Real Part of the Clausius Mossotti (CM) factor (Re_k_w)
#E_rms = rms value of electric field
#K(w) = CM factor which determines the extent of the polarization
def F_DEP(pi, epsilonM, a, Re_k_w, E_rms):
    return 2*pi*epsilonM*(a**3)*Re_k_w*(E_rms**2)


#function to calculate epsilonStar
def epsilonStar(epsilon,sigma, j, w):
    return epsilon-j*(sigma/w)


#epsilon p effective
def epsilonStarPEff(r2, r1, epsilonStarL, k_w):
    f1 = (r2/r1)**2
    f2 = k_w
    return epsilonStarL*((f1 + 2*f2)/(f1 - f2))


#r2 = radius of particle + EDL
r1 = 0  #radius of the particle
delta = 0 #thickness of the double layers
r2 = r1 + delta


#function to calculate the charge density (sigma star)
#sigmaStar = charge density
#phiNote = surface potential
#epsilon = permitivity of the free space / dielectric constant
#K = CM factor
def chargeDensity(r1, phiNote, epsilon, K):
    return ((phiNote * epsilon * (1 + K * r1)) / r1)
    

#function to calculate epsilonL
def epsilonL(sigmaStar,phiNote,kappa):
    return sigmaStar/(phiNote*kappa)


#conductivity of the particle
#sigmaPCore = conductivity of core particle
def sigmaP(sigmaPcore, K_s, r2):
    return sigmaPcore * ((2 * K_s) / r2)


#function to calculate the general Surface Conductance (Ks)
#sigmaStar = charge density
#miuI = ion mobility
def genSurfaceConductance(sigmaStar,miuI):
    return sigmaStar * miuI


