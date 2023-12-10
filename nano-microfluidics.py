#the values which must be given as inputs
var_01 = 0
var_02 = 0
var_03 = 0
var_04 = 0
var_05 = 0
var_06 = 0

#epsilon p effective
def epsilonStarPEff(epsilonStarL, r2, r1, epsilonStarP):
    f1 = (r2/r1)**2
    f2 = (epsilonStarP - epsilonStarL)/(epsilonStarP + 2*epsilonStarL)


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
    

#conductivity of the particle
#sigmaPCore = conductivity of core particle
def sigmaP(sigmaPcore, K_s, r2):
    return sigmaPcore * ((2 * K_s) / r2)


#function to calculate the general Surface Conductance (Ks)
#sigmaStar = charge density
#miuI = ion mobility
def genSurfaceConductance(sigmaStar,miuI):
    return sigmaStar * miuI



sigma_p = sigmaP(var_01,var_02,var_03)
K_s = genSurfaceConductance(var_01,var_02,var_03)
