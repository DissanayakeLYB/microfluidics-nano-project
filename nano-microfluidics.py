#the values which must be given as inputs
var_01 = 0
var_02 = 0
var_03 = 0
var_04 = 0
var_05 = 0
var_06 = 0


#function to calculate the general Surface Conductance (Ks)
#sigmaStar = product of charge density
#miuI = ion mobility

def genSurfaceConductance(sigmaStar,miuI):
    return sigmaStar * miuI
    
#radii required of an particle
r1 = 0
delta = 0
r2 = r1 + delta

def sigmaP(sigmaPcore,K_s,r2):
    return sigmaPcore * (K_s/r2)

#function to calculate the general Surface Conductance (Ks)
#sigmaPCore = conductivity of core particle
def genSurfaceConductance(sigmaStar,miuI):
    return sigmaStar * miuI



sigma_p = sigmaP(var_01,var_02,var_03)
K_s = genSurfaceConductance(var_01,var_02,var_03)

print(f"Required length of the microfluidic channel : {req_channel_length}")
print(f"Required electric field for the microfluidic channel : {req_electric_field}")