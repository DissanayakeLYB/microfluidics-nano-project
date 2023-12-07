#the values which must be given as inputs
var_01 = 0
var_02 = 0
var_03 = 0
var_04 = 0
var_05 = 0
var_06 = 0


#function to calculate the length of the microfludic channels
def lenChannel(a,b,c):
    channel_length = a + b + c
    return channel_length

#function to calculate the electric field required
def elecField(p,q,r):
    electric_field = p + q + r
    return electric_field

req_channel_length = lenChannel(var_01,var_02,var_03)
req_electric_field = elecField(var_01,var_02,var_03)

print(f"Required length of the microfluidic channel : {req_channel_length}")
print(f"Required electric field for the microfluidic channel : {req_electric_field}")