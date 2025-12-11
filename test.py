import math

g_e = 9.8066
frontal_area_1 = 35.4612
frontal_area_2_3 = 20.2683
r_gas = 287
r_e = 6.356766E6
y = 0

# GET EFFECTIVE GRAV
def get_gravity(y):
    g = g_e*(r_e / (r_e+ y))**2
    return g

def get_temperature(y):
    a_1 = -6.5E-3
    a_2 = 3E-3
    a_3 = -4.5E-3
    a_4 = 4E-3
    h = (r_e / (r_e + y)) * y
    
    if 0 <= h and h <= 11000:
        t = 288.16 + (a_1*h)
    elif 11000 <= h and h <= 25000:
        t = 216.66
    elif 25000 <= h and h <= 47000:
        t = 216.66 + (a_2*(h-25000))
    elif 47000 <= h and h <= 53000:
        t = 282.66
    elif 53000 <= h and h <= 79000:
        t = 282.66 + (a_3*(h-53000))
    elif 79000 <= h and h <= 90000:
        t = 165.55
    else:
        t = 165.55 + (a_4*(h-90000))
        
    return t

a_1 = -6.5E-3
a_2 = 3E-3
a_3 = -4.5E-3
a_4 = 4E-3
    
y = 11000
p_base = 101325
p = p_base*((get_temperature(y)/288.16)**(-9.80621/(a_1*r_gas)))

print(p)