import matplotlib.pyplot as plt
import numpy as np
import math

drag_coeff = 0.1
r_earth = 6.356766E6
frontal_area_1 = 35.4612
frontal_area_2_3 = 20.2683

g_o = 9.80621 
r_gas = 287
a_values = [-6.5E-3, 3E-3, -4.5E-3, 4E-3]

altitude_floors = [0, 11000, 25000, 47000, 53000, 79000, 90000]
temperature_floors = [288.16, 216.66, 216.66, 282.66, 282.66, 165.55, 165.55]
a_values_floors = [-6.5E-3, -6.5E-3, 3E-3, 3E-3, -4.5E-3, -4.5E-3, 4E-3]

m_i = 854000
m_s1, m_s2, m_s3 = 88000, 8000, 4000
m_p1, m_p2, m_p3 = 546000, 155000, 35100
m_payload = 17400
i_1, i_2, i_3 = 271.6, 302, 316
bt_1, bt_2, bt_3 = 121.5, 190.0, 223.0

y = 0

##### GET EFFECTIVE GRAVITY #####

def get_gravity(h):
    g = g_o * ((6.371E6) / ((6.371E6)+h))**2
    return g

##### GET TEMPERATURE #####

def get_temperature(h):
    if 0 <= h and h <= 11000:
        t = 288.16 + (a_values[0]*h)
    elif 11000 <= h and h <= 25000:
        t = 216.66
    elif 25000 <= h and h <= 47000:
        t = 216.66 + (a_values[1]*(h-25000))
    elif 47000 <= h and h <= 53000:
        t = 282.66
    elif 53000 <= h and h <= 79000:
        t = 282.66 + (a_values[2]*(h-53000))
    elif 79000 <= h and h <= 90000:
        t = 165.55
    else:
        t = 165.55 + (a_values[3]*(h-90000)) 
    return t

##### GET PRESSURE BASES #####

def get_pressure_bases():
    p_bases = []
    p_base = 101325 # Pa
    for i in (range(len(altitude_floors) - 1)):
        h_base = altitude_floors[i]
        h_top  = altitude_floors[i+1]
        T_base = temperature_floors[i]
        T_top  = get_temperature(h_top)
        p_bases.append(p_base)
        if T_top != T_base: 
            a = a_values_floors[i]   
            p_new = p_base * (T_top/T_base)**(-g_o/(a*r_gas))
        else:                 
            T_iso = T_base   
            p_new = p_base * math.exp(-(g_o/(r_gas*T_iso))*(h_top-h_base))
        p_base = p_new
    p_bases.append(p_base)
    return p_bases
    
##### GET PRESSURES #####

def get_pressure(h):
    pressure_floors = get_pressure_bases()
    if 0 <= h and h <= 11000:
        p = pressure_floors[0] * ( (get_temperature(h)/temperature_floors[0])**(-g_o/(a_values[0]*r_gas)) ) # gradient
    elif 11000 <= h and h <= 25000:
        p = pressure_floors[1] * math.exp(-(g_o/(r_gas*get_temperature(h)))*(h-altitude_floors[1])) # isothermal
    elif 25000 <= h and h <= 47000:
        p = pressure_floors[2] * ( (get_temperature(h)/temperature_floors[2])**(-g_o/(a_values[1]*r_gas)) ) # gradient
    elif 47000 <= h and h <= 53000:
        p = pressure_floors[3] * math.exp(-(g_o/(r_gas*get_temperature(h)))*(h-altitude_floors[3])) # isothermal
    elif 53000 <= h and h <= 79000:
        p = pressure_floors[4] * ( (get_temperature(h)/temperature_floors[4])**(-g_o/(a_values[2]*r_gas)) ) # gradient
    elif 79000 <= h and h <= 90000:
        p = pressure_floors[5] * math.exp(-(g_o/(r_gas*get_temperature(h)))*(h-altitude_floors[5])) # isothermal
    else:
        p = pressure_floors[6] * ( (get_temperature(h)/temperature_floors[6])**(-g_o/(a_values[3]*r_gas)) ) # gradient
    return p

##### GET DENSITY BASES #####

def get_density_bases():
    rho_bases = []
    rho_base = 1.225 # Pa
    for i in (range(len(altitude_floors) - 1)):
        h_base = altitude_floors[i]
        h_top  = altitude_floors[i+1]
        T_base = temperature_floors[i]
        T_top  = get_temperature(h_top)
        rho_bases.append(rho_base)
        if T_top != T_base: 
            a = a_values_floors[i]   
            rho_new = rho_base * (T_top/T_base)**(-((g_o/(a*r_gas)) + 1))
        else:                 
            T_iso = T_base   
            rho_new = rho_base * math.exp(-(g_o/(r_gas*T_iso))*(h_top-h_base))
        rho_base = rho_new
    rho_bases.append(rho_base)
    return rho_bases
    
##### GET DENSITIES #####

def get_density(h):
    density_bases = get_density_bases()
    if 0 <= h and h <= 11000:
        rho = density_bases[0] * ( (get_temperature(h)/temperature_floors[0])**(-((g_o/(a_values[0]*r_gas)) + 1))) # gradient
    elif 11000 <= h and h <= 25000:
        rho = density_bases[1] * math.exp(-(g_o/(r_gas*get_temperature(h)))*(h-altitude_floors[1])) # isothermal
    elif 25000 <= h and h <= 47000:
        rho = density_bases[2] * ( (get_temperature(h)/temperature_floors[2])**(-((g_o/(a_values[1]*r_gas)) + 1))) # gradient
    elif 47000 <= h and h <= 53000:
        rho = density_bases[3] * math.exp(-(g_o/(r_gas*get_temperature(h)))*(h-altitude_floors[3])) # isothermal
    elif 53000 <= h and h <= 79000:
        rho = density_bases[4] * ( (get_temperature(h)/temperature_floors[4])**(-((g_o/(a_values[2]*r_gas)) + 1))) # gradient
    elif 79000 <= h and h <= 90000:
        rho = density_bases[5] * math.exp(-(g_o/(r_gas*get_temperature(h)))*(h-altitude_floors[5])) # isothermal
    else:
        rho = density_bases[6] * ( (get_temperature(h)/temperature_floors[6])**(-((g_o/(a_values[3]*r_gas)) + 1))) # gradient
    return rho

##### (OPTIONAL) GRAPH ENVIRONMENT CONDITIONS #####

def graph_environment_conditions():
    altitudes = np.linspace(0, 90000, 1000)
    temperatures = []
    pressures = []
    densities = []
    for i in range(len(altitudes)):
        temperatures.append(get_temperature(altitudes[i]))
        pressures.append(get_pressure(altitudes[i]))
        densities.append(get_density(altitudes[i]))
    
    plt.plot(densities, altitudes)
    plt.ylabel("Altitude (m)")
    plt.show()

##### GET DRAG FORCE #####

def get_drag_force(stage, u_magnitude, y):
    match stage:
        case 1:
            return 0.5 * drag_coeff * get_density(y) * frontal_area_1 * (u_magnitude**2)
        case _:
            return 0.5 * drag_coeff * get_density(y) * frontal_area_2_3 * (u_magnitude**2)
    
##### TIME STEP SETUP ##### 

def time_step(dt, mass_load, mass_structure, mass_propellant):
    time_step = dt
    total_mass = mass_load + mass_structure + mass_propellant
    
    # DETERMINE STAGE
    if mass_propellant > (m_p2+ m_p3):
        stage = 1
    elif mass_propellant > (m_p2):
        stage = 2
    elif mass_propellant > (0):
        stage = 3
    else:
        stage = 4
    
    
    stage = 1

    mass_propellant = m_p1
    i_sp = i_1
    g = g_o

    velocity = 0

    # CHANGE IN MASS
    fuel_flow = m_p1 / bt_1 
    mass_propellant = mass_propellant - fuel_flow 
    total_mass = mass_load + mass_structure + mass_propellant

    # FIND MASS RATIOS
    lambd_a = mass_load / (mass_structure + mass_propellant)
    epsilon = mass_structure / (mass_structure + mass_propellant)
    r = (1 + lambd_a) / (epsilon + lambd_a)

    theta = 0 # degrees
    u_eq = i_sp * g
    u_e = u_eq * math.log(r)
    u_e_x = u_e * math.sin(math.radians(theta))
    u_e_y = u_e * math.cos(math.radians(theta))

    u_d = time_step * (get_drag_force(1, velocity, y))
    u_d_x = - u_d * math.sin(math.radians(theta))
    u_d_y = - u_d * math.cos(math.radians(theta))

    u_g = u_g_y = time_step * g * math.cos(math.radians(theta))
    u_g_x = 0

##### LOOP #####

mass_structure = m_s1 + m_s2 + m_s3
mass_propellant = m_p1 + m_p2 + m_p3

while (y < 185000):
    dt = 0.01
    time_step(dt, m_payload, mass_structure, mass_propellant)