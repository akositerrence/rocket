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
    h = (r_earth / (r_earth + y)) * y
    match stage:
        case 1:
            return 0.5 * drag_coeff * get_density(h) * frontal_area_1 * (u_magnitude**2)
        case _:
            return 0.5 * drag_coeff * get_density(h) * frontal_area_2_3 * (u_magnitude**2)
    
##### GET STAGES R #####

def get_stages_r(stage):
    # EFFECTIVE PAYLOADS
    ms_effective3 = m_payload
    ms_effective2 = m_s3 + m_p3 + m_payload
    ms_effective1 = m_s2 + m_p2 + m_s3 + m_p3 + m_payload
    ms_effective_list = [ms_effective1, ms_effective2, ms_effective3]
    m_s_list = [m_s1, m_s2, m_s3]
    m_p_list = [m_p1, m_p2, m_p3]
    r_list = []

    for i in range(3):
        m_eq = m_s_list[i] + m_p_list[i]
        lambd = ms_effective_list[i] / m_eq 
        epsilon = m_s_list[i] / m_eq 
        r = (1.0 + lambd) / (epsilon + lambd)
        r_list.append(r)
        
    current_r = r_list[stage]
    return current_r

##### STEP #####

def step(x, y, theta, velocity, velocity_x, velocity_y, time_step, stage, payload_mass, mass_propellant_1, mass_propellant_2, mass_propellant_3):
    if stage == 1 and mass_propellant_1 <= 0:
        stage = 2
    elif stage == 2 and mass_propellant_2 <= 0:
        stage = 3
    elif stage == 3 and mass_propellant_3 <= 0:
        stage = 4   

    match stage:
        case 1:
            current_isp   = i_1
            current_bt    = bt_1
            current_fuel  = mass_propellant_1
            current_total_structural_mass = 88000 + 8000 + 4000
            
            current_total_fuel_mass = current_fuel + mass_propellant_2 + mass_propellant_3
            current_total_mass_before_burn = current_total_structural_mass + current_total_fuel_mass + payload_mass
            
            fuel_flow = m_p1 / current_bt 
            if fuel_flow > current_fuel:
                fuel_flow = current_fuel
            current_fuel = current_fuel - (fuel_flow * time_step)
            
            current_total_fuel_mass = current_fuel + mass_propellant_2 + mass_propellant_3
            
        case 2:
            current_isp   = i_2
            current_bt    = bt_2
            current_fuel  = mass_propellant_2
            current_total_structural_mass = 8000 + 4000
            
            current_total_fuel_mass = current_fuel + mass_propellant_3
            current_total_mass_before_burn = current_total_structural_mass + current_total_fuel_mass + payload_mass
            
            fuel_flow = m_p2 / current_bt 
            if fuel_flow > current_fuel:
                fuel_flow = current_fuel
            current_fuel = current_fuel - (fuel_flow * time_step)
            
            current_total_fuel_mass = current_fuel + mass_propellant_3
            
        case 3:
            current_isp   = i_3
            current_bt    = bt_3
            current_fuel  = mass_propellant_3
            current_total_structural_mass = 4000
            
            current_total_fuel_mass = current_fuel
            current_total_mass_before_burn = current_total_structural_mass + current_total_fuel_mass + payload_mass
            
            fuel_flow = m_p3 / current_bt 
            if fuel_flow > current_fuel:
                fuel_flow = current_fuel
            current_fuel = current_fuel - (fuel_flow * time_step)
            
            current_total_fuel_mass = current_fuel 
            
        case _:  
            current_isp   = 0
            current_bt    = 1 
            current_fuel  = 0
            current_total_structural_mass = 0
            
            fuel_flow = 0
            current_fuel = 0
            
            current_total_fuel_mass = 0
    
    current_total_mass_after_burn = current_total_structural_mass + current_total_fuel_mass + payload_mass
    
    g = get_gravity(y)
    u_eq = current_isp * g
    u_e = u_eq * math.log(current_total_mass_before_burn / current_total_mass_after_burn)
    
    u_e_x = u_e * math.sin(math.radians(theta))
    u_e_y = u_e * math.cos(math.radians(theta))
    
    u_g_x = 0
    u_g_y = - time_step * g * math.cos(math.radians(theta))
    
    u_d = time_step * ( (get_drag_force(stage, velocity, y)) / current_total_mass_after_burn)
    u_d_x = - u_d * math.sin(math.radians(theta))
    u_d_y = - u_d * math.cos(math.radians(theta))
    
    u_x_total = u_e_x + u_g_x + u_d_x
    u_y_total = u_e_y + u_g_y + u_d_y
    
    velocity_x = velocity_x + u_x_total
    velocity_y = velocity_y + u_y_total
    velocity = math.sqrt((velocity_x**2) + (velocity_y**2))
    
    x = x + velocity_x * time_step
    y = y + velocity_y * time_step
    
    if stage == 1:
        mass_propellant_1 = current_fuel
    elif stage == 2:
        mass_propellant_2 = current_fuel
    elif stage == 3:
        mass_propellant_3 = current_fuel
    
    return x, y, theta, velocity, velocity_x, velocity_y, stage, mass_propellant_1, mass_propellant_2, mass_propellant_3

##### INITIALIZATION & LOOP #####

x = 0
y = 0
theta = 0

time_step = 0.01
stage = 1
payload_mass = m_payload
mass_propellant_1 = m_p1
mass_propellant_2 = m_p2
mass_propellant_3 = m_p3
velocity = 0
velocity_x = 0
velocity_y = 0
time = 0

x_postitions = []
y_postitions = []
times = []

while (y >= 0 and y <= 185000):
    time = time + time_step
    times.append(time)
    x_postitions.append(x)
    y_postitions.append(y)
    x, y, theta, velocity, velocity_x, velocity_y, stage, mass_propellant_1, mass_propellant_2, mass_propellant_3 = step(x, y, theta, velocity, velocity_x, velocity_y, time_step, stage, payload_mass, mass_propellant_1, mass_propellant_2, mass_propellant_3)
    
plt.plot(times, y_postitions)
print(max(y_postitions))
plt.xlabel(" Time (s) ")
plt.ylabel(" Altitude (m) ")
plt.show()