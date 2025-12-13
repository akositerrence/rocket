import matplotlib.pyplot as plt
import numpy as np
import math
import csv

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

def step(time, x, y, theta, phi, velocity, velocity_x, velocity_y, 
         time_step, stage, payload_mass, mass_propellant_1, 
         mass_propellant_2, mass_propellant_3, gravity_turn_time,
         gravity_turn_theta_initial, gravity_turn_state, pd_state):
    
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
            dm = fuel_flow * time_step       
            if dm > current_fuel:
                dm = current_fuel
            current_fuel -= dm
            
            current_total_fuel_mass = current_fuel + mass_propellant_2 + mass_propellant_3
            
        case 2:
            current_isp   = i_2
            current_bt    = bt_2
            current_fuel  = mass_propellant_2
            current_total_structural_mass = 8000 + 4000
            
            current_total_fuel_mass = current_fuel + mass_propellant_3
            current_total_mass_before_burn = current_total_structural_mass + current_total_fuel_mass + payload_mass
            
            fuel_flow = m_p2 / current_bt 
            dm = fuel_flow * time_step       
            if dm > current_fuel:
                dm = current_fuel
            current_fuel -= dm
            
            current_total_fuel_mass = current_fuel + mass_propellant_3
            
        case 3:
            current_isp   = i_3
            current_bt    = bt_3
            current_fuel  = mass_propellant_3
            current_total_structural_mass = 4000
            
            current_total_fuel_mass = current_fuel
            current_total_mass_before_burn = current_total_structural_mass + current_total_fuel_mass + payload_mass
            
            fuel_flow = m_p3 / current_bt 
            dm = fuel_flow * time_step       
            if dm > current_fuel:
                dm = current_fuel
            current_fuel -= dm
            
            current_total_fuel_mass = current_fuel 
            
        case _:  
            current_isp   = 0
            current_bt    = 1 
            current_fuel  = 0
            current_total_structural_mass = 0

            current_total_fuel_mass = 0
            current_total_mass_before_burn = 17400
            
            fuel_flow = 0
            current_fuel = 0
            
            current_total_fuel_mass = 0
    
    ############################################################################

    if (gravity_turn_state == False) and (time >= gravity_turn_time):
        phi = math.radians(gravity_turn_theta_initial)
    else:
        phi = 0

    if (y >= 185000 - 2000):
        pd_state = True

    if pd_state == True:
        target_altitude = 185000
        velocity_y_target = 0
        phi_max = math.radians(1)
        
        if (current_isp > 0) and (y != target_altitude):
            altitude_error = y - target_altitude
            vy_error = velocity_y - velocity_y_target
            P = 1e-5
            D = 1e-5

            phi_cmd = (P * altitude_error) + (D * vy_error)

            if phi_cmd > phi_max:
                phi_cmd = phi_max
            elif phi_cmd < -phi_max:
                phi_cmd = -phi_max

            phi = phi_cmd

    ############################################################################
        
    current_total_mass_after_burn = current_total_structural_mass + current_total_fuel_mass + payload_mass
    
    g = get_gravity(y)
    u_eq = current_isp * g
    u_e = u_eq * math.log(current_total_mass_before_burn / current_total_mass_after_burn)
    
    u_e_x = u_e * math.sin(theta + phi)
    u_e_y = u_e * math.cos(theta + phi)
    
    u_g_x = 0
    u_g_y = - time_step * g
    
    u_d = time_step * ( (get_drag_force(stage, velocity, y)) / current_total_mass_after_burn)
    u_d_x = - u_d * math.sin(theta)
    u_d_y = - u_d * math.cos(theta)
    
    u_x_total = u_e_x + u_g_x + u_d_x
    u_y_total = u_e_y + u_g_y + u_d_y
    
    velocity_x = velocity_x + u_x_total
    velocity_y = velocity_y + u_y_total
    velocity = math.sqrt((velocity_x**2) + (velocity_y**2))
    theta = math.atan2(velocity_x, velocity_y)

    x = x + velocity_x * time_step
    y = y + velocity_y * time_step
    
    if stage == 1:
        mass_propellant_1 = current_fuel
    elif stage == 2:
        mass_propellant_2 = current_fuel
    elif stage == 3:
        mass_propellant_3 = current_fuel
    
    return time, x, y, theta, phi, velocity, velocity_x, velocity_y, stage, mass_propellant_1, mass_propellant_2, mass_propellant_3, gravity_turn_time, gravity_turn_theta_initial, gravity_turn_state, pd_state

##### INITIALIZATION & LOOP #####

x = 0
y = 0
theta = 0
phi = 0

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
velocities = []
thetas = []
times = []

stage_events = []
prev_stage = 1

gravity_turn_state = False
gravity_turn_time = 5
gravity_turn_theta_initial = 0.01

pd_state = False
pd_flag = False

with open("bruhlog.csv", "w", newline="") as f: 
    writer = csv.writer(f)
    writer.writerow(["bruh"])

    while (y >= 0 and y <= 10000000):
        time = time + time_step
        times.append(time)
        x_postitions.append(x)
        y_postitions.append(y)
        velocities.append(velocity)
        thetas.append(math.degrees(theta))
        time, x, y, theta, phi, velocity, velocity_x, velocity_y, stage, mass_propellant_1, mass_propellant_2, mass_propellant_3, gravity_turn_time, gravity_turn_theta_initial, gravity_turn_state, pd_state = step(time, x, y, theta, phi, velocity, velocity_x, velocity_y, time_step, stage, payload_mass, mass_propellant_1, mass_propellant_2, mass_propellant_3, gravity_turn_time, gravity_turn_theta_initial, gravity_turn_state, pd_state)
        if stage != prev_stage:
            stage_events.append((time, x, y, stage))
            prev_stage = stage
        if (pd_state == True) and (pd_flag == False):
            pd_flag = True
            stage_events.append((time, x, y, "PD"))
        # writer.writerow([velocity])

#for (t_evt, x_evt, y_evt, stg) in stage_events:
    #plt.scatter([x_evt], [y_evt], zorder=5)                     
    #plt.annotate(f"[{stg}]", (x_evt, y_evt),     
    #             textcoords="offset points", xytext=(8, 8))
    
# plt.text(2.5, 7.0, f'Apogee: {max(y_postitions)}', fontsize=11, color='black')



#### vvvv chatgpt wrote the stuff below this #####
# =========================
#   MULTI-PLOT SUMMARY
# =========================
t = np.array(times)
x_arr = np.array(x_postitions)
y_arr = np.array(y_postitions)
speed_arr = np.array(velocities)     # now speed magnitude because you appended `velocity`
theta_arr = np.array(thetas)

# Safety: align lengths (your arrays can be off by 1 depending on when you append)
n = min(len(t), len(speed_arr), len(x_arr), len(y_arr), len(theta_arr))
t = t[:n]; x_arr = x_arr[:n]; y_arr = y_arr[:n]; speed_arr = speed_arr[:n]; theta_arr = theta_arr[:n]

fig, ax = plt.subplots(2, 2, figsize=(14, 9))

# 1) Speed vs time
ax[0,0].plot(t, speed_arr)
ax[0,0].set_title("Speed vs Time")
ax[0,0].set_xlabel("Time (s)")
ax[0,0].set_ylabel("Speed (m/s)")
ax[0,0].grid(True)

# 2) Altitude vs time
ax[0,1].plot(t, y_arr)
ax[0,1].axhline(185000, linestyle="--", linewidth=1)
ax[0,1].set_title("Altitude vs Time")
ax[0,1].set_xlabel("Time (s)")
ax[0,1].set_ylabel("Altitude y (m)")
ax[0,1].grid(True)

# 3) Theta vs time
ax[1,0].plot(t, theta_arr)
ax[1,0].set_title("Flight Path Angle vs Time")
ax[1,0].set_xlabel("Time (s)")
ax[1,0].set_ylabel("theta (deg)")
ax[1,0].grid(True)
# --- mark PD start on theta vs time ---
pd_events = [(t_evt, x_evt, y_evt, stg) for (t_evt, x_evt, y_evt, stg) in stage_events if str(stg) == "PD"]
if len(pd_events) > 0:
    t_pd, x_pd, y_pd, stg = pd_events[0]

    # vertical line at PD start time
    ax[1,0].axvline(t_pd, linestyle="--", linewidth=1)

    # put a dot at the actual theta value at that time (nearest index)
    i_pd = int(np.argmin(np.abs(t - t_pd)))
    ax[1,0].scatter([t[i_pd]], [theta_arr[i_pd]], zorder=5)

    ax[1,0].annotate("PD start", (t[i_pd], theta_arr[i_pd]),
                     textcoords="offset points", xytext=(8, 8), fontsize=9)


# 4) Trajectory x vs y + stage/PD start markers
ax_xy = ax[1,1]
ax_xy.plot(x_arr, y_arr, linewidth=1.5)
ax_xy.axhline(185000, linestyle="--", linewidth=1)
ax_xy.set_title("Trajectory (x vs y)")
ax_xy.set_xlabel("x (m)")
ax_xy.set_ylabel("y (m)")
ax_xy.grid(True)

for (t_evt, x_evt, y_evt, stg) in stage_events:
    ax_xy.scatter([x_evt], [y_evt], s=40, zorder=5)
    ax_xy.annotate(f"{stg}", (x_evt, y_evt),
                   textcoords="offset points", xytext=(8, 8), fontsize=9)

plt.tight_layout()
plt.show()
