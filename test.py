  
##### TIME STEP SETUP ##### 
def time_step(x, y, dt, mass_load, mass_structure, mass_propellant, velocity, velocity_x, velocity_y , stage, isp):
    time_step = dt
    total_mass = mass_load + mass_structure + mass_propellant

    if stage == 1 and m_p1 <= 0:
        stage = 2
        isp = i_2
    elif (stage == 2 and m_p2 <= 0):
        stage = 3
        isp = i_3
    elif (stage == 3 and m_p3 <= 0): # BALLISTIC
        stage = 4  
        isp = 0.0
        
    g = get_gravity(y)
        
    # DETMERINE BERUN RATE
    if (stage == 1):
        current_stage_propellant = mass_propellant - m_p2 - m_p3
        fuel_flow = current_stage_propellant / bt_1 
    elif (stage == 2):
        current_stage_propellant = mass_propellant - m_p1 - m_p3
        fuel_flow = current_stage_propellant / bt_2
    elif (stage == 3):
        current_stage_propellant = mass_propellant - m_p1 - m_p2
        fuel_flow = current_stage_propellant / bt_3

    # CHANGE IN MASS
    mass_propellant = mass_propellant - (fuel_flow * dt)
    total_mass = mass_load + mass_structure + mass_propellant

    # FIND MASS RATIOS
    lambd_a = mass_load / (mass_structure + mass_propellant)
    epsilon = mass_structure / (mass_structure + mass_propellant)
    r = (1 + lambd_a) / (epsilon + lambd_a)

    theta = 0 # degrees
    u_eq = isp * g
    u_e = u_eq * math.log(r)
    u_e_x = u_e * math.sin(math.radians(theta))
    u_e_y = u_e * math.cos(math.radians(theta))
    
    

    u_d = time_step * ( (get_drag_force(stage, velocity, y)) / total_mass)
    u_d_x = - u_d * math.sin(math.radians(theta))
    u_d_y = - u_d * math.cos(math.radians(theta))

    u_g = u_g_y = - time_step * g * math.cos(math.radians(theta))
    u_g_x = 0
    
    u_x_total = u_e_x + u_d_x + u_g_x
    u_y_total = u_e_y + u_d_y + u_g_y
    velocity_x = velocity_x + u_x_total
    velocity_y = velocity_y + u_y_total
    
    velocity = math.sqrt((velocity_x**2) + (velocity_y**2))
    
    x = x + velocity_x * time_step
    y = y + velocity_y * time_step
    
    return x, y, dt, mass_load, mass_structure, mass_propellant, velocity, velocity_x, velocity_y, stage, isp

##### LOOP #####

mass_structure = m_s1 + m_s2 + m_s3
mass_propellant = m_p1 + m_p2 + m_p3
velocity = 0
velocity_x = 0
velocity_y = 0
x, y, = 0, 0
stage = 1
isp = i_1
x_postitions = []
y_postitions = []
t_time = []
t = 0

while (y < 185000):
    t_time.append(t)
    dt = 0.01 
    t = t + dt
    print(t)
    x_postitions.append(x)
    y_postitions.append(y)
    x, y, dt, mass_load, mass_structure, mass_propellant, velocity, velocity_x, velocity_y, stage, isp = time_step(x, y, dt, m_payload, mass_structure, mass_propellant, velocity, velocity_x, velocity_y, stage, isp)
    
plt.plot(x_postitions, t_time)
plt.xlabel(" Time (s) ")
plt.ylabel(" Altitude (m)")
plt.show()