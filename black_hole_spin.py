# This code is being written as a recreation of the master's code I did in C.
# It is designed to simulate the growth rate and spin change of a supermassive black hole
# In my undergrad I had no concept of user defined functions and the code was monolithic and hard to read/follow
# The recreation is an attempt to tidy up my master's work, and learn python.

# The user defined functions are all storred in the file black_hole_functions.py in this git directory
import black_hole_functions as bh_func
import numpy as np
import matplotlib.pyplot as plt
import math

# set initial parameters
BH_mass = float(input("Input initial mass in solar masses:\n"))
alpha_1 = float(input("Input alpha 1 value:\n"))
alpha_2 = float(input("Input alpha 2 value:\n"))
Kerr_spin = float(input("Input initial Ker spin parameter (between -1 and 1):\n"))

#Convert initial mass from solar masses to kg
BH_mass_kg = BH_mass * 1.989e30
#Convert initial mass from solar masses to 10^8 solar mass units
BH_mass_8 = BH_mass / 1e8
#Calculate the ratios used for alphas 
alpha_1_rat = alpha_1/0.03
alpha_2_rat = alpha_2/0.03

#Calculate the schwarchild radius
schwarzchild_radius = bh_func.calculate_schwarschild_radius(BH_mass_kg)

# initialise while loop
#Counter
ii = 0
accretion_interval = 0
accretion_time=[]
kerr_parameter=[]
BH_mass_track=[]


while ii <= 1000 and BH_mass < 1e10:
    #Record values
    if ii == 0:
        accretion_time.append(0)
        kerr_parameter.append(Kerr_spin)
        BH_mass_track.append(BH_mass)
    else:
        accretion_time.append(accretion_time[ii - 1] + accretion_interval)
        kerr_parameter.append(Kerr_spin)
        BH_mass_track.append(BH_mass)

    #Calculate Black hole/accretion disc parameters
    #efficiancy 
    efficiency, x = bh_func.calculate_efficiency(Kerr_spin)

    #Initial Luminocities
    Luminocity_ratio = bh_func.calculate_Eddington_and_actual_Luminocity(BH_mass_8, efficiency)
    #Self gravitating to schwarzchild ratio
    Self_to_schwar_ratio = bh_func.calculate_self_grav_schwarz_ratio(alpha_1_rat, efficiency, Luminocity_ratio, BH_mass_8)
    #Self gravitating mass
    SG_mass_kg = bh_func.calculate_self_gravitiating_mass(alpha_1_rat, efficiency, Luminocity_ratio, BH_mass_8, Self_to_schwar_ratio)
    #Warp radius to schwarzchild radius ratio
    Warp_to_schwar_ratio = bh_func.calculate_warped_to_schwarz_ratio(alpha_1_rat, alpha_2_rat, Kerr_spin, efficiency, Luminocity_ratio, BH_mass_8) 

    #Determine disc angular momentum based on warping
    #Determine Self Gravitating and warped radii from ratios
    SG_radius = Self_to_schwar_ratio * schwarzchild_radius
    warp_radius = Warp_to_schwar_ratio * schwarzchild_radius

    if warp_radius < SG_radius:
        disc_ang_mom = bh_func.calculate_disc_angular_momentum_1(SG_mass_kg, BH_mass_kg, warp_radius)
    else:
        disc_ang_mom = bh_func.calculate_disc_angular_momentum_2(SG_mass_kg, BH_mass_kg, SG_radius)

    #Calculate accretion event time
    accretion_interval = bh_func.calculate_accretion_time(alpha_1_rat, efficiency, Luminocity_ratio, BH_mass_8)
    #Calculate new kerr spin parameter
    #Calculate black hole angular momentum
    BH_ang_mom = bh_func.calculate_black_hole_angular_momentum(BH_mass_kg, Kerr_spin, schwarzchild_radius)
    #Calculate ISCO angular momentum
    ISCO_ang_mom = bh_func.calculate_ISCO_angular_momentum(BH_mass_kg, x)

    #Determine disc co or counter allignment
    #Generate a random number between -1 and 1
    random_number = np.random.randint(1, 2000, 1)
    random_number = (random_number - 1000)/1000

    ang_mom_ratio = -disc_ang_mom/(2*BH_ang_mom)

    if random_number < ang_mom_ratio: #counter alligned
        BH_ang_mom = BH_ang_mom - (SG_mass_kg * ISCO_ang_mom)
    else:
        BH_ang_mom = BH_ang_mom + (SG_mass_kg * ISCO_ang_mom)



    #Update Black hole vales
    #Update mass
    BH_mass_kg = BH_mass_kg + SG_mass_kg
    #Convert to solar masses
    BH_mass = BH_mass_kg / 1.989e30
    #Convert to 10^8 solar mass
    BH_mass_8 = BH_mass / 1e8
    #Calculate the new schwarchild radius
    schwarzchild_radius = bh_func.calculate_schwarschild_radius(BH_mass_kg)

    Kerr_spin = bh_func.calculate_spin_parameter(BH_ang_mom, BH_mass_kg, schwarzchild_radius)

    if ii % 100 == 0:
        text_print = str(ii)
        print(text_print)
        print("\n")

    ii += 1

#Add the final updates to the records
accretion_time.append(accretion_time[ii - 1] + accretion_interval)
kerr_parameter.append(Kerr_spin)
BH_mass_track.append(BH_mass)

accretion_time = np.array(accretion_time)
kerr_parameter = np.array(kerr_parameter)
BH_mass_track = np.array(BH_mass_track)

plt.plot(accretion_time, kerr_parameter)
plt.show()

#plt.plot(accretion_time, math.log10(BH_mass_track) )
#plt.show()