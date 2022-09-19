import math

#Define constants

#Gravitational constant
GRAVITY = 6.6743e-11 #SI units m^3kg^-1s^-2

#Speed of light in a vacuum
LIGHTSPEED = 299792458 #m s^-1

def calculate_schwarschild_radius(Mass):
    radius_s = (2 * GRAVITY * Mass) / (LIGHTSPEED**2)
    return radius_s

def calculate_black_hole_angular_momentum(Mass, a, Schwarz_rad):
    ang_mom = (1 / math.sqrt(2)) * Mass * a * math.sqrt(GRAVITY * Mass * Schwarz_rad)
    return ang_mom

def calculate_efficiency(a):
    #Define start and end points
    x_1 = 0
    x_2 = 10

    #Initialise while loop counter
    count = 0
    
    while count <= 100:
        #Find the midpoint
        x_mid = (x_1 + x_2)/2

        #Calculate solve to find x value
        eqn_part_1 = (x_mid**(1/2))/3
        eqn_part_2 = 4 - (3*x_mid - 2)**(1/2)
        a_calc = eqn_part_1*eqn_part_2

        if  isinstance(a_calc, complex):
            a_calc = a_calc.real

        if a_calc  == a:
            break
        elif a_calc > a:
            x_1 = x_mid
        else:
            x_2 = x_mid
        
        
        count += 1

    #calculate efficiency
    efficiency = 1-(1-(2/(3*x_mid)))**(1/2)

    return efficiency, x_mid

def calculate_Eddington_and_actual_Luminocity(Mass_8, efficiency):
    Eddington_luminocity = 1.4e46 * Mass_8
    Actual_luminocity = efficiency * Eddington_luminocity

    luminocity_ratio = Actual_luminocity/(0.1*Eddington_luminocity)

    return luminocity_ratio

def calculate_self_grav_schwarz_ratio(alpha_1_rat, efficiency, luminocity_ratio, mass_8):
    part_1 = alpha_1_rat ** (14/27)
    part_2 = (efficiency/0.1)**(2/27)
    part_3 = luminocity_ratio**(-8/27)

    ratio = 1.13e3 * part_1 * part_2 * part_3 * (mass_8 ** (-26/27))

    return ratio

def calculate_self_gravitiating_mass(alpha_1_rat, efficiency, luminocity_ratio, mass_8, self_grav_schwar_ratio):
    part_1 = alpha_1_rat**(-4/5)
    part_2 = (efficiency/0.1)**(-3/5)
    part_3 = luminocity_ratio**(3/5)

    SG_mass_grams = 2.94e34 * part_1 * part_2 * part_3 * (mass_8**(11/5)) * (self_grav_schwar_ratio**(7/5))

    #Convert to kg
    SG_mass_kg = SG_mass_grams/1000

    return SG_mass_kg

def calculate_warped_to_schwarz_ratio(alpha_1_rat, alpha_2_rat, a, efficiency, luminocity_rat, mass_8):
    part_1 = alpha_1_rat**(1/8)
    part_2 = alpha_2_rat**(-5/8)
    part_3 = a**(5/8)
    part_4 = (efficiency/0.1)**(1/4)
    part_5 = luminocity_rat**(-1/4)
    part_6 = mass_8**(1/8)

    ratio = 990 * part_1 * part_2 * part_3 * part_4 * part_5 * part_6    

    return ratio

def calculate_disc_angular_momentum_1(self_grav_mass, BH_mass, warp_radius):
    disc_ang_mom = self_grav_mass * ((GRAVITY * BH_mass * warp_radius)**(1/2))
    return disc_ang_mom

def calculate_disc_angular_momentum_2(self_grav_mass, BH_mass, self_grav_radius):
    disc_ang_mom = self_grav_mass * ((GRAVITY * BH_mass * self_grav_radius)**(1/2))
    return disc_ang_mom

def calculate_accretion_time(alpha_1_rat, efficiency, luminocity_rat, mass_8):
    part_1 = alpha_1_rat**(-2/27)
    part_2 = (efficiency/0.1)**(22/27)
    part_3 = luminocity_rat**(-22/27)

    accretion_time_years = 1.12e6 * part_1 * part_2 * part_3 *(mass_8**(-4/27))

    return accretion_time_years

def calculate_ISCO_angular_momentum(bh_mass_kg, x):
    part_1 = 2/(3*math.sqrt(3))
    part_2 = (GRAVITY * bh_mass_kg)/LIGHTSPEED
    part_3 = 1 + 2*math.sqrt(3 * x - 2)

    ISCO_ang_mom = part_1 * part_2 * part_3

    return ISCO_ang_mom

def calculate_spin_parameter(BH_ang_mom, BH_mass_kg, schwarzchild_radius):
    top_bit = math.sqrt(2) * BH_ang_mom
    bottom_bit = BH_mass_kg * math.sqrt(GRAVITY * BH_mass_kg * schwarzchild_radius)

    new_spin_param = top_bit / bottom_bit

    return new_spin_param