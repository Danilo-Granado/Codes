#Numeric method for solving for exit temperature,
#given the known variables below, and the air properties

#Change props.csv according to your substance properties table
#Also pay attention to exponentials in the table (e.g. u*10**7)

#Used both barycentric and krogh interpolations, produced very similar results

from cmath import pi
import pandas as pd
import numpy as np
from scipy.interpolate import krogh_interpolate

#File path for csv with relevant material properties
file_path = "FT Projects/props.csv"
data = pd.read_csv(file_path)

#Known variables (Temperatures in Celsius, corrected to Kelvin)
Tin = 600 + 273.15
Tinf = 25 + 273.15
diam = 0.5
length = 6
vinf = 5
mdot = 0.5

#Properties from csv converted to lists
Temp_list = np.array(data["Temp"].tolist(), dtype=float)
Dens_list = np.array(data["Dens"].tolist(), dtype=float)
Cp_list = np.array(data["Cp"].tolist(), dtype=float)
Kf_list = np.array(data["K"].tolist(), dtype=float)
u_list = np.array(data["u"].tolist(), dtype=float)
Pr_list = np.array(data["Pr"].tolist(), dtype=float)

#Steps:
#1 - Estimate exit temperature (Tout) and calculate average temperature (Tavg)
#1.1 - Find properties for given Tavg
#2 - Calculate Re and Nu (Reynolds and Nusselt's adimensionals)
#3 - Calculate hi (Inner convective coefficient)
#4 - Calculate film temperature (Tfilm)
#4.1 - Find properties for given Tfilm
#5 - Calulate Re and Nu
#6 - Calculate he (Outer convective coefficient)
#7 - Calculate U (Inverse of the sum of the inverses of the termic resistances)
#8 - Verify q1 == q2 (10% tolerance)

#Step 1 - Estimate exit temperature (Tout) (First guess) and calculate Tavg
Tout = (Tin + Tinf)*2/3
Tavg = (Tin + Tout)/2

val = 0

while True:
    #Step 1 - Estimate exit temperature (Tout) and calculate Tavg
    if val == 1:
        Tout = (Tout-273.15)*1.00001 + 273.15
    elif val == -1:
        Tout = (Tout-273.25)*0.99999 + 273.15
    else:
        Tout = Tout

    Tavg = (Tin + Tout)/2

    print("Current Tout = " + str(Tout))
    #Step 1.1 - Find properties for given Tavg
    Kf = krogh_interpolate(Temp_list, Kf_list, Tavg)*10**-3
    u = krogh_interpolate(Temp_list, u_list, Tavg)*10**-7
    Pr = krogh_interpolate(Temp_list, Pr_list, Tavg)

    #Step 2 - Calculate Re and Nu
    Re = (4*mdot)/(u*pi*diam)
    Nu = 0.023*(Re**0.8)*(Pr**0.3)

    #Step 3 - Calculate hi
    hi = (Nu*Kf)/diam

    #Step 4 - Calculate Tfilm
    Tfilm = (Tinf + Tout)/2

    #Step 4.1 - Find properties for given Tfilm 
    Dens = krogh_interpolate(Temp_list, Dens_list, Tfilm)
    Cp = krogh_interpolate(Temp_list, Cp_list, Tfilm)
    Kf = krogh_interpolate(Temp_list, Kf_list, Tfilm)*10**-3
    u = krogh_interpolate(Temp_list, u_list, Tfilm)*10**-7
    Pr = krogh_interpolate(Temp_list, Pr_list, Tfilm)

    #Step 5 - Calculate Re and Nu
    Re = (diam*vinf*Dens)/u
    #print(Re)

    if Re > 0.4 and Re <= 4:
        CNu = 0.989
        m = 0.33
    elif Re > 4 and Re <= 40:
        Cnu = 0.911
        m = 0.385
    elif Re > 40 and Re <= 4000:
        Cnu = 0.683
        m = 0.466
    elif Re > 4000 and Re <= 40000:
        Cnu = 0.193
        m = 0.618
    elif Re > 40000 and Re <= 4000000:
        Cnu = 0.027
        m = 0.805
    else:
        continue

    Nu = Cnu*Re**m*Pr**(1/3)
    #print(Nu)

    #Step 6 - Calculate he
    he = (Nu * Kf)/diam
    #print(he)

    #Calculate U and dTm
    Ur = 1/((1/hi)+(1/he))
    dTm = ((Tin - Tinf) - (Tout - Tinf))/np.log((Tin-Tinf)/(Tout-Tinf))
    
    #Verify q1 == q2 (10% tolerance, but using 0.1% to test results)
    q1 = Ur * dTm * pi * diam * length
    q2 = mdot * Cp * (Tin-Tout)
    err = 0.001

    if (q1-q2)/q2 > err:
        val = -1
        continue
    elif (q1-q2)/q2 < -err:
        val = 1
        continue
    else:
        print("q1 = " + str(q1))
        print("q2 = " + str(q2))
        print("Final Tout = " + str(Tout-273.15))
        break