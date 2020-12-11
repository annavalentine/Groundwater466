
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, math
import imageio

def hi():
    print("hello")

def plot_piper(Ca, Mg, Na, K, SO4, HCO3, CO3, Cl):
    
    #read in our awesome figure
    img = imageio.imread("../PiperDiagramColors.jpg")
    
    Ca_N, Mg_N, Na_N, K_N, SO4_N, HCO3_N, CO3_N, Cl_N, Cat_T, An_T = normalities(Ca, Mg, Na, K, SO4, HCO3, CO3, Cl)
    Ca_P, Mg_P, NaK_P, SO4_P, CO3_P, Cl_P = meq(Ca_N, Mg_N, Na_N, K_N, SO4_N, HCO3_N, CO3_N, Cl_N, Cat_T, An_T)
    
    meqs = [Ca, Mg, Na, K, SO4, HCO3, CO3, Cl]
    norms = [Ca_N, Mg_N, Na_N, K_N, SO4_N, HCO3_N, CO3_N, Cl_N]
    p_meq = [Ca_P, Mg_P, NaK_P, SO4_P, CO3_P, Cl_P]
    
    df(meqs, norms, p_meq, Cat_T, An_T)
    
    #Plot the image
    plt.figure(figsize=(20,20) )
    
    plot_cations(Ca_P, Mg_P, NaK_P)
    plot_anions(Cl_P, SO4_P, CO3_P)
    
    plt.axis('off')

    plt.imshow(img)
    
    
def df(meqs, norms, p_meq, Cat_T, An_T):
    
    Ca_M = 40.078
    Mg_M = 24.305
    Na_M = 22.9 
    K_M = 39.0983 
    SO4_M = 96.06
    HCO3_M = 61.10168
    CO3_M = 60.009
    Cl_M = 35.453
    
    MM = [Ca_M, Mg_M, Na_M, K_M, SO4_M, HCO3_M, CO3_M, Cl_M]
    
    names = ['Ca', 'Mg', 'Na', 'K', 'SO4', 'HCO3', 'CO3', 'Cl']
    cols = ['mg/L', 'Molecular Mass', 'Normality']
    
    df = pd.DataFrame([meqs, MM, norms],  index = cols, columns = names)
    df = df.round(2)
    print(df.to_markdown())
    
    print("")
    print("Total Anions:", An_T)
    print("Total Cations:", Cat_T)
    print("")
    
    names2 = ['Ca', 'Mg', 'Na + K', 'SO4', 'CO3 + HCO3', 'Cl']
    col2 = ['% mEq/L']
    df2 = pd.DataFrame([p_meq],  index = col2, columns = names2)
    df2 = df2.round(2)
    print(df2.to_markdown())
    

# Calculate Normalities:
def normalities(Ca, Mg, Na, K, SO4, HCO3, CO3, Cl):
    Ca_M = 40.078
    Mg_M = 24.305
    Na_M = 22.9 
    K_M = 39.0983 

    SO4_M = 96.06
    HCO3_M = 61.10168
    CO3_M = 60.009
    Cl_M = 35.453
    
    #Calculate Normalities
    #Normality
    Ca_N = (Ca*2) / Ca_M
    Mg_N = (Mg*2) / Mg_M
    Na_N = Na / Na_M
    K_N = K / K_M
    Cat_T = Ca_N + Mg_N + Na_N + K_N 

    SO4_N = (SO4*2) / SO4_M
    HCO3_N = HCO3 / HCO3_M
    CO3_N = (CO3*2) / CO3_M
    Cl_N = Cl / Cl_M
    An_T = SO4_N + HCO3_N + CO3_N + Cl_N
    
    return Ca_N, Mg_N, Na_N, K_N, SO4_N, HCO3_N, CO3_N, Cl_N, Cat_T, An_T
   

def meq(Ca_N, Mg_N, Na_N, K_N, SO4_N, HCO3_N, CO3_N, Cl_N, Cat_T, An_T):
    Ca_P = Ca_N / Cat_T
    Mg_P = Mg_N / Cat_T
    NaK_P = (Na_N + K_N) / Cat_T

    Cl_P = Cl_N / An_T
    SO4_P = SO4_N / An_T
    CO3_P = (HCO3_N + CO3_N) / An_T
    
    return Ca_P, Mg_P, NaK_P, SO4_P, CO3_P, Cl_P

def plot_cations(Ca_P, Mg_P, NaK_P):
    
    #Start with Ca
    Ca_P = (1 - Ca_P)
    Ca_point = (Ca_P*100*3.16) + 120
    
    #starting y
    Ca_y1 = (1.55*(280)) + 124
    Ca_y2 = (1.55*(Ca_point)) + (124)
    
    Ca_diff = Ca_y2 - Ca_y1
    
    #Find Mg_Y
    Mg_y = (-Mg_P*251)+554
    Mg_x_min = (Mg_y-740)/-1.55
    
    #Find Xi
    xi = (Mg_y - (120 - Ca_diff))/1.55
    
    #Na_K
    Na_point = (NaK_P*100*3.16)+120
    
    Na_y1 = (-1.55*(120)) + (740)
    Na_y2 = (-1.55*(Na_point)) + (740)
    
    Na_diff = Na_y2 - Na_y1
    
    #Bound by 
    xb = ((740- Na_diff)- (-732)) / (3.1)
    
    #Plotting
    Na_x = np.linspace(xi, xb, 100)
    Ca_x = np.linspace(Ca_point, xi, 300)
    
    Ca_y = (1.55*(Ca_x) + (120 - Ca_diff))  #Ca_line
    Na_y= (-1.55*(Na_x) + (740 - Na_diff)) #Na Line
    
    
    plt.plot(Na_x, Na_y, linewidth = 4, color = 'black')
    plt.plot(Ca_x, Ca_y, linewidth = 4, linestyle = 'dashed', color = 'black')
    plt.hlines(Mg_y, Mg_x_min, xi, linewidth = 4, linestyle = 'dashed', colors = 'black')
    
    b2 = (740 - Na_diff)
    
    return b2 
    

def plot_anions(Cl_P, SO4_P, CO3_P):
    
    #Find Cl 
    Cl_point = (Cl_P*100*3.16)+508
    
    #starting y
    Cl_y1 = (-1.55*(120)) + (740)
    Cl_y2 = (-1.55*(Cl_point)) + (740)
    
    Cl_diff = Cl_y2 - Cl_y1
    
    #Find SO4 Y
    SO4_y = (-SO4_P*251)+554
    SO4_x_max = (SO4_y+732)/1.55
    
    #Find Xi
    xi = (SO4_y - (740 - Cl_diff))/ -1.55
    
    #CO3
    CO3_P = (1 - CO3_P)
    CO3_point = (CO3_P*100*3.16) + 508
    
    #starting y
    CO3_y1 = (1.55*(280)) + 124
    CO3_y2 = (1.55*(CO3_point)) + (124)
    
    CO3_diff = CO3_y2 - CO3_y1
    
    #Boundary
    xb = ((120-CO3_diff)- 740) / (-3.1)
    
    #Plotting
    CO3_x = np.linspace(xi, xb, 100)
    Cl_x = np.linspace(Cl_point, xi, 300)
    
    CO3_y = (1.55*(CO3_x) + (120 - CO3_diff))  #CO3_line
    Cl_y = (-1.55*(Cl_x) + (740 - Cl_diff)) #Cl Line
    
    
    plt.plot(CO3_x, CO3_y, linewidth = 4, color = 'black')
    plt.plot(Cl_x, Cl_y, linewidth = 4, linestyle = 'dashed', color = 'black')
    plt.hlines(SO4_y, SO4_x_max, xi, linewidth = 4, linestyle = 'dashed', colors = 'black')
    
    b1 = (120 - CO3_diff)
    return b1
 
