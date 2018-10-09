# -*- coding: utf-8 -*-
"""
This script is designed to evaluate the aggregation kinetics of amyloid-beta
in physiological conditions by fitting experimental data to the Finke-Watzky 
2-step aggregation model. 

The parameters obtained will then be injected into a COMSOL model for 
amyloid-beta production, aggregation and degradation in the brain and inside a 
microfluidic device.

Created on Fri Jun  1 17:45:19 2018

@author: Thomas Galeon
"""

import csv
import scipy
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score
from matplotlib.ticker import FuncFormatter
import os


################################## VARIABLES ##################################

# INITIAL CONDITIONS
# Experimental initial condition cO in uM
c0 = 50
# Physiological concentration cphysio in nM
cphysio = 1


# DEFAULT AB-42 KINETIC RATES FOR PARAMETER SWEEP

# Nucleation rate [h^-1]
k1_n = 4.20*10**(-2)
# Clearance rate [h^-1]
kc_n = 7.6*10**(-2)
# Production rate [nM.h^-1]
kp_n = 6.7*10**(-3)

# Production rate in Familial Alzheimer's disease [nM.h^-1]
kp_FAD = 13.4*10**(-3)
# Clearance rate in Sporadic Alzheimer's disease [h^-1]
kc_SAD = 5.3*10**(-2)
# Production rate in Sporadic Alzheimer's disease [nM.h^-1]
kp_SAD = 6.6*10**(-3)


# DISPLAY SETTINGS
width = 3
height = width / 1.618

############################## MAIN FUNCTION ##################################

def main():
    # Set plot settings for optimal display
    plt.rc('font', serif='Arial')    
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=7)
    plt.rc('legend', fontsize=5) 
    
    # Uncomment desired portion of code to obtain images
    # CTRL+4 to comment block, CTRL+5 to uncomment block
    
    #   AGGREGATION KINETICS FIT TO EXPERIMENTAL DATA  
    #AB42
    experimental_fit(0,2)
    #AB40
    experimental_fit(3,0)

# =============================================================================
#     #   AGGREGATION KINETICS INSIDE THE BRAIN
#     global cphysio
#     #   AB42
#     cphysio = 0.1
#     brain(0,2)
#     #   AB40
#     cphysio = 1
#     brain(3,0)
# =============================================================================
    
    
# =============================================================================
#     #   PARAMETER SWEEP
#     clearance()
#     production()
#     kinetics()
# =============================================================================

    return

############################ DATA ACQUISITION #################################

# Open data file
my_path = os.path.abspath(os.path.dirname(__file__))
file = os.path.join(my_path, 'DATA\CSV_Fit_Neurotoxicity_of_AD_AB_peptides_Fig_1A.csv')
data_file = open(file)

def csv_list(file, column):
    '''Builds two x and y arrays from the first two columns of the csv file '''
    read = csv.reader(file)
    x = []
    y = []
    next(read)
    for row in read:
        if row[column]:
            x.append(float(row[column]))
            y.append(float(row[column+1]))
    file.seek(0)
    return np.array(x),np.array(y)

############################## DATA ANALYSIS ##################################

def aggregation(t, k1, k2):
    '''Returns fraction of aggregated amyloid-beta at time t given the nucleation
    rate k1 and the autocatalytic growth rate k2 in experimental conditions'''
    a0 = c0
    return (1 - (k1/(k2*a0)+1)/(1+k1*np.exp((k1+k2*a0)*t)/(k2*a0)))

def brain_agg(t, k1, k2):
    '''Returns fraction of aggregated amyloid-beta at time t given the nucleation
    rate k1 [h-1] and the autocatalytic growth rate k2 [uM.h-1] inside the brain'''
    # cphysio is in nM > convert to uM
    a0 = cphysio/1000
    return (1 - (k1/(k2*a0)+1)/(1+k1*np.exp((k1+k2*a0)*t)/(k2*a0)))

def FW(x_data, y_data):
   ''' Identifies optimal parameters to fit data with FW mechanism '''
   return(scipy.optimize.curve_fit(aggregation, x_data, y_data, bounds=(0, 1)))
   
def experimental_fit(column, peptide):
    ''' Diplays experimental fit to data contained in the column. Peptide is 
    equal to 0 if plotting AB40 and 2 if plotting AB42.'''
    
    # Fitting data
    x_data, y_data = csv_list(data_file, column)
    popt, pcov = FW(x_data, y_data)
    y_filt = aggregation(x_data,*popt)
    r2 =r2_score(y_data, y_filt)
    
    # Extracting aggregation kinetic rates from results
    k1 = ("{:.2e}".format(popt[0]))
    k2 = ("{:.2e}".format(popt[1]))
    
    
    # Plotting data
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=.15, bottom=.2, right=.95, top=.97)   
   
    ax.set_ylabel("% Aggregated Protein")
    ax.set_xlabel("Time (hours)")
    
    # AB40
    if column == 0:
        color = 'r'
        plt.text(2.9, 0.2, r'$k_1 = $' + str(k1) + r' $h^{-1}$' + '\n' + '$k_2 = $' + str(k2) + ' $\mu M^{-1} h^{-1}$ \n\n$[M]_0 = %s $ $\mu M $' %c0, fontname='Arial', fontsize=6)
    
    # AB42
    else:
        color = 'b'
        plt.text(0, 0.2, r'$k_1 = $' + str(k1) + r' $h^{-1}$' + '\n' + '$k_2 = $' + str(k2) + ' $\mu M^{-1} h^{-1}$ \n\n$[M]_0 = %s $ $\mu M $' %c0, fontname='Arial', fontsize=6)
  
    plt.plot(x_data, y_data, '%s s' %color, markersize = 2.5, alpha = 0.5, label = r'$A\beta_{4%s}$ Data' %peptide)
    plt.plot(x_data, y_filt, '%s' %color, linewidth = 1.0, label = r'$A\beta_{4%s}$ Fit, $R^{2}=$ %.4f' %(peptide,r2))  
    plt.legend()
    #plt.title('AB4%s Aggregation Kinetics fitted with FW mechanism' %peptide)
    
    # Save and display image
    fig.set_size_inches(width, height)
    fig.savefig('AB4%s Experimental Aggregation' %peptide + '.png', dpi=1200)
    plt.show()
    
    return

def brain(column, peptide):
    """ Displays aggregations kinetics inside the brain by extrapolating 
    kinetic rates from experimental data contained in column. Peptide is 
    equal to 0 if plotting AB40 and 2 if plotting AB42. """
    
    # Fitting experimental data
    x_data, y_data = csv_list(data_file, column)
    popt, pcov = FW(x_data, y_data)
    k1 = ("{:.2e}".format(popt[0]))
    k2 = ("{:.2e}".format(popt[1]))
    
    
    # The parameter c is a constant adjusting for the different aggregation
    # timescales for AB40 and AB42
    # AB 40
    if column == 0:
        c = 1       
    #AB42
    else:
        c = 500
        
    # Calculating aggregate fraction as a function of time    
    time = np.linspace(0, 150*c, 100)
    y_ded = brain_agg(time, *popt)
    
    
    # Plotting figures
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=.15, bottom=.2, right=.95, top=.97)  
    
    if column == 0:
        color = 'r'
        ax.set_xlabel("Time (hours)")   
    else:
        color = 'b'
        #convert to days
        h_to_y = 24*365.25
        time = time/h_to_y
        ax.set_xlabel("Time (years)")
        
    ax.set_ylabel("% Aggregated Protein")
     
    plt.plot(time, y_ded, '%s' %color, linewidth = 2.0, label = r'$A\beta_{4%s}$' %peptide)  

    if column == 0:
        plt.text(75, 0.2, r'$k_1 = $' + str(k1) + ' $h^{-1}$' + '\n' + '$k_2$ = ' + str(k2) + ' $\mu M^{-1} h^{-1}$ \n\n$[M]_0 = %s $ nM' %cphysio, fontname = 'Arial', fontsize=6)
    else:
        plt.text(4, 0.2, r'$k_1 = $' + str(k1) + ' $h^{-1}$' + '\n' + '$k_2$ = ' + str(k2) + ' $\mu M^{-1} h^{-1}$ \n\n$[M]_0 = %s $ nM' %cphysio, fontname = 'Arial', fontsize=6)
    
    plt.legend()
   
     #plt.title('AB4%s Aggregation Kinetics fitted with FW mechanism' %peptide)
    
    # Save and display image
    fig.set_size_inches(width, height)
    fig.savefig('AB4%s Physiological Aggregation' %peptide + '.png', dpi=1200)
    plt.show()
    return



def clearance():
    """" Parameter sweep for clearance rate kc using default normal AB42 parameters"""
    # kc is the interval considered
    kc = np.linspace(4*10**(-2),9*10**(-2),1000)
    
    # Aggregate fraction calculated analytically
    a_n = 3*k1_n*kp_n/(kc*(kc+k1_n))
    # Changing units to uM
    a_n = a_n*1000
    
    
    # For a parameter sweep for clearance rate using SAD AB42 parameters:
    #a_SAD = 3*k1*kp_SAD/(kc*(kc+k1))
    #a_SAD = a_SAD*1000
    
    # Plotting figure
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=.15, bottom=.2, right=.95, top=.97)  
    ax.set_ylabel(r'Aggregated $A\beta_{42}$ in ISF (pM)')
    ax.set_xlabel('Clearance rate ($h^{-1}$)')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.0%}'.format(y))) 
    
    plt.plot(kc, a_n, 'black', linewidth = 1.0)
    ax.axvline(x=kc_n, color = 'g', linestyle='dashed', alpha=0.5)
    ax.text(x=kc_n + 0.001, y=200, s='Normal', alpha=1, color='g', fontsize = 6)
    
    ax.axvline(x=kc_SAD, color = 'r', linestyle='dashed', alpha=0.5)
    ax.text(x=kc_SAD + 0.001, y=200, s='SAD', alpha=1, color='r', fontsize = 6)
    
    #plt.text(0.06, 1250, 'k1 = ' + str(k1) + ' $h^{-1}$\n$k_{prod}$ = ' + str(kp_n*1000) + 'pM $h^{-1}$', fontname='Arial', fontsize=6)
     
    # Save and display image
    fig.set_size_inches(width, height)
    fig.savefig('Variations in clearance rate.png', dpi=1200)
    plt.show()
    return
    
def production():
    """" Parameter sweep for production rate kp using default normal AB42 parameters"""
    # kp is the interval considered
    kp = np.linspace(4*10**(-3),16*10**(-3),1000)
    
    # Aggregate fraction calculated analytically
    a_n = 3*k1_n*kp/(kc_n*(kc_n+k1_n))
    # Changing units to uM
    a_n = a_n*1000
    kp = kp * 1000
    
    # Plotting figure
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=.15, bottom=.2, right=.95, top=.97)  
    ax.set_ylabel(r'Aggregated $A\beta_{42}$ in ISF (pM)')
    ax.set_xlabel('Production rate ($pM.h^{-1}$)') 
    
    plt.plot(kp, a_n, 'black', linewidth = 1.0)
    ax.axvline(x=kp_n*1000, color = 'g', linestyle='dashed', alpha=0.5)
    ax.text(x=kp_n*1000 + 0.2, y=150, s='Normal', alpha=1, color='g', fontsize = 6)
    
    ax.axvline(x=kp_FAD*1000, color = 'r', linestyle='dashed', alpha=0.5)
    ax.text(x=kp_FAD*1000 + 0.2, y=150, s='FAD', alpha=1, color='r', fontsize = 6)

    # Save and display image      
    fig.set_size_inches(width, height)
    fig.savefig('Variations in production rate.png', dpi=1200)
    plt.show()
    return

def kinetics():
    """" Parameter sweep for production rate k1 using default normal AB42 parameters"""
    # k1_v is the interval considered
    k1_v = np.linspace(2*10**(-5),10*10**(-2))
    
    # Aggregate fraction calculated analytically
    a_n = 3*k1_v*kp_n/(kc_n*(kc_n+k1_v))
    # Changing units to uM
    a_n = a_n*1000
    
    
    # Plotting figure
    fig, ax = plt.subplots()
    fig.subplots_adjust(left=.15+0.01, bottom=.2, right=.95, top=.97)  
    ax.set_ylabel(r'Aggregated $A\beta_{42}$ in ISF (pM)')
    ax.set_xlabel('Nucleation rate ($h^{-1}$)') 
    
    plt.loglog(k1_v, a_n, 'black', linewidth = 1.0)
    ax.axvline(x=k1_n, color = 'g', linestyle='dashed', alpha=0.5)
    ax.text(x=k1_n - 0.033, y=1, s='Baseline' + '\n' + r'      $A\beta_{42}$', alpha=1, color='g', fontsize = 6)
    
    k1_40 = 4.20*10**(-5)
    ax.axvline(x=k1_40, color = 'b', linestyle='dashed', alpha=0.5)
    ax.text(x=k1_40 + 0.00001, y=1.8, s=r'$A\beta_{40}$', alpha=1, color='b', fontsize = 6)

    # Save and display image         
    fig.set_size_inches(width, height)
    fig.savefig('Variations in nucleation rate.png', dpi=1200)
    plt.show()
    return

main()

data_file.close()