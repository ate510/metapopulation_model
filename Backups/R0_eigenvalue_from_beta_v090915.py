## Author: Anne Ewing
## Date: 07/08/15
## Function: calculate R0 from beta (want R0 to be close to 1.2)

### packages/modules ###
import csv
import sys
import datetime as date
#import networkx as nx
import numpy as np

### local modules ###
sys.path.append('/home/anne/Dropbox/Anne_Bansal_Lab')

# import functions #
import population_parameters as func

# functions #
###################################################
def calc_R_matrix (beta, gamma, alpha):
# calculate R matrix from beta value
    
    # calc contact matrix C
    C = func.calc_contact_matrix(filename_germ_contact_data, filename_germ_pop_data, alpha)
    
    # assign components of C matrix
    C_cc = C.item((0, 0))
    C_ca = C.item((0, 1))
    C_ac = C.item((1, 0))
    C_aa = C.item((1, 1))
    
    # multiply components of matrix C by alpha or (1-alpha)
    mx_11 = (C_cc * alpha)
    mx_12 = (C_ca * alpha)
    mx_21 = (C_ac * (1 - alpha))
    mx_22 = (C_aa * (1 - alpha))
    
    # calc R matrix
    # SB needs to check - might be (beta / (beta + gamma))
    R = ((beta / gamma) * (np.matrix([[mx_11, mx_12], [mx_21, mx_22]])))
    #R = ((beta / (beta + gamma)) * (np.matrix([[mx_11, mx_12], [mx_21, mx_22]])))
    
    return R
    
###################################################    
def calc_R0 (R):
#determine R0 from largest eigenvalue of R matrix

    # calc eigenvalues from R matrix
    # should return in descending order (largest will be first)
    evals, matrix_v = np.linalg.eig(R) # matrix_v --> normalized eigenvectors
    
    # largest eigenvalue is R0
    R0 = max(evals)
       
    return R0   
    
###################################################
if __name__ == "__main__":
    
    # READ POPULATION DATA FROM FILE
    popdata = csv.reader(open('Dropbox/Anne_Bansal_lab/SDI_Data/totalpop_age.csv', 'r'),delimiter = ',')
    dict_popdata, ages, years = func.import_popdata(popdata, 0, 1, 2)
    dict_childpop, dict_adultpop = func.pop_child_adult (dict_popdata, years)
    
    # READ GERMAN CONTACT DATA FROM FILE
    filename_germ_contact_data = 'Dropbox/Anne_Bansal_lab/Contact_Data/polymod_germany_contact_matrix_Mossong_2008.csv'
    filename_germ_pop_data = 'Dropbox/Anne_Bansal_lab/UNdata_Export_2008_Germany_Population.csv'

    # DEFINE DISEASE PARAMETERS
    # if beta/beta+gamma
    #beta = 0.0165
    # if beta/gamma
    beta = 0.015
    gamma = 0.2
    
    # CALCULATE ALPHA
    year = 2010
    alpha = func.calc_alpha(year, dict_childpop, dict_adultpop)

    # CALCULATE R MATRIX
    R = calc_R_matrix (beta, gamma, alpha)
    
    # CALCULATE R0
    R0 = calc_R0(R)
    
    print R0