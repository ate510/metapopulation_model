import networkx as nx # to work with the population contact network
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
import operator

### local modules ###
#sys.path.append('/home/anne/Dropbox/Anne_Bansal_Lab')

### functions ###
import population_parameters as pop_func
import chain_binomial as bin_func
    
### program ###
###################################################  
# import data #

# READ METRO NETWORK FROM FILE
filename_metropop = 'Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/metedges.txt'
d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
filename_air_network = 'Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/air_traffic_edgelist.txt'
air_network = bin_func.read_edgelist_anne(filename_air_network)

# READ US population data
us_popdata = csv.reader(open('Dropbox/Anne_Bansal_lab/SDI_Data/totalpop_age.csv', 'r'),delimiter = ',')
dict_popdata, ages, years = pop_func.import_popdata(us_popdata, 0, 1, 2)
dict_childpop, dict_adultpop = pop_func.pop_child_adult (dict_popdata, years)

# READ Germany contact data
#filename_germ_contact_data = 'Dropbox/Anne_Bansal_lab/Contact_Data/polymod_germany_contact_matrix_Mossong_2008.csv'
filename_germ_within_group_contact_data = 'Dropbox/Anne_Bansal_lab/Contact_Data/within_group_polymod_germany_contact_matrix_Mossong_2008.csv'
filename_germ_all_contact_data = 'Dropbox/Anne_Bansal_lab/Contact_Data/all_ages_polymod_germany_contact_matrix_Mossong_2008.csv'

# READ Germany population data
filename_germ_pop_data = 'Dropbox/Anne_Bansal_lab/UNdata_Export_2008_Germany_Population.csv'
   
     
# DEFINE POPULATION PARAMETERS
year = 2010
alpha = pop_func.calc_alpha(year, dict_childpop, dict_adultpop)
d_metro_age_pop = pop_func.calc_metro_age_pop(filename_metropop, alpha)
ch_travelers_r = 0.0 # fraction of children who travel
ad_travelers_s = 1

# CONTACT MATRIX
q_c, q_a, p_c, p_a, _, _ = pop_func.calc_p(filename_germ_within_group_contact_data, filename_germ_pop_data, filename_germ_all_contact_data)
params_to_calc_C = []
params_to_calc_C.append(p_c)
params_to_calc_C.append(p_a)
params_to_calc_C.append(q_c)
params_to_calc_C.append(q_a)
params_to_calc_C.append(alpha)
C = pop_func.calc_contact_matrix_pqa(p_c, p_a, q_c, q_a, alpha)
#print C
                        
# DEFINE DISEASE PARAMETERS
R0 = 1.2
gamma = 0.5 # recovery rate based on (1/gamma) day infectious period
test_gammas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
beta = 0.037
test_betas = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
#beta = calculate_beta(R0, gamma, air_network)
#beta = 0
#beta = 0.037 #gamma 0.5
#jbeta = 0.005
num_metro_zeros = 1 # set how many metros to select patients from to start with
num_child_zeros = 1
num_adult_zeros = 0
time_end = 150

# DEFINE TRAVEL PARAMETERS
theta_susc = 1
theta_infc = 1
theta_recv = 1
#d_theta = {}
#d_theta[(s)] = 1
#d_theta[(i)] = 1
#d_theta[(r)] = 1

# DEFINE EXP PARAMETERS
manip_exp_params = []
experiment = 'yes' #turn exp on or off
manip_exp_params.append(experiment)
intervention = 'red_C_cc'
manip_exp_params.append(intervention)
timing = 'real_time'
length_before = 7 # days
length_after = 7 # days
manip_exp_params.append(timing)
manip_exp_params.append(length_before)
manip_exp_params.append(length_after)
beta_exp = 0.02
manip_exp_params.append(beta_exp)

# RUN EPIDEMIC SIMULATIONS
#num_sims = 250 # if debugging, reduce this number to something small like 10
#num_sims = 100
#average_epidemic_size = bin_func.chain_binomial_monte_carlo(R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_sims, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, C) 

# plot beta v epidemic size
bin_func.beta_test(test_betas, gamma, R0, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params)

#plot gamma v epidemic size
#bin_func.gamma_test(beta, test_gammas, R0, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params)
