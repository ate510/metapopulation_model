import networkx as nx # to work with the population contact network
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
import operator

### local modules ###
sys.path.append('/home/anne/Dropbox/Anne_Bansal_Lab')

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

## UNIT TEST ##
unit_test = 'no_travel'
## NO TRAVEL ##
ch_travelers_r = 0.0 # fraction of children who travel
ad_travelers_s = 0

# CONTACT MATRIX
q_c, q_a, p_c, p_a, _, _ = pop_func.calc_p(filename_germ_within_group_contact_data, filename_germ_pop_data, filename_germ_all_contact_data)
C = pop_func.calc_contact_matrix_pqa(p_c, p_a, q_c, q_a, alpha)
#print C
                        
# DEFINE DISEASE PARAMETERS
R0 = 1.2
gamma = 0.5 # recovery rate based on (1/gamma) day infectious period
beta = 0.037
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

# RUN EPIDEMIC SIMULATIONS
#average_epidemic_size = bin_func.chain_binomial_monte_carlo_unit_tests(unit_test, R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, C)
#average_epidemic_size, adult_epi_size, child_epi_size, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult = bin_func.chain_binomial_monte_carlo_unit_tests_csv(unit_test, R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, C)
#TEST#
average_epidemic_size, adult_epi_size, child_epi_size, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult = bin_func.chain_binomial_monte_carlo_unit_tests_csv_test(unit_test, R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, C)
bin_func.write_csv_file_no_travel(num_metro_zeros, num_child_zeros, num_adult_zeros, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult)

# OUTPUT AR
print "\nAverage Large Epidemic Size = ", round(100*average_epidemic_size,2), '%.\n'
print "\nAverage Adult Epidemic Size = ", round(100*adult_epi_size,2), '%.\n'
print "\nAverage Child Epidemic Size = ", round(100*child_epi_size,2), '%.\n'
