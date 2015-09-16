## Author: Anne Ewing
## Date: 06/30/15
## Function: functions to calculate alpha, eta, and epsilon

### packages/modules ###
import csv
import sys
import datetime as date
import networkx as nx
import numpy as np



# functions #
###################################################
def import_metropop (datafile, metrocol, popcol):
#create function to import pop of each metro area
# SB will send me data
##Names in script: d_pop_for_metro
##args: (metropop, 0, 1)
    
    dict_metropop, metro_list = {}, []
    for row in datafile:
        metro_id = int(row[metrocol])
        metro_list.append(metro_id)
        pop = float(row[popcol])
        dict_metropop[(metro_id)] = pop
    return dict_metropop, list(set(metro_list))

###################################################
def import_popdata (datafile, yrcol, agecol, popcol):
#create a function that imports national pop data into dictionaries and lists
##Names in script: d_yr_age, ages, years = {}, [], []
##args: (popdata, 0, 1, 2)
  
    dict_popdata, age_list, yr_list = {}, [], []
    for row in datafile:
        year = int(row[yrcol])
	yr_list.append(year)
        age = str.lower(row[agecol])
	age_list.append(age)
        pop = float(row[popcol])
        dict_popdata[(year, age)] = pop
    return dict_popdata, list(set(age_list)), list(set(yr_list))

####################################################
def import_contact_matrix (datafile, age_part_col, age_contact_col, num_contact_col):
# data format - 3 columns: participant age group, contact age group, avg number of daily contacts
    
    dict_contacts, age_list = {}, []
    for row in datafile:
        age_1 = str(row[age_part_col])
        age_list.append(age_1)
        age_2 = str(row[age_contact_col])
        contacts = float(row[num_contact_col])
        dict_contacts[(age_1, age_2)] = contacts
    
    return dict_contacts, list(set(age_list))
    
###################################################
def aggregate_contact_matrix (age_list, dict_contacts):
    
    child = ['05-09', '10-14', '15-19']
    C_ac = ((dict_contacts[age_a, age_c]) * part_a)# don't have part data        
	  	  
###################################################
def pop_child_adult (dict_popdata, years, child_list, adult_list):
# create a function to combine children and adults into 2 groups
##Names in script: d_childpop, d_adultpop = {}, {}
##args: (d_pop_for_yr_age, years, child, adult)
    
    dict_childpop, dict_adultpop = {}, {}
    for y in years:
        childpop = sum([dict_popdata[y, a] for a in child_list])
        #childpop = ((d_yr_age[y, '5-9 years']) + (d_yr_age[y, '10-14 years']) + (d_yr_age[y, '15-19 years']))
        dict_childpop[y] = childpop
        adultpop = sum([dict_popdata[y, a] for a in adult_list])
        #adultpop = ((d_yr_age[y, '20-29 years']) + (d_yr_age[y, '30-39 years']) + (d_yr_age[y, '40-49 years']) + (d_yr_age[y, '50-59 years']) + (d_yr_age[y, '60-69 years']))
        dict_adultpop[y] = adultpop
        
    return dict_childpop, dict_adultpop

###################################################
def pop_child_frac (year, dict_childpop, dict_adultpop):
# determine alpha = the fraction of the U.S. population that is children 
# in functions_chain_binomial_sim_pop_disease_model.py
    # there is a function (zippop_child_frac) to calculate alpha across zip3s
    # may be useful if ever calc alpha per metro area
    # would also need function (import_zip3_popdata)
##Names in script: a = value
##args: (year, d_childpop, d_adultpop)
    
    childpop = dict_childpop[year]
    adultpop = dict_adultpop[year]
    alpha = ((childpop) / (childpop + adultpop))
    
    return alpha
    
###################################################
def assign_demographics (metro_id_list, dict_metropop, alpha):
#assign individuals to metro and to age
      
        dict_metro_for_indv, dict_age_for_indv = {}, {}
        total_pop = sum([dict_metropop[metro_id] for metro_id in metro_id_list])
        total_pop_child = (total_pop * alpha)
        total_pop_adult = (total_pop - total_pop_child)
        total_child_ids = range(1, (total_pop_child + 1)) #range(i, j-1)
        total_adult_ids = range((total_pop_child + 2), (total_pop_child + total_pop_adult + 3))
        for i in total_child_ids:
            dict_age_for_indv[i] = "child"
        for i in total_adult_ids:
            dict_age_for_indv[i] = "adult"
        #for metro_id in metro_id_list:
        #    pop = dict_metropop[metro_id]
        #    childpop = (pop * alpha)
        #    #d[metro_id] = childpop
        #    child_ids = range(1, (childpop + 1))
        #    adultpop = (pop - childpop)
        #    adult_ids = range(1, (adultpop + 1))
        #dict_metro_for_indv[x] = metro_id
        #dict_age_for_indv[x] = "child"
        #dict_age_for_indv[x] = "adult"
        
        return dict_age_for_indv
      
###################################################    
def weighted_avg_q (age_list, d_mean_contacts, d_num_part):
# calc avg q value (avg # of contacts) for one age group
    # from POLYMOD table 1 mean number of contacts for age group
    # and number of participants in each age group

    qxN = [(d_mean_contacts[age] * d_num_part[age]) for age in age_list]
    N_list = [d_num_part[age] for age in age_list]
    avg_q = (sum(qxN))/(sum(N_list)) # divide by total # participants
    
    return avg_q
    
###################################################
def calc_eta (avg_q_child, avg_q_adult):
# calc eta (n - ratio of # of contacts) 
    # from avg # of contacts (q) for each age group 
    
    n = (avg_q_adult / avg_q_child)
    
    return n

####################################################
def read_edgelist_anne (filename):
#import edgelist data without using networkx
   
    G = nx.Graph()
    file = open(filename, 'rU')
    reader = csv.reader(file)
    for row in reader:
        data = row[0].split('\t')
        G.add_edge(int(data[0]),int(data[1]), weight=float(data[2]))
        
    return G

####################################################
def calc_prob_travel (network, alpha, ch_travelers_r, dict_metropop):
# calculate probability of child or adult traveling from metro i to metro j
    
    dict_prob_ch_travel, dict_prob_ad_travel = {}, {}
    r = ch_travelers_r #fraction of travelers who are children
    edges = network.edges()
    for (i, j) in edges:
        w_ij = network[i][j]['weight']
        N_i = dict_metropop[i]
        prob_child_travels = ((r / alpha) * (w_ij / N_i))
        prob_adult_travels = (((1 - r) / (1 - alpha)) * (w_ij / N_i))
        dict_prob_ch_travel[(i, j)] = prob_child_travels
        dict_prob_ad_travel[(i, j)] = prob_adult_travels
    
    return dict_prob_ch_travel, dict_prob_ad_travel
    
####################################################
def calc_contact_matrix (a, n, E, q_c):
#Apollini 2014 Eq 3

    C_cc = ((a - E) / (a **2))
    C_ca = (E / (a * (1 - a)))
    C_ac = (E / (a * (1 - a)))
    C_aa = (((n * (1 - a)) - E) / ((1 - a) **2))
    # matrix
    C_matrix = np.matrix([[C_cc, C_ca], [C_ac, C_aa]])
    C = (C_matrix * q_c)
    
    return C