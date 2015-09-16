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
def import_germany_pop_data (datafile, agecol, popcol):
    
    dict_age_pop = {}
    for row in datafile:
        age = str(row[agecol])
        pop = float(row[popcol])
        dict_age_pop[(age)] = pop
    
    return dict_age_pop
            
###################################################
def aggregate_contact_matrix ():
# using Germany data

    germ_pop_data = csv.reader(open('Dropbox/Anne_Bansal_lab/UNdata_Export_2008_Germany_Population.csv', 'r'),delimiter = ',')
    headers = germ_pop_data.next()
    
    dict_age_pop = import_germany_pop_data(germ_pop_data, 4, 8)
    group_keys = [key for key in dict_age_pop if '-' in key] # grab only group keys
    first_key = sorted([int(key[0:2]) for key in group_keys])
    child_first_key = [key for key in first_key if key >= 5 and key <=15] #first_key[2:5]
    adult_first_key = [key for key in first_key if key >=20 and key <=65] #first_key[5:15]
    child_keys = [key for key in group_keys if int(key[0:2]) in child_first_key]
    adult_keys = [key for key in group_keys if int(key[0:2]) in adult_first_key]
    
    #dict_agg_pop = {}
    #dict_agg_pop['c'] = sum([dict_age_pop[c] for c in child_keys])
    #dict_agg_pop['a'] = [dict_age_pop[a] for a in adult_keys]
    germ_contact_data = csv.reader(open('Dropbox/Anne_Bansal_lab/Contact_Data/polymod_germany_contact_matrix_Mossong_2008.csv', 'r'),delimiter = ',')
    headers = germ_contact_data.next()
    
    dict_contacts, ages = import_contact_matrix(germ_contact_data, 0, 1, 2)
    
    #contacts_59 = (dict_contacts[59, 
    
    #child_contacts = (dict_contacts[59, 2024])*(dict_age_pop[59])
    
    #child = ['05-09', '10-14', '15-19']
    #C_ij = ((dict_contacts[age_a, age_c]) * part_a)# don't have part data        
	  	  
###################################################
def pop_child_adult (dict_popdata, years):
# create a function to combine children and adults into 2 groups
##Names in script: d_childpop, d_adultpop = {}, {}
##args: (d_pop_for_yr_age, years, child, adult)
    
    dict_childpop, dict_adultpop = {}, {}
    child_list = ['5-9 years', '10-14 years', '15-19 years']
    adult_list = ['20-29 years', '30-39 years', '40-49 years', '50-59 years', '60-69 years'] # added 60-69 bc of feedback EL received
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
        total_adult_ids = range((total_pop_child + 1), (total_pop_child + total_pop_adult + 3))
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
        
        return total_pop, dict_age_for_indv

###################################################      
def contacts_per_agegroup (age, contacts, participants):
# return dictionaries with data from Table 1 in POLYMOD (avg # contacts per agegroup)
    
    d_mean_contacts = dict(zip(age, contacts))
    d_num_part = dict(zip(age, participants))
    
    return d_mean_contacts, d_num_part
    
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
#def calc_epsilon



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