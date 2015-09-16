## Author: Anne Ewing
## Date: 06/30/15
## Function: functions to calculate alpha, eta, and epsilon

### packages/modules ###
import csv
import sys
import datetime as date
import networkx as nx #(Mac - no module named networkx)??
import numpy as np



# functions #
####################################################
#def import_metrolist (datafile, metrocol1, metrocol2):
## create function to create list of metros with edge data
#
#    metros_1, metros_2 = [], []
#    for row in datafile:
#        metro1 = int(row[metrocol1])
#        metros_1.append(metro1)
#        metro2 = int(row[metrocol2])
#        metros_2.append(metro2)
#    return list(set(metros_1)), list(set(metros_2))
    
###################################################
def import_metropop (filename, metrocol, popcol):
#create function to import pop of each metro area
##Names in script: d_pop_for_metro
##args: (metropop, 0, 1)
    
    # import US metro pop data
    datafile = csv.reader(open(filename, 'r'),delimiter = '\t')
        
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
    ## to remove 'years' in '10-14 years'
     #   if len(age) > 6:
     #       age = (age[:-6])
	    #age_list.append(age)
     #   else:
     #       age_list.append(age)
        pop = float(row[popcol])
        dict_popdata[(year, age)] = pop
    return dict_popdata, list(set(age_list)), list(set(yr_list))

####################################################
def import_contact_matrix (filename, age_part_col, age_contact_col, num_contact_col):
# data format - 3 columns: participant age group, contact age group, avg number of daily contacts
    
    # import germany contact data
    datafile = csv.reader(open(filename, 'r'),delimiter = ',')
    headers = datafile.next()
    
    dict_contacts, age_list = {}, []
    for row in datafile:
        age_1 = str(row[age_part_col])
        age_list.append(age_1)
        age_2 = str(row[age_contact_col])
        contacts = float(row[num_contact_col])
        dict_contacts[(age_1, age_2)] = contacts
    
    return dict_contacts, list(set(age_list))
            
###################################################
def import_sameage_contact_matrix (datafile, age_part_col, num_contact_col):
# age participant and age contact are same
# use this data to calculate p

    dict_contacts, age_list = {}, []
    for row in datafile:
        age = str(row[age_part_col])
        age_list.append(age)
        contacts = float(row[num_contact_col])
        dict_contacts[(age)] = contacts
        
    return dict_contacts, list(set(age_list))

###################################################
def sort_contact_keys_by_age (dict_contacts, age_list):
# dict_contacts - keys: age1, age2, value: contacts
# grab and sort ages of participants from first key (age1) in dict_contacts

    child_ages = set([key[0] for key in dict_contacts if int(key[0][0:2]) >= ch_1 and int(key[0][0:2]) <= ch_2])
    #child_ages = ['05-09', '10-14', '15-19']
    adult_ages = set([key[0] for key in dict_contacts if int(key[0][0:2]) >= ad_1 and int(key[0][0:2]) <= ad_2])
    #adult_ages = ['30-34', '35-39', '25-29', '40-44', '50-54', '65-69', '55-59', '60-64', '20-24', '45-49']

    return child_ages, adult_ages
    
###################################################
def import_germany_pop_data (filename_germ_pop_data, agecol, popcol):
#organize population data from germany census into dictionary
#keys are a mix of single ages, grouped ages, and other
    
    # import population data
    datafile = csv.reader(open(filename_germ_pop_data, 'r'),delimiter = ',')
    headers = datafile.next()
    
    dict_age_pop = {}
    for row in datafile:
        age = str(row[agecol])
        pop = float(row[popcol])
        dict_age_pop[(age)] = pop
    
    return dict_age_pop 
    
################################################### 
def organize_germ_contact_data(filename_germ_contact_data):
    
    # import germany contact data
    #germ_contact_data = csv.reader(open('Dropbox/Anne_Bansal_lab/Contact_Data/polymod_germany_contact_matrix_Mossong_2008.csv', 'r'),delimiter = ',')
    #headers = germ_contact_data.next()
    # organize contact data into a dictionary and list of ages
    dict_contacts, germ_contact_ages = import_contact_matrix(filename_germ_contact_data, 0, 1, 2)
    # sort keys of contact dictionary by age    
    contact_child_ages, contact_adult_ages = sort_contact_keys_by_age(dict_contacts, germ_contact_ages)
    # make a dictionary to link contact keys to pop keys (formatted differently, same age group breakdown)
    contact_key_dict = {} #key: integer value = first age of group, #value: key for that group for contact dict
    for age in germ_contact_ages:
        new_key = int(age[0:2])
        contact_key_dict[new_key] = age
        
    return dict_contacts, contact_child_ages, contact_adult_ages, contact_key_dict
                        
################################################### 
def organize_germ_pop_data (filename_germ_pop_data, ch_1, ch_2, ad_1, ad_2):
    
    # import population data
    #germ_pop_data = csv.reader(open('Dropbox/Anne_Bansal_lab/UNdata_Export_2008_Germany_Population.csv', 'r'),delimiter = ',')
    #headers = germ_pop_data.next()
   
    # organize population data from germany census into dictionary (key: age)
    dict_germ_pop = import_germany_pop_data(filename_germ_pop_data, 4, 8)
    
    # sort for only grouped keys and sort into child and adult
    group_keys = [key for key in dict_germ_pop if '-' in key] # grab only group keys
    first_key = sorted([int(key[0:2]) for key in group_keys]) # grab only first age in group key in order to assign child v adult
    child_first_key = [key for key in first_key if key >= 5 and key <=15] #first_key[2:5]
    adult_first_key = [key for key in first_key if key >=20 and key <=65] #first_key[5:15]
    
    #once ages id'd as child v. adult, grab full group key again into sorted lists
    child_keys = [key for key in group_keys if int(key[0:2]) in child_first_key]
    adult_keys = [key for key in group_keys if int(key[0:2]) in adult_first_key]
    
    # second dictionary to link pop keys to contact keys using common integers
    all_keys = child_keys + adult_keys
    pop_key_dict = {}
    for x in all_keys:
        new_key = int(x[0:2])
        pop_key_dict[new_key] = x    
        
    return dict_germ_pop, child_first_key, adult_first_key, child_keys, adult_keys, pop_key_dict
                                        
###################################################
def aggregate_contacts (filename_germ_contact_data, filename_germ_pop_data):
# using Germany data

    # import age to contact dictionary, lists of contact keys for children and adults, and integer age to contact formatted age dict
    dict_contacts, child_ages, adult_ages, contact_key_dict = organize_germ_contact_data(filename_germ_contact_data)
        
    # import age to pop dictionary, and dictionary with integer age to pop-formatted age 
    dict_germ_pop, child_first_key, adult_first_key, child_keys, adult_keys, pop_key_dict = organize_germ_pop_data(filename_germ_pop_data, ch_1, ch_2, ad_1, ad_2) 
    
    # calc weighted average of across group contacts (using pop of contacts) for individual groups
    # ie. dict_agg_contacts['05-09'] = # contacts with adults
    dict_agg_contacts = {} #key: age of part, value: number of contacts with opp age group
    for age in child_ages:
        # multiply contacts by pop of contact age group
        contacts_top = sum([((dict_contacts[age, (contact_key_dict[int_age])]) * (dict_germ_pop[(pop_key_dict[int_age])])) for int_age in adult_first_key])
        pop_bottom = sum([dict_germ_pop[pop_age] for pop_age in adult_keys])
        weighted_avg = (contacts_top / pop_bottom)    
        dict_agg_contacts[age] = weighted_avg
    for age in adult_ages:
        contacts_top = sum([((dict_contacts[age, (contact_key_dict[int_age])]) * (dict_germ_pop[(pop_key_dict[int_age])])) for int_age in child_first_key])
        pop_bottom = sum([dict_germ_pop[pop_age] for pop_age in child_keys])
        weighted_avg = (contacts_top / pop_bottom)    
        dict_agg_contacts[age] = weighted_avg
        
    contacts_child = sum([((dict_agg_contacts[(contact_key_dict[int_age])]) * (dict_germ_pop[(pop_key_dict[int_age])])) for int_age in child_first_key])
    child_pop = sum([dict_germ_pop[pop_age] for pop_age in child_keys])
    weighted_avg = (contacts_child / child_pop)   
    C_ij = weighted_avg
    
    contacts_adult = sum([((dict_agg_contacts[(contact_key_dict[int_age])]) * (dict_germ_pop[(pop_key_dict[int_age])])) for int_age in adult_first_key])
    adult_pop = sum([dict_germ_pop[pop_age] for pop_age in adult_keys])
    weighted_avg = (contacts_adult / adult_pop)   
    C_ji = weighted_avg
    
    return C_ij, C_ji, child_pop, adult_pop
    
###################################################
def simmetrize_contacts (C_ij, C_ji, child_pop, adult_pop):
    
    Cca = (0.5 * (((C_ij * child_pop) + (C_ji * adult_pop)) / child_pop))
    
    return Cca
	  	  	  	  
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
def calc_alpha (year, dict_childpop, dict_adultpop):
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
def calc_metro_age_pop (filename_metropop, alpha):
#assign metro pop for child and adult based on alpha

    d_pop_for_metro, metro_ids = import_metropop(filename_metropop, 2, 3) #key: metro_id, value: popsize
    #x = [d_pop_for_metro[met_id] for met_id in metro_ids] # list with populations of each metro area
    #print len(set(x)) # = 225 - means all pops are unique
    #print len(metro_ids) # = 225
    
    d_pop_for_metro_age = {}
    for met_id in metro_ids:
        metro_pop = d_pop_for_metro[met_id]
        child_pop = int(metro_pop * alpha)
        adult_pop = (metro_pop - child_pop)
        d_pop_for_metro_age[(met_id, 'child')] = child_pop
        d_pop_for_metro_age[(met_id, 'adult')] = adult_pop
  
    return d_pop_for_metro_age

###################################################      
def contacts_per_agegroup ():
# return dictionaries with data from Table 1 in POLYMOD (avg # contacts per agegroup)
    
    age = [5, 10, 15, 20, 30, 40, 50, 60]
    #child = age[0:3]
    #adult = age[3:8]
    contacts = [14.81, 18.22, 17.58, 13.57, 14.14, 13.83, 12.30, 9.21] 
    participants = [661, 713, 685, 879, 815, 908, 906, 728]
    
    d_mean_contacts = dict(zip(age, contacts))
    d_num_part = dict(zip(age, participants))
    
    return d_mean_contacts, d_num_part
    
###################################################
#def calc_p ():
#    
#    #import data for same age contacts
#    germ_contact_sameage_data = csv.reader(open('Dropbox/Anne_Bansal_lab/Contact_Data/same_age_polymod_germany_contact_matrix_Mossong_2008.csv', 'r'),delimiter = ',')
#    headers = germ_contact_sameage_data.next()
#    
#    dict_sameage_contacts, ages = import_sameage_contact_matrix(germ_contact_sameage_data, 0, 1)
#    
#    for age in ages:
#        
#        pxN = [(dict_sameage_contacts[age] * dict_age_pop[age]) for age in ages]
#    
#    return p_child, p_adult
         
###################################################    
def weighted_avg_q (age_list):
# calc avg q value (avg # of contacts) for one age group
    # from POLYMOD table 1 mean number of contacts for age group
    # and number of participants in each age group

    d_mean_contacts, d_num_part = contacts_per_agegroup()

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
def calc_epsilon (avg_q_child, avg_q_adult, Cca, alpha, n):

    E_c = ((1 / avg_q_child) * Cca)
    #E_a = ((1 / avg_q_adult) * Cca)
    E = (E_c * alpha)
    #E_test = (E_a * n * (1 - alpha))
    
    return E

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
        w_ij_daily = (w_ij / 365) # beta and gamma are in # of days
        N_i = dict_metropop[i]
        N_j = dict_metropop[j]
        prob_child_travels_i_j = ((r / alpha) * (w_ij_daily / N_i))
        prob_adult_travels_i_j = (((1 - r) / (1 - alpha)) * (w_ij_daily / N_i))
        prob_child_travels_j_i = ((r / alpha) * (w_ij_daily / N_j))
        prob_adult_travels_j_i = (((1 - r) / (1 - alpha)) * (w_ij_daily / N_j))
        dict_prob_ch_travel[(i, j)] = prob_child_travels_i_j
        dict_prob_ad_travel[(i, j)] = prob_adult_travels_i_j
        dict_prob_ch_travel[(j, i)] = prob_child_travels_j_i
        dict_prob_ad_travel[(j, i)] = prob_adult_travels_j_i
    
    return dict_prob_ch_travel, dict_prob_ad_travel
    
####################################################
def calc_contact_matrix (filename_germ_contact_data, filename_germ_pop_data, a):
#Apollini 2014 Eq 3

    age = [5, 10, 15, 20, 30, 40, 50, 60]
    child = age[0:3]
    adult = age[3:8]

    #calc avg # of contacts for each age group
    q_c = weighted_avg_q(child)
    q_a = weighted_avg_q(adult)
    
    #calc ratio of avg # of contacts
    n = calc_eta(q_c, q_a)
    
    #calc contacts with opp age group to calculate epsilon (E - avg frac of contacts across age groups)
    C_ij, C_ji, child_pop, adult_pop = aggregate_contacts(filename_germ_contact_data, filename_germ_pop_data)
    Cca = simmetrize_contacts(C_ij, C_ji, child_pop, adult_pop)
    E = calc_epsilon(q_c, q_a, Cca, a, n)

    # calc components of matrix C
    C_cc = ((a - E) / (a **2))
    C_ca = (E / (a * (1 - a)))
    C_ac = (E / (a * (1 - a)))
    C_aa = (((n * (1 - a)) - E) / ((1 - a) **2))
    # matrix
    C_matrix = np.matrix([[C_cc, C_ca], [C_ac, C_aa]])
    C = (C_matrix * q_c)
    
    return C

####################################################    
def calc_contact_matrix_pqa (p_c, p_a, q_c, q_a, a):
# calculate contact matrix C using eq. 1 in Apollini 2014
   
    C_cc = ((p_c * q_c) / a)
    C_ca = (((1 - p_a) * q_a) / a)
    C_ac = (((1 - p_c) * q_c) / (1 - a))
    C_aa = ((p_a * q_a) / (1 - a))
    C = np.matrix([[C_cc, C_ca], [C_ac, C_aa]])
    
    return C
    
####################################################
# flexible parameters
# change these

#age that starts first and last group of each age
ch_1 = 05
ch_2 = 15
ad_1 = 20
ad_2 = 65

####################################################
if __name__ == "__main__":
    
    filename_metropop = 'Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/metedges.txt'
    filename_germ_contact_data = 'Dropbox/Anne_Bansal_lab/Contact_Data/polymod_germany_contact_matrix_Mossong_2008.csv'
    filename_germ_pop_data = 'Dropbox/Anne_Bansal_lab/UNdata_Export_2008_Germany_Population.csv'
    
    C_ij, C_ji, child, adult = aggregate_contacts(filename_germ_contact_data, filename_germ_pop_data)
    print C_ij
    

    
