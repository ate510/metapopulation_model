## Author: Anne Ewing
## Date: 06/25/15
## Function: Outline for project population and disease model

### packages/modules ###
import csv
import sys
import networkx as nx
import numpy as np

### local modules ###
sys.path.append('/home/anne/Dropbox/Anne_Bansal_Lab')


### functions ###

import population_parameters as func
#import chain_binomial as sbfunc

## POPULATION ##

# 1
## pop of metropolitan areas and pop sizes (US)
## format: 6 columns (met_id, pop_size, met_id, pop_size, travelers in one direction, travelers in other direction)

#metropop = csv.reader(open('/home/anne/Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/metedges.csv', 'r'),delimiter = ',')
#
#d_pop_for_metro, metro_ids = func.import_metropop(metropop, 0, 1) #key: metro_id, value: popsize


# 2
## connectivity btwn metro areas & # of travelers btwn them (US)
## format: edgelist - (met_id1, met_id2, # travelers)
# ex: G = networkx.read_edgelist()

metro_edges_travelers = 'Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/air_traffic_edgelist.txt'
G = func.read_edgelist_anne(metro_edges_travelers)
#print G.edges()
#print G[1][2]['weight']

# 3a
## w/in met area, split into children/adults
### use same fraction of child/adults for each metro area
### determine fraction from entire US, fraction children out of only children & adult
### fraction children = alpha (a) = (#ch)/(#ch + #ad)
#### alt option - (find dif fraction for each metro area)

## data ##
# AGE SPLIT #
# POLYMOD contact study (Ref 22 in Apollini 2013)
## Children: (5 - 19) 
## Adults: (20 - 69) 
#popdata
#totalpop_age.csv

# import csv file of pop data
popdata = csv.reader(open('Dropbox/Anne_Bansal_lab/SDI_Data/totalpop_age.csv', 'r'),delimiter = ',')

# import data into dicts
d_pop_for_yr_age, ages, years = func.import_popdata(popdata, 0, 1, 2)

print ages

#group ages into children and adults

#d_childpop, d_adultpop = func.pop_child_adult(d_pop_for_yr_age, years)

# set value of alpha from U.S. pop data
#year = 2010 # decided to use pop from 2010 bc most recent
#a = func.calc_alpha(year, d_childpop, d_adultpop)
#print a # = 0.239 for 2010 with 60-69; was 0.27 without 60-69

#metro = func.assign_metro(a)

#totpop of metros = 10817
# ids of ppl = [0 - 10816]



# 3b
## contact btwn children and adults
## 'C' = equation 3 in Apollini 2014
### alpha (a) = fraction of ch  --> calc in 3a
### (n) = ratio of ad/ch avg # contacts (q_a/q_c)
### (E) = avg fraction of contacts across age groups
#### --> from Euro data in Table 2 in Apollini 2013 + Ref 22

#Apollini 2013, Table 2, "Europe (average values)"
#n = 0.79
E = 0.097 #avg fraction of contacts across age groups

#q_1 = #see additional file of Apollini 2013
# Table 1 Mossong POLYMOD
# weighted avg for children and adults avg # of contacts

#C_ij, C_ji, child_pop, adult_pop = func.weighted_avg_germ_pop()
#
#C_ca = func.simmetrize_contacts(C_ij, C_ji, child_pop, adult_pop)

# now in function
age = [5, 10, 15, 20, 30, 40, 50, 60]
#contacts = [14.81, 18.22, 17.58, 13.57, 14.14, 13.83, 12.30, 9.21] 
#participants = [661, 713, 685, 879, 815, 908, 906, 728]

child = age[0:3]
adult = age[3:8]
#
#avg_q_ch = func.weighted_avg_q(child)
#avg_q_ad = func.weighted_avg_q(adult)
#
#n = func.calc_eta(avg_q_ch, avg_q_ad) # ratio of avg # contacts (adult / child)
##print n # n = 0.752
#
##Ec, Ea = func.calc_epsilon(avg_q_ch, avg_q_ad, C_ca, a, n)
##print Ec
##print Ea
#
#####################################
#
#C = func.calc_contact_matrix(a)
#
#############################
#
### DISEASE ##
#
### Infc_list = [patient_Zero]
### Susc_list = [everyone_else]
### for t (time steps) in (1, 100):
#    ## S --> I ?
#    ## for s in susc:
#        ## infect nodes with prob (1-e^-B(# of infec contacts)) # infected degree = #infected contacts
#    
#    ## I --> R ?
#    ## for i in infc:
#        ## recover with prob u
#        ## (u = 0.1 - infectious period is 10 days)
#    
