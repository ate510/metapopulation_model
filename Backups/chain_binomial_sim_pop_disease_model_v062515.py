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

import functions_chain_binomial_sim_pop_disease_model_v062515 as func
#import chain_binomial as sbfunc

## POPULATION ##

# 1
## pop of metropolitan areas and pop sizes (US)
## format: 2 columns (met_id, pop_size)

#metropop = csv.reader(open('/home/anne/Dropbox/Anne_Bansal_lab/SDI_Data/_______.csv', 'r'),delimiter = ',')
#
#d_pop_for_metro = func.import_metropop(metropop, 0, 1)

# 2
## connectivity btwn metro areas & # of travelers btwn them (US)
## format: edgelist - (met_id1, met_id2, # travelers)
# ex: G = networkx.read_edgelist()

metro_edges_travelers = 'Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/air_traffic_edgelist.txt'
G = func.read_edgelist_anne(metro_edges_travelers)
print G.edges()
#contact_network = sbfunc.read_contact_network(metro_edges_travelers)

#metro_edges_travelers = csv.reader(open('/home/anne/Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/air_traffic_edgelist.txt', 'rb'), delimiter='\t')
#metro_edges_travelers = open('/home/anne/Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/air_traffic_edgelist.txt')
#
#
#for edge in metro_edges_travelers:
#    print edge
#
#metro_edges_travelers.close()
#
#G = nx.read_edgelist('/home/anne/Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/air_traffic_edgelist.txt', delimiter='\t', data=True)
#G = nx.read_edgelist(metro_edges_travelers, data=True)
#print G.nodes()
#G = nx.read_edgelist(metro_edges_travelers, nodetype=int, data=True)
#G = nx.read_weighted_edgelist(metro_edges_travelers, delimiter=None, create_using=None, nodetype=None, encoding='utf-8')
#print G.edges()

# 3a
## w/in met area, split into children/adults
### use same fraction of child/adults for each metro area
### determine fraction from entire US, fraction children out of only children & adult
### fraction children = alpha (a) = (#ch)/(#ch + #ad)
#### alt option - (find dif fraction for each metro area)

## data ##
# AGE SPLIT #
# POLYMOD contact study (Ref 22 in Apollini 2013)
## Children: (5 - 19) ?
## Adults: (20 - 59) ?
#popdata
#totalpop_age.csv

# import csv file of pop data
popdata = csv.reader(open('/home/anne/Dropbox/Anne_Bansal_lab/SDI_Data/totalpop_age.csv', 'r'),delimiter = ',')
zippopdata = csv.reader(open('/home/anne/Dropbox/Anne_Bansal_lab/SDI_Data/allpopstat_zip3_season_cl_nocodes.csv', 'r'),delimiter = ',')

# import data into dicts

d_pop_for_yr_age, ages, years = func.import_popdata(popdata, 0, 1, 2)

#group ages into children and adults
child = ['5-9 years', '10-14 years', '15-19 years']
adult = ['20-29 years', '30-39 years', '40-49 years', '50-59 years']
d_childpop, d_adultpop = func.pop_child_adult(d_pop_for_yr_age, years, child, adult)
#print d_childpop[2010]
#print d_adultpop[2010]

# set value of alpha from U.S. pop data
year = 2010 # decided to use pop from 2010 bc most recent
a = func.pop_child_frac(year, d_childpop, d_adultpop)
#print alpha # = 0.27 for 2010

#### calc alpha for each zip3 #####
# sensitivity analysis - look at variance across zip3s
season = 9
d_zip_popdata, zip_popdata = func.import_zip3_popdata(zippopdata, 0, 2, 3, 4)
#d_zip_popdata[(zip3, seas, age)] = pop
d_zip_alphas, alphas = func.zippop_child_frac(zip_popdata, season, d_zip_popdata)
#print len(alphas)

mean_alpha = np.mean(alphas)
std_alpha = np.std(alphas)
var_alpha = np.var(alphas)
#print var_alpha
#print std_alpha
#print mean_alpha

#print len(zip_popdata)

# 3b
## contact btwn children and adults
## 'C' = equation 3 in Apollini 2014
### alpha (a) = fraction of ch  --> calc in 3a
### (n) = ratio of ad/ch avg # contacts 
### (E) = avg fraction of contacts across age groups
#### --> from Euro data in Table 2 in Apollini 2013 + Ref 22

#Apollini 2013, Table 2, "Europe (average values)"
n = 0.79
E = 0.097

#q_1 = #see additional file of Apollini 2013
# Table 1 Mossong POLYMOD
# weighted avg for children and adults avg # of contacts

age = [5, 10, 15, 20, 30, 40, 50]
child = age[0:3]
adult = age[3:7]
contacts = [14.81, 18.22, 17.58, 13.57, 14.14, 13.83, 12.30]
participants = [661, 713, 685, 879, 815, 908, 906]

d_mean_contacts = dict(zip(age, contacts))
d_num_part = dict(zip(age, participants))

avg_q_ch = func.weighted_avg_q(child, d_mean_contacts, d_num_part)
avg_q_ad = func.weighted_avg_q(adult, d_mean_contacts, d_num_part)
#print avg_q_ad
#print d_num_part

#Apollini 2014 Eq 3
C_cc = ((a - E) / (a **2))
C_ca = (E / (a * (1 - a)))
C_ac = (E / (a * (1 - a)))
C_aa = (((n * (1 - a)) - E) / ((1 - a) **2))

# matrix
C_matrix = np.matrix([[C_cc, C_ca], [C_ac, C_aa]])
#print C_matrix
C = (C_matrix * avg_q_ch)
#print C

############################

## DISEASE ##

## Infc_list = [patient_Zero]
## Susc_list = [everyone_else]
## for t (time steps) in (1, 100):
    ## S --> I ?
    ## for s in susc:
        ## infect nodes with prob (1-e^-B(# of infec contacts)) # infected degree = #infected contacts
    
    ## I --> R ?
    ## for i in infc:
        ## recover with prob u
        ## (u = 0.1 - infectious period is 10 days)
    
