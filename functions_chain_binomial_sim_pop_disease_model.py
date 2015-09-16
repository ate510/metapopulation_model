## Author: Anne Ewing
## Date: 06/25/15
## Function: Group of functions for project model

### packages/modules ###
import csv
import sys
import datetime as date
import networkx as nx


# functions #

#create function to import pop of each metro area
##Names in script: d_pop_for_metro
##args: (metropop, 0, 1)

def import_metropop (datafile, metrocol, popcol):
    
    dict_metropop, metro_list = {}, []
    for row in datafile:
        metro_id = int(row[metrocol])
        metro_list.append(metro_id)
        pop = float(row[popcol])
        dict_metropop[(metro_id)] = pop
    return dict_metropop, list(set(metro_list))

#create a function that imports pop data into dictionaries and lists
##Names in script: d_yr_age, ages, years = {}, [], []
##args: (popdata, 0, 1, 2)

def import_popdata (datafile, yrcol, agecol, popcol):
    
    dict_popdata, age_list, yr_list = {}, [], []
    for row in datafile:
        year = int(row[yrcol])
	yr_list.append(year)
        age = str.lower(row[agecol])
	age_list.append(age)
        pop = float(row[popcol])
        dict_popdata[(year, age)] = pop
    return dict_popdata, list(set(age_list)), list(set(yr_list))
    
#create a function to import pop data at zip3 level  (v032215 incidence_funcitons_RRdata)
##Name in script: d_zip_popdata
##args: (popdata, 0, 2, 3, 4)
# in script: zip3_top_bottom_20_percentchange_variance_notshifted_v062315.py
def import_zip3_popdata (datafile, zip3col, popcol, seascol, agecol):
# if we come back to this function (not in population_parameters), 
    # double check what ages the groups A and C include

	dict_popdata, zips_popdata = {}, []
	for row in datafile:
            zip3 = int(row[zip3col])
            zips_popdata.append(zip3)
            pop = float(row[popcol])
            seas = int(row[seascol])
            age = str(row[agecol])
            dict_popdata[(zip3, seas, age)] = pop
	#dict_popdata_tot[(zip3, seas)] = ((dict_popdata[zip3, seas, 'A']) + (dict_popdata[zip3, seas, 'C']) + (dict_popdata[zip3, seas, 'O']))
	return dict_popdata, list(set(zips_popdata))
    
# create a function to combine children and adults into 2 groups
##Names in script: d_childpop, d_adultpop = {}, {}
##args: (d_pop_for_yr_age, years, child, adult)

def pop_child_adult (dict_popdata, years, child_list, adult_list):
    
    dict_childpop, dict_adultpop = {}, {}
    for y in years:
        childpop = sum([dict_popdata[y, a] for a in child_list])
        #childpop = ((d_yr_age[y, '5-9 years']) + (d_yr_age[y, '10-14 years']) + (d_yr_age[y, '15-19 years']))
        dict_childpop[y] = childpop
        adultpop = sum([dict_popdata[y, a] for a in adult_list])
        #adultpop = ((d_yr_age[y, '20-29 years']) + (d_yr_age[y, '30-39 years']) + (d_yr_age[y, '40-49 years']) + (d_yr_age[y, '50-59 years']))
        dict_adultpop[y] = adultpop
        
    return dict_childpop, dict_adultpop

# determine alpha = the fraction of the U.S. population that is children 
##Names in script: a = value
##args: (year, d_childpop, d_adultpop)

def pop_child_frac (year, dict_childpop, dict_adultpop):
    
    childpop = dict_childpop[year]
    adultpop = dict_adultpop[year]
    alpha = ((childpop) / (childpop + adultpop))
    
    return alpha
    
# determine alpha = the fraction of the U.S. population that is children for each zip3
##Names in script: d_zip_alphas
##args: (zip_popdata, season, d_zip_popdata)

def zippop_child_frac (zip3_list, season, dict_zippop):
    
    dict_zip_alphas, alpha_list = {}, []
    for zip3 in zip3_list:
        childpop = dict_zippop[zip3, season, 'C']
        adultpop = dict_zippop[zip3, season, 'A']
        alpha = ((childpop) / (childpop + adultpop))
        dict_zip_alphas[zip3] = alpha
        alpha_list.append(alpha)
        
    return dict_zip_alphas, alpha_list
    
def weighted_avg_q (child, d_mean_contacts, d_num_part):

    qxN = [(d_mean_contacts[age] * d_num_part[age]) for age in child]
    N_list = [d_num_part[age] for age in child]
    avg_q = (sum(qxN))/(sum(N_list))
    
    return avg_q

########################################
#import edgelist data without using networkx
def read_edgelist_anne (filename):
    
    G = nx.Graph()
    file = open(filename, 'rU')
    reader = csv.reader(file)
    for row in reader:
        data = row[0].split('\t')
        G.add_edge(int(data[0]),int(data[1]), weight=float(data[2]))
        
    return G
    
