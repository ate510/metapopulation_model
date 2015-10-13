# plot histograms of percent epidemic
# x axis: number of sims
# y axis: number of subpops above epidemic threshold x%


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
def import_sim_data (filename, simcol, timecol, metrocol, agecol, totcol):
# import data from csv file storing all sim results
# (filename, 0, 1, 2, 3, 5)
    
    #import csv file
    datafile = csv.reader(open(filename, 'r'),delimiter = ',')
    headers = datafile.next() #omit headers
    
    sim_list, time_list, metro_list = [], [], []
    dict_age_infc, dict_total_infc = {}, {}
    for row in datafile:
        #year = int(row[yrcol])
        #yr_list.append(year)
        sim = int(row[simcol])
        sim_list.append(sim)
        time_step = int(row[timecol])
        time_list.append(time_step)
        metro_id = int(row[metrocol])
        metro_list.append(metro_id)
        age = str(row[agecol])
        #currently_infected = float(row[currentcol])
        total_infected = float(row[totcol])
        dict_age_infc[(sim, time_step, metro_id, age)] = total_infected
        #dict_total_infc[(sim, time_step, metro_id)] = ((dict_age_infc[(sim, time_step, metro_id, 'A')]) + (dict_age_infc[(sim, time_step, metro_id, 'C')]))
        
    return dict_age_infc, list(set(sim_list)), list(set(time_list)), list(set(metro_list))
    
###################################################
def infc_final_size (time_list, dict_age_infc):
    
    final_time_step = (time_list[-1])
    dict_final_infc = {}
    for (sim, time_step, metro_id, age) in dict_age_infc:
    #for time_step in time_list:
        if time_step == final_time_step:
            dict_final_infc[(sim, metro_id, age)] = dict_age_infc[(sim, final_time_step, metro_id, age)]
            
    return dict_final_infc

###################################################
def infc_agg_ages (dict_final_infc):
    
    dict_total_infc = {}
    for (sim, metro_id, age) in dict_final_infc:
        dict_total_infc[(sim, metro_id)] = ((dict_final_infc[(sim, metro_id, 'A')]) + (dict_final_infc[(sim, metro_id, 'C')]))
    
    return dict_total_infc
    
###################################################
def divide_by_pop (d_metropop, metro_ids, dict_total_infc):
    
    dict_percent_epidemic = {}
    #population_size = sum([d_metropop[x] for x in metro_ids]) #calculate total population over all metro areas, not just those infected
    for (sim, metro_id) in dict_total_infc:
        dict_percent_epidemic[(sim, metro_id)] = ((dict_total_infc[(sim, metro_id)])/(d_metropop[metro_id])) #divide each metro area # infected by # people - get fraction for each metro 
    return dict_percent_epidemic
          
###################################################
def histogram (threshold, sim_list, metro_list, dict_percent_epidemic):
# metro_ids = metro_list

    # number of metro areas above x% for given sim
    
    dict_plot_metros, dict_number_metros = {}, {} #create dictionary with sim number as key and list of metros over the threshold as value, and dictionary with number of metros as value
    for sim in sim_list: # loop thru sims
        metros_over_threshold = [] # for each sim, create empty list of metros 
        for metro in metro_list: # loop through all metros
            epi_size = dict_percent_epidemic[(sim, metro)] # define epi size as the fraction of the metro population infected at last time step
            if epi_size > threshold: # if the epi size is greater than the set threshold
                metros_over_threshold.append(metro) # add metro where epi size exceeds threshold to list
        dict_plot_metros[(sim)] = metros_over_threshold # after looping through all metros, assign list to dictionary with sim as key
        dict_number_metros[(sim)] = (len(metros_over_threshold)) # for histogram, don't care about which metros are over, just how many
        
    xvalues = [dict_number_metros[(sim)] for sim in sim_list]
    #xvalues = [dict_cum_incid_perc[(wk, zip3)] for zip3 in zip_list]
    #print xvalues
    
    percent_threshold = (threshold * 100)
    
    num_bins = 100

    ## the histogram of the data
    n, bins, patches = plt.hist(xvalues, num_bins, facecolor='green', alpha=0.5)
    plt.xlabel('Number of Metros with more than %s%% of Population Infected' % (percent_threshold))
    #plt.xlabel('Cum Incid %')
    plt.ylabel('Number of Sims')
    #plt.ylabel('Proportion of Zip Threes')
    plt.title(r'Histogram Threshold: %s%% Epidemic' % (percent_threshold))
    #plt.title(r'Histogram of Cum Incid: %s' % (week))
    #plt.show()
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Threshold_Histograms/histogram_num_metros_above_threshold_%s.png' % (percent_threshold))
    plt.close()

          

###################################################
if __name__ == "__main__":

    # READ CSV FILE    
    sim_results = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Saved_CSV_Outputs/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_numsims_100.csv'
    
    #Import population info
    filename_metropop = 'Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/metedges.txt'
    d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
    
    # THRESHOLD
    threshold_list = [0.03, .05, .10, .15, .2, .5, .7, .9, 1]
    
    d_age_infc, sim_list, time_list, metro_list = import_sim_data(sim_results, 0, 1, 2, 3, 5)
    d_final_infc = infc_final_size(time_list, d_age_infc)
    #print d_final_infc
    d_tot_infc = infc_agg_ages(d_final_infc)
    #print d_tot_infc
    d_perc_epi = divide_by_pop(d_metropop, metro_ids, d_tot_infc)
    #print d_perc_epi
    
    #print len(metro_ids)
    #print len(metro_list) # length of both = 225
    
    # HISTOGRAM
    for t in threshold_list:
        histogram(t, sim_list, metro_list, d_perc_epi)