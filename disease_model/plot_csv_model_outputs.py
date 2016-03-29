# Author: Anne Ewing
## Date: 2/25/16
## Function: script to plot csv outputs from model results

import networkx as nx # to work with the population contact network
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
import operator

### local modules ###
sys.path.append('/home/anne/Dropbox/Anne_Bansal_lab')

###################################################
def import_csv_output (filename, metzerocol, timecol, metidcol, agecol, currentcol, totcol, intervention, dict_currentflu, dict_currentflu_adult, dict_currentflu_child, dict_totflu_adult, dict_totflu_child):
#create function to import ILI time series from model output
    
    # import US metro pop data
    datafile = csv.reader(open(filename, 'r'),delimiter = ',')
    headers = datafile.next()
        
    for row in datafile:
        metro_zero = int(row[metzerocol])
        time_step = int(row[timecol])
        metro_id = int(row[metidcol])
        #metro_list.append(metro_id)
        age = str(row[agecol]) #'A' or 'C'
        current_flu = float(row[currentcol])
        tot_flu = float(row[totcol])
        dict_currentflu[(intervention, metro_zero, time_step, metro_id, age)] = current_flu
        if age == 'A':
            dict_currentflu_adult[(intervention, metro_zero, time_step, metro_id)] = current_flu
            dict_totflu_adult[(intervention, metro_zero, time_step, metro_id)] = tot_flu
        elif age == 'C':
            dict_currentflu_child[(intervention, metro_zero, time_step, metro_id)] = current_flu
            dict_totflu_child[(intervention, metro_zero, time_step, metro_id)] = tot_flu
        
    #return dict_currentflu, dict_currentflu_adult, dict_currentflu_child, dict_totflu_adult, dict_totflu_child #, list(set(metro_list))

###################################################
def sum_national_curve (dict_currentflu_adult, dict_currentflu_child, intervention, metro_zero, national_incidence_time_series, child_incidence_time_series, adult_incidence_time_series):
# this function will sum flu at each time step to find the national peak    
    
    # aggregate time steps and metro ids across interventions
    time_series_child = list(set([key[2] for key in sorted(dict_currentflu_child)]))
    time_series_adult = list(set([key[2] for key in sorted(dict_currentflu_adult)]))
    time_series = list(set(time_series_child + time_series_adult))
    infc_met_ids_child = list(set([key[3] for key in sorted(dict_currentflu_child)]))
    infc_met_ids_adult = list(set([key[3] for key in sorted(dict_currentflu_adult)]))
    infc_met_ids = list(set(infc_met_ids_child + infc_met_ids_adult))
    
    #try to create a list of all possible tuples in the graph
    list_tuples_graph = []
    for t in time_series:
        for met_id in infc_met_ids:
            list_tuples_graph.append((intervention, metro_zero, t, met_id))
    
    #assign missing keys zero
    for tpl in list_tuples_graph:
        if tpl not in dict_currentflu_child:
            dict_currentflu_child[tpl] = dict_currentflu_child.get(tpl, 0)
        if tpl not in dict_currentflu_adult:
            dict_currentflu_adult[tpl] = dict_currentflu_adult.get(tpl, 0)
    
    #national_incidence_time_series, child_incidence_time_series, adult_incidence_time_series = {}, {}, {} #create dictionaries for national time series (all ages, child, and adult)
    #for (1, t, met_id) in dict_currentflu_adult:
    #    if (1, t, met_id) not in dict_currentflu_child:
    #        dict_currentflu_child[(1, t, met_id)] = dict_currentflu_child.get((1, t, met_id), 0)
    #for (1, t, met_id) in dict_currentflu_adult:
    #    if (1, t, met_id) not in dict_currentflu_child:
    #        dict_currentflu_child[(1, t, met_id)] = dict_currentflu_child.get((1, t, met_id), 0)
    for t in time_series: 
        sum_metro_child = sum([dict_currentflu_child[(intervention, metro_zero, t, met_id)] for met_id in infc_met_ids_child])
        sum_metro_adult = sum([dict_currentflu_adult[(intervention, metro_zero, t, met_id)] for met_id in infc_met_ids_adult])
        child_incidence_time_series[(intervention, metro_zero, t)] = sum_metro_child
        adult_incidence_time_series[(intervention, metro_zero, t)] = sum_metro_adult
        sum_ages = sum_metro_child + sum_metro_adult
        national_incidence_time_series[(intervention, metro_zero, t)] = sum_ages
    
    return time_series
    #return national_incidence_time_series, child_incidence_time_series, adult_incidence_time_series, list_tuples_graph, time_series
      
###################################################
def plot_national_cases (national_incidence_time_series, child_incidence_time_series, adult_incidence_time_series, interventions, time_series, beta, metro_zero, timing):
# time series for currently infected cases
# one plot for child one for adult
# one line for each metro area

    lncolor = ['hotpink', 'red', 'orange', 'green', 'blue', 'darkviolet'] #,'cyan','gold']
    
    #national - all ages
    for interv, c in zip(interventions, lncolor):
        #time_series = range(0, time_end)
        #time_series = list(set([key[2] for key in sorted(d_metro_infected_child)]))
        graphxax = time_series
        graphyax = [national_incidence_time_series[(interv, metro_zero, t)] for t in time_series]
        plt.plot(graphxax, graphyax, color = c, label = interv)
        
    plt.plot((134, 134), (0, 500000), color = 'grey')
    plt.plot((141, 141), (0, 500000), color = 'black')
    plt.plot((148, 148), (0, 500000), color = 'grey')

        
    #plt.xlim([100, 160])
   # plt.ylim([0, 500000])
    plt.xlabel('Time Step')
    plt.ylabel('Current Cases')
    plt.title('National')
    plt.legend(loc = 'upper left', fontsize = 'small')
    #plt.xticks(range(0, 10), wklab)
    #plt.show()
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/national_time_series_compare_interventions_beta_%s_metrozero_%s_timing_%s.png' % (beta, metro_zero, timing))
    plt.close()

 #national - children
    for interv, c in zip(interventions, lncolor):
        graphxax = time_series
        graphyax = [child_incidence_time_series[(interv, metro_zero, t)] for t in time_series]
        plt.plot(graphxax, graphyax, color = c, label = interv)

    #plt.plot((134, 134), (0, 500000), color = 'grey')
    #plt.plot((141, 141), (0, 500000), color = 'black')
    #plt.plot((148, 148), (0, 500000), color = 'grey')
                
    #plt.xlim([100, 160])
    #plt.ylim([0, 500000])
    plt.xlabel('Time Step')
    plt.ylabel('Current Cases')
    plt.title('National - Child')
    plt.legend(loc = 'upper left', fontsize = 'small')
    #plt.xticks(range(0, 10), wklab)
    #plt.show()
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/child_national_time_series_compare_interventions_beta_%s_metrozero_%s_timing_%s.png' % (beta, metro_zero, timing))
    plt.close()
    
    #national - adult
    for interv, c in zip(interventions, lncolor):
        graphxax = time_series
        graphyax = [adult_incidence_time_series[(interv, metro_zero, t)] for t in time_series]
        plt.plot(graphxax, graphyax, color = c, label = interv)

    #plt.plot((134, 134), (0, 500000), color = 'grey')
    #plt.plot((141, 141), (0, 500000), color = 'black')
    #plt.plot((148, 148), (0, 500000), color = 'grey')
    #                
    #plt.xlim([100, 160])
    #plt.ylim([0, 500000])
    plt.xlabel('Time Step')
    plt.ylabel('Current Cases')
    plt.title('National - Adult')
    plt.legend(loc = 'upper left', fontsize = 'small')
    #plt.xticks(range(0, 10), wklab)
    plt.show()
    #plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/adult_national_time_series_compare_interventions_beta_%s_metrozero_%s_timing_%s.png' % (beta, metro_zero, timing))
    #plt.close()

###################################################
if __name__ == "__main__":
    
    # CREATE DICTIONARIES TO IMPORT DATA INTO
    d_current, d_current_adult, d_current_child, d_tot_adult, d_tot_child = {}, {}, {}, {}, {} # need dicts outside of function becuase importing all data from separate files into these dicts
    # lenght of d_current should be 892350 with only baseline data

    interventions = ['baseline', 'red_C_cc', 'red_C_aa', 'red_C_all', 'inc_child_trav', 'inc_all_trav', 'combo']
    #interventions = ['baseline', 'red_C_cc']
    beta = 0.03 # from simulation, for notation purposes only in figure title
    metro_zero = 1
    #timing = 'extreme'
    timing = 'actual'
    
    # IMPORT MODEL RESULTS
    filename_baseline = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_none_travel_none.csv'
    filename_red_C_cc = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_cc_travel_none.csv'
    filename_red_C_aa = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_aa_travel_none.csv'
    filename_red_C_all = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_all_travel_none.csv'
    filename_inc_child_trav = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_none_travel_inc_child_trav.csv'
    filename_inc_all_trav = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_none_travel_inc_all_trav.csv'
    filename_combo = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_all_travel_inc_all_trav_timing_actual.csv' #red_C_all and inc_all_trav
    
    #combo interventions
    #filename_red_C_cc_inc_child_trav = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_cc_travel_inc_child_trav.csv'

    
    #extreme timing
    #filename_baseline = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_none_travel_none_timing_extreme.csv'
    #filename_red_C_cc = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_cc_travel_none_timing_extreme.csv'
    #filename_red_C_aa = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_aa_travel_none_timing_extreme.csv'
    #filename_red_C_all = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_all_travel_none_timing_extreme.csv'
    #filename_inc_child_trav = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_none_travel_inc_child_trav_timing_extreme.csv'
    #filename_inc_all_trav = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_none_travel_inc_all_trav_timing_extreme.csv'
    #
    # UPDATE DICTIONARIES WITH DATA
    import_csv_output(filename_baseline, 0, 1, 2, 3, 4, 5, interventions[0], d_current, d_current_adult, d_current_child, d_tot_adult, d_tot_child)
    import_csv_output(filename_red_C_cc, 0, 1, 2, 3, 4, 5, interventions[1], d_current, d_current_adult, d_current_child, d_tot_adult, d_tot_child)
    import_csv_output(filename_red_C_aa, 0, 1, 2, 3, 4, 5, interventions[2], d_current, d_current_adult, d_current_child, d_tot_adult, d_tot_child)
    import_csv_output(filename_red_C_all, 0, 1, 2, 3, 4, 5, interventions[3], d_current, d_current_adult, d_current_child, d_tot_adult, d_tot_child)
    import_csv_output(filename_inc_child_trav, 0, 1, 2, 3, 4, 5, interventions[4], d_current, d_current_adult, d_current_child, d_tot_adult, d_tot_child)
    import_csv_output(filename_inc_all_trav, 0, 1, 2, 3, 4, 5, interventions[5], d_current, d_current_adult, d_current_child, d_tot_adult, d_tot_child)
    import_csv_output(filename_combo, 0, 1, 2, 3, 4, 5, interventions[6], d_current, d_current_adult, d_current_child, d_tot_adult, d_tot_child)
    
    #d_current_adult, d_current_child, d_tot_adult, d_tot_child, met_ids = import_csv_output(filename_baseline, 0, 1, 2, 3, 4, 5)    
    #d_current_adult_cc, d_current_child_cc, d_tot_adult_cc, d_tot_child_cc, met_ids_1 = import_csv_output(filename_red_C_cc, 0, 1, 2, 3, 4, 5)
    #d_current_adult_aa, d_current_child_aa, d_tot_adult_aa, d_tot_child_aa, met_ids_2 = import_csv_output(filename_red_C_aa, 0, 1, 2, 3, 4, 5)
    #d_current_adult_C_all, d_current_child_C_all, d_tot_adult_C_all, d_tot_child_C_all, met_ids_3 = import_csv_output(filename_red_C_all, 0, 1, 2, 3, 4, 5)
    #d_current_adult_child_trav, d_current_child_child_trav, d_tot_adult_child_trav, d_tot_child_child_trav, met_ids_4 = import_csv_output(filename_inc_child_trav, 0, 1, 2, 3, 4, 5)
    #d_current_adult_all_trav, d_current_child_all_trav, d_tot_adult_all_trav, d_tot_child_all_trav, met_ids_5 = import_csv_output(filename_inc_all_trav, 0, 1, 2, 3, 4, 5)

    
    # SUM DATA TO NATIONAL LEVEL
    national_time_series, child_time_series, adult_time_series = {}, {}, {} # create dictionaries for national time series (all ages, child, and adult)
    time_steps = sum_national_curve(d_current_adult, d_current_child, interventions[0], metro_zero, national_time_series, child_time_series, adult_time_series)
    _ = sum_national_curve(d_current_adult, d_current_child, interventions[1], metro_zero, national_time_series, child_time_series, adult_time_series)
    _ = sum_national_curve(d_current_adult, d_current_child, interventions[2], metro_zero, national_time_series, child_time_series, adult_time_series)
    _ = sum_national_curve(d_current_adult, d_current_child, interventions[3], metro_zero, national_time_series, child_time_series, adult_time_series)
    _ = sum_national_curve(d_current_adult, d_current_child, interventions[4], metro_zero, national_time_series, child_time_series, adult_time_series)
    _ = sum_national_curve(d_current_adult, d_current_child, interventions[5], metro_zero, national_time_series, child_time_series, adult_time_series)
    _ = sum_national_curve(d_current_adult, d_current_child, interventions[6], metro_zero, national_time_series, child_time_series, adult_time_series)

    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult, d_current_child, (interventions[0]))
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_cc, d_current_child_cc, met_ids_1, (interventions[1]))
    ##national_time_series.update(
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_aa, d_current_child_aa, met_ids_2, (interventions[2]))
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_C_all, d_current_child_C_all, met_ids_3, (interventions[3]))
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_child_trav, d_current_child_child_trav, met_ids_4, (interventions[4]))
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_all_trav, d_current_child_all_trav, met_ids_5, (interventions[5]))    
    #
    
    # PLOT
    plot_national_cases(national_time_series, child_time_series, adult_time_series, interventions, time_steps, beta, metro_zero, timing)
