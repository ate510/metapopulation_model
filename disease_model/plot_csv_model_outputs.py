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
def import_csv_output (filename, metzerocol, timecol, metidcol, agecol, currentcol, totcol):
#create function to import ILI time series from model output
    
    # import US metro pop data
    datafile = csv.reader(open(filename, 'r'),delimiter = ',')
    headers = datafile.next()
        
    dict_currentflu, dict_currentflu_adult, dict_currentflu_child, dict_totflu_adult, dict_totflu_child, metro_list = {}, {}, {}, {}, {}, []
    for row in datafile:
        metro_zero = int(row[metzerocol])
        time_step = int(row[timecol])
        metro_id = int(row[metidcol])
        metro_list.append(metro_id)
        age = str(row[agecol]) #'A' or 'C'
        current_flu = float(row[currentcol])
        tot_flu = float(row[totcol])
        dict_currentflu[(metro_zero, time_step, metro_id, age)] = current_flu
        if age == 'A':
            dict_currentflu_adult[(metro_zero, time_step, metro_id)] = current_flu
            dict_totflu_adult[(metro_zero, time_step, metro_id)] = tot_flu
        elif age == 'C':
            dict_currentflu_child[(metro_zero, time_step, metro_id)] = current_flu
            dict_totflu_child[(metro_zero, time_step, metro_id)] = tot_flu
        
    return dict_currentflu_adult, dict_currentflu_child, dict_totflu_adult, dict_totflu_child, list(set(metro_list))

###################################################
def sum_national_curve (dict_currentflu_adult, dict_currentflu_child, metro_ids, intervention):
# this function will sum flu at each time step to find the national peak    
    
    time_series_child = list(set([key[1] for key in sorted(dict_currentflu_child)]))
    time_series_adult = list(set([key[1] for key in sorted(dict_currentflu_adult)]))
    time_series = list(set(time_series_child + time_series_adult))
    infc_met_ids_child = list(set([key[2] for key in sorted(dict_currentflu_child)]))
    infc_met_ids_adult = list(set([key[2] for key in sorted(dict_currentflu_adult)]))
    infc_met_ids = list(set(infc_met_ids_child + infc_met_ids_adult))
    
    #try to create a list of all possible tuples in graph
    list_tuples_graph = []
    for t in time_series:
        for met_id in infc_met_ids:
            list_tuples_graph.append((1, t, met_id))
    
    #assign missing keys zero
    for tpl in list_tuples_graph:
        if tpl not in dict_currentflu_child:
            dict_currentflu_child[tpl] = dict_currentflu_child.get(tpl, 0)
        if tpl not in dict_currentflu_adult:
            dict_currentflu_adult[tpl] = dict_currentflu_adult.get(tpl, 0)
    
    national_incidence_time_series, child_incidence_time_series, adult_incidence_time_series = {}, {}, {} #create dictionaries for national time series (all ages, child, and adult)
    #for (1, t, met_id) in dict_currentflu_adult:
    #    if (1, t, met_id) not in dict_currentflu_child:
    #        dict_currentflu_child[(1, t, met_id)] = dict_currentflu_child.get((1, t, met_id), 0)
    #for (1, t, met_id) in dict_currentflu_adult:
    #    if (1, t, met_id) not in dict_currentflu_child:
    #        dict_currentflu_child[(1, t, met_id)] = dict_currentflu_child.get((1, t, met_id), 0)
    for t in time_series: 
        sum_metro_child = sum([dict_currentflu_child[(1, t, met_id)] for met_id in infc_met_ids_child])
        sum_metro_adult = sum([dict_currentflu_adult[(1, t, met_id)] for met_id in infc_met_ids_adult])
        child_incidence_time_series[(1, t, intervention)] = sum_metro_child
        adult_incidence_time_series[(1, t, intervention)] = sum_metro_adult
        sum_ages = sum_metro_child + sum_metro_adult
        national_incidence_time_series[(1, t, intervention)] = sum_ages
    
    return national_incidence_time_series, child_incidence_time_series, adult_incidence_time_series, list_tuples_graph, time_series
      
###################################################
def plot_national_cases (national_incidence_time_series, child_incidence_time_series, adult_incidence_time_series, interventions, time_series):
# time series for currently infected cases
# one plot for child one for adult
# one line for each metro area

    lncolor = ['hotpink', 'red', 'orange', 'gold', 'green', 'blue'] #,'cyan', 'darkviolet']
    
 #national - all ages
    for interv, c in zip(interventions, lncolor):
        #time_series = range(0, time_end)
        #time_series = list(set([key[2] for key in sorted(d_metro_infected_child)]))
        graphxax = time_series
        graphyax = [national_incidence_time_series[(1, t, interv)] for t in time_series]
        plt.plot(graphxax, graphyax, color = c, label = interv)
        
    #plt.xlim([0, 9])
    #plt.ylim([0, 1000000])
    plt.xlabel('Time Step')
    plt.ylabel('National Current Cases')
    plt.legend(loc = 'upper left', fontsize = 'small')
    #plt.xticks(range(0, 10), wklab)
    plt.show()
    #plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots_Jan_2016/Control/Current_Cases/Child/chain_binomial_current_cases_child_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (alpha, ch_travelers_r, R0, gamma, beta, metro_zero)) #use this line to save fig in jan 2016 folder for base params
    #plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots_Jan_2016/Experiments/%s/Current_Cases/Child/chain_binomial_current_cases_child_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (intervention, alpha, ch_travelers_r, R0, gamma, beta, metro_zero)) # use this line to save fig in jan 2016 folder for experimental params
    #plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Current_Cases/Child/chain_binomial_current_cases_child_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
    #plt.close()

##adult
#    for met_id in metro_ids:
#        #time_series = range(0, time_end)
#        time_series = list(set([key[2] for key in sorted(d_metro_infected_adult)]))
#        graphxax = time_series
#        graphyax = [d_metro_infected_adult[(1, t, )] for t in time_series]
#        plt.plot(graphxax, graphyax)
#        
#    #plt.xlim([0, 40])
#    #plt.ylim([0, 1000])
#    plt.xlabel('Time Step')
#    plt.ylabel('Current Cases - Adult')
#    #plt.xticks(range(0, 10), wklab)
#    #plt.show()
#    #plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots_Jan_2016/Control/Current_Cases/Adult/chain_binomial_current_cases_adult_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (alpha, ch_travelers_r, R0, gamma, beta, metro_zero)) #use this line to save fig in jan 2016 folder for base params
#    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots_Jan_2016/Experiments/%s/Current_Cases/Adult/chain_binomial_current_cases_adult_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (intervention, alpha, ch_travelers_r, R0, gamma, beta, metro_zero)) # use this line to save fig in jan 2016 folder for experimental params
##    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Current_Cases/Adult/chain_binomial_current_cases_adult_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
#    plt.close()

###################################################
if __name__ == "__main__":
    
    # IMPORT MODEL RESULTS
    #filename_red_C_cc = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/1. Red_C_cc by 58%/Thru time_step 500/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_cc_travel_none.csv'
    filename_baseline = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_none_travel_none.csv'
    filename_red_C_cc = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_cc_travel_none.csv'
    filename_red_C_aa = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_aa_travel_none.csv'
    filename_red_C_all = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_red_C_all_travel_none.csv'
    filename_inc_child_trav = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_none_travel_inc_child_trav.csv'
    filename_inc_all_trav = 'Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Exp_Results_Feb_2016/chain_binomial_output_nummetrozeros_1_numchildzeros_1_numadultzeros_0_disease_none_travel_inc_all_trav.csv'
    d_current_adult, d_current_child, d_tot_adult, d_tot_child, met_ids = import_csv_output(filename_baseline, 0, 1, 2, 3, 4, 5)    
    d_current_adult_cc, d_current_child_cc, d_tot_adult_cc, d_tot_child_cc, met_ids_1 = import_csv_output(filename_red_C_cc, 0, 1, 2, 3, 4, 5)
    d_current_adult_aa, d_current_child_aa, d_tot_adult_aa, d_tot_child_aa, met_ids_2 = import_csv_output(filename_red_C_aa, 0, 1, 2, 3, 4, 5)
    d_current_adult_C_all, d_current_child_C_all, d_tot_adult_C_all, d_tot_child_C_all, met_ids_3 = import_csv_output(filename_red_C_all, 0, 1, 2, 3, 4, 5)
    d_current_adult_child_trav, d_current_child_child_trav, d_tot_adult_child_trav, d_tot_child_child_trav, met_ids_4 = import_csv_output(filename_inc_child_trav, 0, 1, 2, 3, 4, 5)
    d_current_adult_all_trav, d_current_child_all_trav, d_tot_adult_all_trav, d_tot_child_all_trav, met_ids_5 = import_csv_output(filename_inc_all_trav, 0, 1, 2, 3, 4, 5)



        
    #interventions = ['baseline', 'red_C_cc', 'red_C_aa', 'red_C_all', 'inc_child_trav', 'inc_all_trav']
    interventions = ['baseline']
    
    
    # SUM DATA TO NATIONAL LEVEL
    national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult, d_current_child, met_ids, (interventions[0]))
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_cc, d_current_child_cc, met_ids_1, (interventions[1]))
    ##national_time_series.update(
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_aa, d_current_child_aa, met_ids_2, (interventions[2]))
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_C_all, d_current_child_C_all, met_ids_3, (interventions[3]))
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_child_trav, d_current_child_child_trav, met_ids_4, (interventions[4]))
    #national_time_series, child_time_series, adult_time_series, list_tuples, time_steps = sum_national_curve(d_current_adult_all_trav, d_current_child_all_trav, met_ids_5, (interventions[5]))    
    #
    
    # PLOT
    plot_national_cases(national_time_series, child_time_series, adult_time_series, interventions, time_steps)
