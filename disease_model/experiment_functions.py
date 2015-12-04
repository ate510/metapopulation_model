# Author: Anne Ewing
## Date: 11/16/15
## Function: script with experimental parameters and functions

import networkx as nx # to work with the population contact network
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
import operator

### local modules ###
sys.path.append('/home/anne/Dropbox/Anne_Bansal_lab')

### functions ###
import population_parameters as pop_func

#######################################
def list_exp_names ():
    
    time_intervals = ['cum_incd', 'real_time']
    exp_names = ['red_C_cc', 'red_all_child', 'red_beta']
    
    return time_intervals, exp_names

#######################################
def set_time_start_and_end (time_units, length_before, length_after):
    # length_before = how many weeks or time_steps before xmas
    # length_after = how many weeks or time_steps after xmas
    #wk 50 on cum_incid_graph 20%
    #wk 1 - 30%
    #import CSV file with wk# and cum incidence for each season
    ## take average of cumulative incidence across seasons
    #use time_step for xmas
    ## if epi trajectory is reasonable, just use time_step = 1 day
    ### can use actual day of xmas - take average across seasons
    
    # calc cum incid percent
    #cum_ili = np.cumsum(ili_list) # or use tot infc time series
    #sum_ili = sum(ili_list) # total numb infc at end
    #d_cum_percent = ((cum_ili)/(sum_ili)) * 100
    
    #if timing == 'cum_incd': #measure in weeks
    #    xmas = 51
    #    wk_start = xmas - length_before # (51 - 2 = 49)
    #    cum_incd = d_cum_incd[(wk_start)] # 49
    #    time_step = [time for time in d_cum_incd_model[(time)] if d_cum_incd_model[(time)] == cum_incd]
    #    #time_start = #time_step when cum incidence = 20%
    #    #time_end = 
    
    time_start = 2 #80
    time_end = 3 #94
    
    return time_start, time_end
    
#######################################
#def set_time_end (time_units, ):
#
#    
#######################################
def reduce_C_cc (C):
    # reduce child to child contacts only
    
    C_cc = C.item((0, 0))
    C_ca = C.item((0, 1))
    C_ac = C.item((1, 0))
    C_aa = C.item((1, 1))
    red_percent = .4
    red_C_cc = (C_cc - (C_cc * red_percent))
    C_exp = np.matrix([[red_C_cc, C_ca], [C_ac, C_aa]])
    #C.item((0,0)) = red_C_cc #reassign C_cc in contact matrix
    
    return C_exp
    
#######################################
#def reduce_all_C_child ():
    # reduce all child contacts