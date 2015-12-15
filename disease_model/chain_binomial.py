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
import experiment_functions as exp_func
    

###################################################    
def infected_contact_child(C, Infc_A, Infc_C, metro_id, d_metro_age_pop):
# returns infected degree based on contact probabilities for adults

    C_cc = C.item((0, 0)) #((row, column))
    C_ac = C.item((1, 0))
    child_pop = d_metro_age_pop[(metro_id, 'child')]
    adult_pop = d_metro_age_pop[(metro_id, 'adult')]
    
    return ((C_cc*(Infc_C/child_pop))+(C_ac*(Infc_A/adult_pop)))
        
###################################################    
def infected_contact_adult(C, Infc_A, Infc_C, metro_id, d_metro_age_pop):
# returns infected degree based on contact probabilities for children


    C_ca = C.item((0, 1))
    C_aa = C.item((1, 1))
    child_pop = d_metro_age_pop[(metro_id, 'child')]
    adult_pop = d_metro_age_pop[(metro_id, 'adult')]

    return ((C_ca*(Infc_C/child_pop))+(C_aa*(Infc_A/adult_pop)))
    #call this function in a larger one with population sizes

###################################################
def lambda_child_calc(C, Infc_A, Infc_C, metro_id, d_metro_age_pop, beta):
    
    infc_cont_child = infected_contact_child(C, Infc_A, Infc_C, metro_id, d_metro_age_pop)
    
    lambda_child = ((beta) * (infc_cont_child))
    
    return lambda_child

###################################################    
def lambda_adult_calc(C, Infc_A, Infc_C, metro_id, d_metro_age_pop, beta):
# calculates force of infection (lambda)

    infc_cont_adult = infected_contact_adult(C, Infc_A, Infc_C, metro_id, d_metro_age_pop)

    lambda_adult = ((beta) * (infc_cont_adult))
    
    return lambda_adult
    
####################################################
#def susc_infc_binomial_prob_child (met_id, d_Infc, C, beta, d_metro_age_pop):
## returns probability of infection for binomial function    
#    # Si --> Ii
#    # determine how many susceptibles get infected in each metro area
#  
#    Infc_C = d_Infc[(met_id, 'C')] #number children infected in that metro area
#    Infc_A = d_Infc[(met_id, 'A')] #number adults infected in that metro area
#    infected_degree_of_s_C = infected_contact_child(C, Infc_A, Infc_C, met_id, d_metro_age_pop) #number of infected contacts
#    
#    return (susc_infc_event_child(infected_degree_of_s_C, beta))
#    
#    #if susc_infc_event_child(infected_degree_of_s_C, beta):
#    #    new_cases = new_cases + 1
#        
####################################################
#def susc_infc_binomial_prob_adult (met_id, d_Infc, C, beta):
## returns probability of infection for binomial function    
#
#    # Si --> Ii
#    # determine how many susceptibles get infected in each metro area
#
#    Infc_C = d_Infc[(met_id, 'C')]
#    Infc_A = d_Infc[(met_id, 'A')]  
#    infected_degree_of_s_A = infected_contact_adult(C, Infc_A, Infc_C)
#    
#    return (susc_infc_event_adult(infected_degree_of_s_A, beta))
 
              
####################################################    
#def infc_recv_event(gamma):
## returns true or false for recovery event with probability (gamma)
##   assumes exponential distribution for infectious period
## same for children and adults
#
#    return (rnd.random() < gamma)

###################################################
def SIR_initial_pops(metro_ids, d_metro_age_pop):
        
    d_Susc, d_Infc, d_Recv = {}, {}, {}
    for met_id in metro_ids:
        child_pop = d_metro_age_pop[(met_id, 'child')] #
        adult_pop = d_metro_age_pop[(met_id, 'adult')] 
        Susc_C = child_pop #value = pop size
        Susc_A = adult_pop
        Infc_C = 0
        Infc_A = 0 # number of infected adults
        Recv_C = 0
        Recv_A = 0 # number of recovered adults (empty for now)
        d_Susc[(met_id, 'C')] = Susc_C
        d_Susc[(met_id, 'A')] = Susc_A
        d_Infc[(met_id, 'C')] = Infc_C
        d_Infc[(met_id, 'A')] = Infc_A
        d_Recv[(met_id, 'C')] = Recv_C
        d_Recv[(met_id, 'A')] = Recv_A    
        
    return d_Susc, d_Infc, d_Recv
        
###################################################
def update_SI(metro_zero, met_id, d_Susc, d_Infc, num_new_infc_child, num_new_infc_adult, d_metro_infected_child, d_metro_infected_adult, time_step):
# updates dictionary values for d_Susc, d_Infc
    
    #children
    d_Susc[(met_id, 'C')] -= num_new_infc_child
    d_Infc[(met_id, 'C')] += num_new_infc_child
    num_infc_child = d_Infc[(met_id, 'C')]
    
    #adults
    d_Susc[(met_id, 'A')] -= num_new_infc_adult
    d_Infc[(met_id, 'A')] += num_new_infc_adult
    num_infc_adult = d_Infc[(met_id, 'A')]

    # record number of current infecteds in metro area at this time step 
    # not just those that were newly infected at this time step
    d_metro_infected_child[(metro_zero, met_id, time_step)] = num_infc_child
    d_metro_infected_adult[(metro_zero, met_id, time_step)] = num_infc_adult
    
###################################################
def update_IR(metro_zero, met_id, d_Infc, d_Recv, new_recov_child, new_recov_adult, d_metro_infected_child, d_metro_infected_adult, time_step):
# updates dictionary values for d_Infc, d_Recv, and d_metro_infected_child, _adult   
     
    #child
    d_Infc[(met_id, 'C')] -= new_recov_child
    d_Recv[(met_id, 'C')] += new_recov_child
    num_infc_child = d_Infc[(met_id, 'C')] # number of children infected at this time step in this metro area

    #adults
    d_Infc[(met_id, 'A')] -= new_recov_adult
    d_Recv[(met_id, 'A')] += new_recov_adult
    num_infc_adult = d_Infc[(met_id, 'A')]
    
    # record number of current infecteds in metro area at this time step 
    # not just those that were newly infected at this time step
    d_metro_infected_child[(metro_zero, met_id, time_step)] = num_infc_child
    d_metro_infected_adult[(metro_zero, met_id, time_step)] = num_infc_adult
    #d_metro_tot_infected_child[(met_id)] = 
    
###################################################
def travel_btwn_metros(air_network, d_Susc, d_Infc, d_Recv, d_prob_travel_C, d_prob_travel_A, theta_susc, theta_infc, theta_recv):
    
    edges = air_network.edges() 
    for (i, j) in edges:
        # Si <-> Sj
        
        # children
        
        # travel metro i --> j
        ch_susc_i = (d_Susc[(i, 'C')]) # binomial number of trials
        ch_susc_prob_i_j = ((d_prob_travel_C[(i, j)]) * (theta_susc)) # binomial probability
        # select number of children who travel
        ch_travel_i_j = ((ch_susc_i) * (ch_susc_prob_i_j))
        # multiply ch_susc_i by prob = number of children travel from i to j
        
        # travel metro j --> i
        ch_susc_j = (d_Susc[(j, 'C')])
        ch_susc_prob_j_i = ((d_prob_travel_C[(j, i)]) * (theta_susc))
        # select number of children who travel
        ch_travel_j_i = ((ch_susc_j) * (ch_susc_prob_j_i))
        
        # update pop sizes
        net_travel = ((ch_travel_j_i) - (ch_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)
        d_Susc[(i, 'C')] = ((d_Susc[(i, 'C')]) + net_travel)
        d_Susc[(j, 'C')] = ((d_Susc[(j, 'C')]) - net_travel)
	        
        #adults
        
        # travel metro i --> j
        ad_susc_i = (d_Susc[(i, 'A')])
        ad_susc_prob_i_j = ((d_prob_travel_A[(i, j)]) * (theta_susc))
        # select number of adults who travel
        ad_travel_i_j = ((ad_susc_i) * (ad_susc_prob_i_j))
        
        # travel metro j --> i
        ad_susc_j = (d_Susc[(j, 'A')])
        ad_susc_prob_j_i = ((d_prob_travel_A[(j, i)]) * (theta_susc))
        # select number of adults who travel
        ad_travel_j_i = ((ad_susc_j) * (ad_susc_prob_j_i))
        
        # update pop sizes
        net_travel = ((ad_travel_j_i) - (ad_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)	        
        d_Susc[(i, 'A')] = ((d_Susc[(i, 'A')]) + net_travel)
        d_Susc[(j, 'A')] = ((d_Susc[(j, 'A')]) - net_travel)
               
               
        # Ii <-> Ij
        
        # children
        
        # travel metro i --> j
        ch_infc_i = (d_Infc[(i, 'C')]) # binomial number of trials
        ch_infc_prob_i_j = ((d_prob_travel_C[(i, j)]) * (theta_infc)) # binomial probability
        # select number of children who travel
        ch_travel_i_j = ((ch_infc_i) * (ch_infc_prob_i_j))
        
        # travel metro j --> i
        ch_infc_j = (d_Infc[(j, 'C')])
        ch_infc_prob_j_i = ((d_prob_travel_C[(j, i)]) * (theta_infc))
        # select number of children who travel
        ch_travel_j_i = ((ch_infc_j) * (ch_infc_prob_j_i))
        
        # update pop sizes
        net_travel = ((ch_travel_j_i) - (ch_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)
        d_Infc[(i, 'C')] = ((d_Infc[(i, 'C')]) + net_travel)
        d_Infc[(j, 'C')] = ((d_Infc[(j, 'C')]) - net_travel)
	        
        #adults
        
        # travel metro i --> j
        ad_infc_i = (d_Infc[(i, 'A')])
        ad_infc_prob_i_j = ((d_prob_travel_A[(i, j)]) * (theta_infc))
        # select number of adults who travel
        ad_travel_i_j = ((ad_infc_i) * (ad_infc_prob_i_j))
        
        # travel metro j --> i
        ad_infc_j = (d_Infc[(j, 'A')])
        ad_infc_prob_j_i = ((d_prob_travel_A[(j, i)]) * (theta_infc))
        # select number of adults who travel
        ad_travel_j_i = ((ad_infc_j) * (ad_infc_prob_j_i))
        
        # update pop sizes
        net_travel = ((ad_travel_j_i) - (ad_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)	        
        d_Infc[(i, 'A')] = ((d_Infc[(i, 'A')]) + net_travel)
        d_Infc[(j, 'A')] = ((d_Infc[(j, 'A')]) - net_travel)


        # Ri <-> Rj
        
        # children
        
        # travel metro i --> j
        ch_recv_i = (d_Recv[(i, 'C')]) # binomial number of trials
        ch_recv_prob_i_j = ((d_prob_travel_C[(i, j)]) * (theta_recv))
        # select number of children who travel
        ch_travel_i_j = ((ch_recv_i) * (ch_recv_prob_i_j))
        
        # travel metro j --> i
        ch_recv_j = (d_Recv[(j, 'C')])
        ch_recv_prob_j_i = ((d_prob_travel_C[(j, i)]) * (theta_recv))
        # select number of children who travel
        ch_travel_j_i = ((ch_recv_j) * (ch_recv_prob_j_i))
        
        # update pop sizes
        net_travel = ((ch_travel_j_i) - (ch_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)
        d_Recv[(i, 'C')] = ((d_Recv[(i, 'C')]) + net_travel)
        d_Recv[(j, 'C')] = ((d_Recv[(j, 'C')]) - net_travel)
	        
        #adults
        
        # travel metro i --> j
        ad_recv_i = (d_Recv[(i, 'A')])
        ad_recv_prob_i_j = ((d_prob_travel_A[(i, j)]) * (theta_recv))
        # select number of adults who travel
        ad_travel_i_j = ((ad_recv_i) * (ad_recv_prob_i_j))
        
        # travel metro j --> i
        ad_recv_j = (d_Recv[(j, 'A')])
        ad_recv_prob_j_i = ((d_prob_travel_A[(j, i)]) * (theta_recv))
        # select number of adults who travel
        ad_travel_j_i = ((ad_recv_j) * (ad_recv_prob_j_i))
        
        # update pop sizes
        net_travel = ((ad_travel_j_i) - (ad_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)	        
        d_Recv[(i, 'A')] = ((d_Recv[(i, 'A')]) + net_travel)
        d_Recv[(j, 'A')] = ((d_Recv[(j, 'A')]) - net_travel)
        
        
###################################################
def chain_binomial_one_simulation(d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, metro_zero, d_metro_age_pop, d_metropop, metro_ids, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, air_network, ch_travelers_r, ad_travelers_s, num_metro_zeros, num_child_zeros, num_adult_zeros, params_to_calc_C, manip_exp_params):
# given a per contact, per time step tranmssion probability (beta) and
#  a per time step recovery probability (gamma), use the SIR chain binomial
#  model to simulate ONE outbreak on the population (provided by contact_network)
#  and return the total number of infected individuals (num_infected)
        
    #create dicts with initial pops of S, I, and R for each metro and age
    # keys: (metro, 'age'), value: pop in #
    d_Susc, d_Infc, d_Recv = SIR_initial_pops(metro_ids, d_metro_age_pop)
    
    time_step = 0 # clock counter keeping track of current time
    d_nat_infected_child, d_nat_infected_adult = {}, {} # this dictionary will keep track of how many infected in current time step
    #record number currently infected at each metro at each time step


    d_new_cases_child, d_new_cases_adult = {}, {}
    # record new cases for each time step (agg over metro areas)
    d_tot_new_cases_child, d_tot_new_cases_adult = {}, {}
    
    # infect patient_zeros and upated Susc and Infc lists
    metro_zeros = []
    metro_zeros.append(metro_zero)
    #second patient_zero selection for indv - select fixed number of patient_zeros - #children, #adults
    for met_id in metro_zeros:
        update_SI(met_id, met_id, d_Susc, d_Infc, num_child_zeros, num_adult_zeros, d_metro_infected_child, d_metro_infected_adult, time_step)
        num_infc_child = d_Infc[(met_id, 'C')]
        num_infc_adult = d_Infc[(met_id, 'A')]
        d_metro_tot_infected_child[(met_id, met_id, time_step)] = num_infc_child
        d_metro_tot_infected_adult[(met_id, met_id, time_step)] = num_infc_adult
        d_new_cases_child[(met_id, time_step)] = num_infc_child
        d_new_cases_adult[(met_id, time_step)] = num_infc_adult
       
    metros_not_zeros = [metro for metro in metro_ids if metro not in metro_zeros] 
    for met_id in metros_not_zeros:
        d_metro_infected_child[(metro_zero, met_id, time_step)] = 0
        d_metro_infected_adult[(metro_zero, met_id, time_step)] = 0
        d_metro_tot_infected_child[(metro_zero, met_id, time_step)] = 0
        d_metro_tot_infected_adult[(metro_zero, met_id, time_step)] = 0
        d_new_cases_child[(met_id, time_step)] = 0
        d_new_cases_adult[(met_id, time_step)] = 0
    
    num_infc_child = sum([d_Infc[(met_id, 'C')] for met_id in metro_ids])
    num_infc_adult = sum([d_Infc[(met_id, 'A')] for met_id in metro_ids])
    d_nat_infected_child[(time_step)] = num_infc_child
    d_nat_infected_adult[(time_step)] = num_infc_adult
    
    #experiment parameters
    C, _ = orig_params_to_exp(params_to_calc_C)
    experiment, intervention, timing, length_before, length_after, beta_exp = manip_exp_params
    exp_param_list, _ = def_exp_params(experiment, intervention, timing, length_before, length_after, beta_exp, C)
    experiment, time_exp_start, time_exp_end, beta_exp, C_cc_red = exp_param_list
    
    # while there are infected individuals
    # go to next time step
    num_infected = ((d_nat_infected_child[time_step]) + (d_nat_infected_adult[time_step])) 
    while num_infected > 0 and time_step < time_end:
        
        print "time %s has %s child, %s adult, %s total infections" % (time_step, d_nat_infected_child[time_step], d_nat_infected_adult[time_step], num_infected)
        
        time_step += 1 #update clock
    
        # TRAVEL #
    
        # create two dictionaries with probabilities of travel for each age group, keys being tuples of cities: (i, j) and (j, i)
        d_prob_travel_C, d_prob_travel_A = pop_func.calc_prob_travel(air_network, alpha, ch_travelers_r, ad_travelers_s, d_metropop)
        
        #update population sizes for S, I, R for each metro
        travel_btwn_metros(air_network, d_Susc, d_Infc, d_Recv, d_prob_travel_C, d_prob_travel_A, theta_susc, theta_infc, theta_recv)
        
        # DISEASE #
        
        #set control values for contact matrix and beta
        C, beta = orig_params_to_exp(params_to_calc_C)
        #print C
        
        #create list of interventions
        _, exp_inter_list = def_exp_params(experiment, intervention, timing, length_before, length_after, beta_exp, C)
        
        #set experimental parameters
        if time_step in range(time_exp_start, time_exp_end):
            while exp_inter_list:
                inter_apply = exp_inter_list.pop()
                if inter_apply == 'red_C_cc':
                    C = C_cc_red
                elif inter_apply == 'red_beta':
                    beta = beta_exp
            
        for met_id in metro_ids:
            
            

            # Ii --> Ri
            # determine how many child / adult infected will recover in each metro area
            
            #child
            Infc_C = d_Infc[(met_id, 'C')] # number of infected children in metro area
            prob = gamma # probability of recovery = gamma
            new_recov_child = ((Infc_C) * (prob))
            
            #adult
            Infc_A = d_Infc[(met_id, 'A')]
            prob = gamma
            new_recov_adult = ((Infc_A) * (prob))
            # multiply Infc_A by prob = new_recov_adult
            # don't round
            
            # subtract from Ii, add to Ri            
            update_IR(metro_zero, met_id, d_Infc, d_Recv, new_recov_child, new_recov_adult, d_metro_infected_child, d_metro_infected_adult, time_step)
               
            # Si --> Ii
            # determine how many susceptibles get infected in each metro area
            
            # child
            Susc_C = d_Susc[(met_id, 'C')] # number of susc children in metro area
            #assign C back to original C?
            prob = lambda_child_calc(C, Infc_A, Infc_C, met_id, d_metro_age_pop, beta) 
            #print "time %s contact matrix is %s" % (time_step, C)
            #else:
            #    prob = lambda_child_calc(C, Infc_A, Infc_C, met_id, d_metro_age_pop, beta)
                 
                                           
            #if experiment == 'yes':
            #    if intervention == 'red_C_cc':
            #        if time_step in range(time_exp_start, time_exp_end):
            #            prob = lambda_child_calc(C_cc_red, Infc_A, Infc_C, met_id, d_metro_age_pop, beta)
            #        else:    
            #            prob = lambda_child_calc(C, Infc_A, Infc_C, met_id, d_metro_age_pop, beta)
            #    elif intervention == 'red_beta':
            #        if time_step in range(time_exp_start, time_exp_end):
            #            prob = lambda_child_calc(C, Infc_A, Infc_C, met_id, d_metro_age_pop, beta_exp)
            #        else:    
            #            prob = lambda_child_calc(C, Infc_A, Infc_C, met_id, d_metro_age_pop, beta) # calc probability of infection
            #else:
            #    prob = lambda_child_calc(C, Infc_A, Infc_C, met_id, d_metro_age_pop, beta) # calc probability of infection
            #    
            # determine how many are infected (coin flip 'Susc_C' number of times, with probability 'prob' each flip will result in infected)
            new_cases_child = ((Susc_C) * (prob)) # determine how many are infected
            d_new_cases_child[(met_id, time_step)] = new_cases_child
            previous_time_step = (time_step - 1)
            previous_cases = d_metro_tot_infected_child[(metro_zero, met_id, previous_time_step)]
            d_metro_tot_infected_child[(metro_zero, met_id, time_step)] = previous_cases + new_cases_child
                         
            #adult
            Susc_A = d_Susc[(met_id, 'A')] # number of susc adults in metro area
            prob = lambda_adult_calc(C, Infc_A, Infc_C, met_id, d_metro_age_pop, beta) # calc probability of infection
            new_cases_adult = ((Susc_A) * (prob))
            d_new_cases_adult[(met_id, time_step)] = new_cases_adult
            previous_time_step = (time_step - 1)
            previous_cases = d_metro_tot_infected_adult[(metro_zero, met_id, previous_time_step)]
            d_metro_tot_infected_adult[(metro_zero, met_id, time_step)] = previous_cases + new_cases_adult
            #don't sum over total all time steps - either grab time-1 or sum over new cases
                    
            #subtract from Si, add to Ii
            update_SI(metro_zero, met_id, d_Susc, d_Infc, new_cases_child, new_cases_adult, d_metro_infected_child, d_metro_infected_adult, time_step)

        #record how many total infected across metro_ids at this time step
        num_infc_child = sum([d_Infc[(met_id, 'C')] for met_id in metro_ids])
        num_infc_adult = sum([d_Infc[(met_id, 'A')] for met_id in metro_ids])
        d_nat_infected_child[(time_step)] = num_infc_child
        d_nat_infected_adult[(time_step)] = num_infc_adult
        num_infected = ((d_nat_infected_child[time_step]) + (d_nat_infected_adult[time_step])) 
        d_tot_new_cases_adult[(time_step)] = sum([d_new_cases_adult[(met_id, time_step)] for met_id in metro_ids])
        d_tot_new_cases_child[(time_step)] = sum([d_new_cases_child[(met_id, time_step)] for met_id in metro_ids])        
        #print d_tot_new_cases_child[(time_step)]
        # go back thru while loop
        # next time step
        # travel again
        # S --> I --> R
        
        
        
    # Note num_newly_infected is the incidence time series      
    return d_new_cases_child, d_new_cases_adult, d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, sum(d_tot_new_cases_child.values()), sum(d_tot_new_cases_adult.values()) # return total number infected in outbreak

###################################################
def chain_binomial_monte_carlo(R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params):
# run many (num_sims) instances of chain binomial simulation and return average epidemic size    

    d_metro_infected_child, d_metro_infected_adult = {}, {}
    #record number total infected at each metro area
    d_metro_tot_infected_child, d_metro_tot_infected_adult = {}, {}
    # record previous new cases plus current new cases at each time step and each metro
    
    #d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
    population_size = sum([d_metropop[x] for x in metro_ids])
    threshold = 0.10 # 10% of population size is our threshold for a large epidemic

    d_metro_age_pop = pop_func.calc_metro_age_pop(filename_metropop, alpha)

    large_epidemic_sizes= [] # will keep list of outbreak sizes that are large epidemics
    for metro in metro_ids:

        new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, outbreak_size_child, outbreak_size_adult = chain_binomial_one_simulation(d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, metro, d_metro_age_pop, d_metropop, metro_ids, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, air_network, ch_travelers_r, ad_travelers_s, num_metro_zeros, num_child_zeros, num_adult_zeros, params_to_calc_C, manip_exp_params)

        outbreak_size = (outbreak_size_child + outbreak_size_adult)
        
        ## SB said don't worry about this for now ##
        # figure out if this is small outbreak or large epidemic
        #if outbreak_size > threshold * population_size: # if outbreak reached more than 10% of the population
            #large_epidemic_sizes.append(outbreak_size)

        # call it a large epidemic and save its size
        large_epidemic_sizes.append(outbreak_size)
        
        ##output - plot for each sim
        #plot_new_cases(metro_ids, time_end, new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, sim, alpha, ch_travelers_r, R0, gamma, beta)
        #plot_current_cases(metro_ids, time_end, incidence_time_series_metro_child, incidence_time_series_metro_adult, sim, alpha, ch_travelers_r, R0, gamma, beta)
            
    # calculate average large epidemic size, and how frequent they were
    if large_epidemic_sizes:
        average_epidemic_size = np.mean(large_epidemic_sizes)/float(population_size)
        epi_size_fractions = [x / float(population_size) for x in large_epidemic_sizes] #divide cases by pop to get percent of pop infected
        standard_deviation = np.std(epi_size_fractions)
        #probability_epidemic = len(large_epidemic_sizes)/float(num_sims)
    else:
        average_epidemic_size = 0
        #probability_epidemic = 0
    
    
    return average_epidemic_size, standard_deviation, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult

###################################################
def chain_binomial_monte_carlo_plots_csv(R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, intervention, manip_exp_params):
# run many (num_sims) instances of chain binomial simulation and return average epidemic size    

    d_metro_infected_child, d_metro_infected_adult = {}, {}
    #record number total infected at each metro area
    d_metro_tot_infected_child, d_metro_tot_infected_adult = {}, {}
    # record previous new cases plus current new cases at each time step and each metro
    
    #d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
    population_size = sum([d_metropop[x] for x in metro_ids])
    threshold = 0.10 # 10% of population size is our threshold for a large epidemic

    d_metro_age_pop = pop_func.calc_metro_age_pop(filename_metropop, alpha)

    large_epidemic_sizes, adult_epidemic_sizes, child_epidemic_sizes = [], [], [] # will keep list of outbreak sizes that are large epidemics
    for metro in metro_ids:
        
        new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, outbreak_size_child, outbreak_size_adult = chain_binomial_one_simulation(d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, metro, d_metro_age_pop, d_metropop, metro_ids, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, air_network, ch_travelers_r, ad_travelers_s, num_metro_zeros, num_child_zeros, num_adult_zeros, params_to_calc_C, manip_exp_params)
        
        #d_sim[sim] = incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult
        
        outbreak_size = (outbreak_size_child + outbreak_size_adult)
        
        ## SB said don't worry about this for now ##
        # figure out if this is small outbreak or large epidemic
        #if outbreak_size > threshold * population_size: # if outbreak reached more than 10% of the population
            #large_epidemic_sizes.append(outbreak_size)

        # call it a large epidemic and save its size
        large_epidemic_sizes.append(outbreak_size)
        adult_epidemic_sizes.append(outbreak_size_adult)
        child_epidemic_sizes.append(outbreak_size_child)
        #output - plot for each sim
        #plot_new_cases(metro_ids, time_end, new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, metro, alpha, ch_travelers_r, R0, gamma, beta, intervention)
        #plot_current_cases(metro_ids, time_end, incidence_time_series_metro_child, incidence_time_series_metro_adult, metro, alpha, ch_travelers_r, R0, gamma, beta)
            
    # calculate average large epidemic size, and how frequent they were
    if large_epidemic_sizes:
        average_epidemic_size = np.mean(large_epidemic_sizes)/float(population_size)
        epi_size_fractions = [x / float(population_size) for x in large_epidemic_sizes] #divide cases by pop to get percent of pop infected
        standard_deviation = np.std(epi_size_fractions)
        #probability_epidemic = len(large_epidemic_sizes)/float(num_sims)
    else:
        average_epidemic_size = 0
        #probability_epidemic = 0
        
    # calculate average adult epi size
    if adult_epidemic_sizes:
        avg_adult_epi_size = np.mean(adult_epidemic_sizes)/float(population_size-(population_size*alpha)) # adult pop = total pop - child pop
    else:
        avg_adult_epi_size = 0
    
    # calculate average child epi size
    if child_epidemic_sizes:
        avg_child_epi_size = np.mean(child_epidemic_sizes)/float(population_size*alpha) #child pop = alpha * total pop
    else:
        avg_child_epi_size = 0
    
    return average_epidemic_size, standard_deviation, avg_adult_epi_size, avg_child_epi_size, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult

###################################################
def chain_binomial_monte_carlo_unit_tests(unit_test, R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params):
# run many (num_sims) instances of chain binomial simulation and return average epidemic size    

    d_metro_infected_child, d_metro_infected_adult = {}, {}
    #record number total infected at each metro area
    d_metro_tot_infected_child, d_metro_tot_infected_adult = {}, {}
    # record previous new cases plus current new cases at each time step and each metro
    
    #d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
    population_size = sum([d_metropop[x] for x in metro_ids])
    threshold = 0.10 # 10% of population size is our threshold for a large epidemic

    d_metro_age_pop = pop_func.calc_metro_age_pop(filename_metropop, alpha)

    large_epidemic_sizes= [] # will keep list of outbreak sizes that are large epidemics
    for metro_zero in metro_ids:

        new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, outbreak_size_child, outbreak_size_adult = chain_binomial_one_simulation(d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, metro_zero, d_metro_age_pop, d_metropop, metro_ids, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, air_network, ch_travelers_r, ad_travelers_s, num_metro_zeros, num_child_zeros, num_adult_zeros, params_to_calc_C, manip_exp_params)

        outbreak_size = (outbreak_size_child + outbreak_size_adult)
        
        ## SB said don't worry about this for now ##
        # figure out if this is small outbreak or large epidemic
        #if outbreak_size > threshold * population_size: # if outbreak reached more than 10% of the population
            #large_epidemic_sizes.append(outbreak_size)

        # call it a large epidemic and save its size
        large_epidemic_sizes.append(outbreak_size)
        
        #output - plot for each metro_zero
        plot_new_cases_unit_tests(unit_test, metro_ids, time_end, new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, metro_zero, alpha, ch_travelers_r, R0, gamma, beta)
        #plot_current_cases(metro_ids, time_end, incidence_time_series_metro_child, incidence_time_series_metro_adult, sim)
        #write_csv_file(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult)
            
    # calculate average large epidemic size, and how frequent they were
    if large_epidemic_sizes:
        average_epidemic_size = np.mean(large_epidemic_sizes)/float(population_size)
        #probability_epidemic = len(large_epidemic_sizes)/float(num_sims)
    else:
        average_epidemic_size = 0
        #probability_epidemic = 0
    
    
    return average_epidemic_size

###################################################
def chain_binomial_monte_carlo_unit_tests_csv(unit_test, R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params):
# run many (num_sims) instances of chain binomial simulation and return average epidemic size    

    d_metro_infected_child, d_metro_infected_adult = {}, {}
    #record number total infected at each metro area
    d_metro_tot_infected_child, d_metro_tot_infected_adult = {}, {}
    # record previous new cases plus current new cases at each time step and each metro
    
    #d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
    population_size = sum([d_metropop[x] for x in metro_ids])
    threshold = 0.10 # 10% of population size is our threshold for a large epidemic

    d_metro_age_pop = pop_func.calc_metro_age_pop(filename_metropop, alpha)

    large_epidemic_sizes, adult_epidemic_sizes, child_epidemic_sizes = [], [], [] # will keep list of outbreak sizes that are large epidemics
    for metro_zero in metro_ids:

        new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, outbreak_size_child, outbreak_size_adult = chain_binomial_one_simulation(d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, metro_zero, d_metro_age_pop, d_metropop, metro_ids, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, air_network, ch_travelers_r, ad_travelers_s, num_metro_zeros, num_child_zeros, num_adult_zeros, params_to_calc_C, manip_exp_params)

        outbreak_size = (outbreak_size_child + outbreak_size_adult)
        
        ## SB said don't worry about this for now ##
        # figure out if this is small outbreak or large epidemic
        #if outbreak_size > threshold * population_size: # if outbreak reached more than 10% of the population
            #large_epidemic_sizes.append(outbreak_size)

        # call it a large epidemic and save its size
        large_epidemic_sizes.append(outbreak_size)
        adult_epidemic_sizes.append(outbreak_size_adult)
        child_epidemic_sizes.append(outbreak_size_child)
        
        #output - plot for each metro_zero
        #plot_new_cases_unit_tests(unit_test, metro_ids, time_end, new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, metro_zero, alpha, ch_travelers_r, R0, gamma, beta)
        #plot_current_cases(metro_ids, time_end, incidence_time_series_metro_child, incidence_time_series_metro_adult, sim)
        #write_csv_file(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult)
            
    # calculate average large epidemic size, and how frequent they were
    if large_epidemic_sizes:
        average_epidemic_size = np.mean(large_epidemic_sizes)/float(population_size)
        #probability_epidemic = len(large_epidemic_sizes)/float(num_sims)
    else:
        average_epidemic_size = 0
        #probability_epidemic = 0
        
    # calculate average adult epi size
    if adult_epidemic_sizes:
        avg_adult_epi_size = np.mean(adult_epidemic_sizes)/float(population_size-(population_size*alpha))
    else:
        avg_adult_epi_size = 0
    
    # calculate average child epi size
    if child_epidemic_sizes:
        avg_child_epi_size = np.mean(child_epidemic_sizes)/float(population_size*alpha)
    else:
        avg_child_epi_size = 0
    
    return average_epidemic_size, avg_adult_epi_size, avg_child_epi_size, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult

###################################################
def chain_binomial_monte_carlo_unit_tests_csv_test(unit_test, R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params):
# run many (num_sims) instances of chain binomial simulation and return average epidemic size    

    d_metro_infected_child, d_metro_infected_adult = {}, {}
    #record number total infected at each metro area
    d_metro_tot_infected_child, d_metro_tot_infected_adult = {}, {}
    # record previous new cases plus current new cases at each time step and each metro
    
    #d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
    population_size = sum([d_metropop[x] for x in metro_ids])
    threshold = 0.10 # 10% of population size is our threshold for a large epidemic

    d_metro_age_pop = pop_func.calc_metro_age_pop(filename_metropop, alpha)

    large_epidemic_sizes, adult_epidemic_sizes, child_epidemic_sizes = [], [], [] # will keep list of outbreak sizes that are large epidemics
    metro_zero = 1

    new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, outbreak_size_child, outbreak_size_adult = chain_binomial_one_simulation(d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, metro_zero, d_metro_age_pop, d_metropop, metro_ids, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, air_network, ch_travelers_r, ad_travelers_s, num_metro_zeros, num_child_zeros, num_adult_zeros, params_to_calc_C, manip_exp_params)

    outbreak_size = (outbreak_size_child + outbreak_size_adult)
    
    time_list = range(1, time_end+1)
    
        
    metro_zero_outbreak_adult = sum([new_cases_incidence_time_series_metro_adult[(metro_zero, time_step)] for time_step in time_list])

    metro_zero_outbreak_child = sum([new_cases_incidence_time_series_metro_child[(metro_zero, time_step)] for time_step in time_list])
    
    metro_zero_outbreak = (metro_zero_outbreak_adult + metro_zero_outbreak_child)
        
        ## SB said don't worry about this for now ##
        # figure out if this is small outbreak or large epidemic
        #if outbreak_size > threshold * population_size: # if outbreak reached more than 10% of the population
            #large_epidemic_sizes.append(outbreak_size)

        # call it a large epidemic and save its size
    large_epidemic_sizes.append(metro_zero_outbreak)
    adult_epidemic_sizes.append(metro_zero_outbreak_adult)
    child_epidemic_sizes.append(metro_zero_outbreak_child)
        
        #output - plot for each metro_zero
        #plot_new_cases_unit_tests(unit_test, metro_ids, time_end, new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, metro_zero, alpha, ch_travelers_r, R0, gamma, beta)
        #plot_current_cases(metro_ids, time_end, incidence_time_series_metro_child, incidence_time_series_metro_adult, sim)
        #write_csv_file(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult)
            
    # calculate average large epidemic size, and how frequent they were
    metro_zero_pop = d_metropop[metro_zero]
    if large_epidemic_sizes:
        average_epidemic_size = np.mean(large_epidemic_sizes)/float(metro_zero_pop)
        #probability_epidemic = len(large_epidemic_sizes)/float(num_sims)
    else:
        average_epidemic_size = 0
        #probability_epidemic = 0
        
    # calculate average adult epi size
    if adult_epidemic_sizes:
        avg_adult_epi_size = np.mean(adult_epidemic_sizes)/float(metro_zero_pop-(metro_zero_pop*alpha))
    else:
        avg_adult_epi_size = 0
    
    # calculate average child epi size
    if child_epidemic_sizes:
        avg_child_epi_size = np.mean(child_epidemic_sizes)/float(metro_zero_pop*alpha)
    else:
        avg_child_epi_size = 0
    
    return average_epidemic_size, avg_adult_epi_size, avg_child_epi_size, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult
               
###################################################
def chain_binomial_monte_carlo_unit_tests_no_plot(unit_test, R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params):
# run many (num_sims) instances of chain binomial simulation and return average epidemic size    

    d_metro_infected_child, d_metro_infected_adult = {}, {}
    #record number total infected at each metro area
    d_metro_tot_infected_child, d_metro_tot_infected_adult = {}, {}
    # record previous new cases plus current new cases at each time step and each metro
    
    #d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
    population_size = sum([d_metropop[x] for x in metro_ids])
    threshold = 0.10 # 10% of population size is our threshold for a large epidemic

    d_metro_age_pop = pop_func.calc_metro_age_pop(filename_metropop, alpha)

    large_epidemic_sizes= [] # will keep list of outbreak sizes that are large epidemics
    for metro in metro_ids:
        #incidence_time_series, outbreak_size = chain_binomial_one_simulation(d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, sim, d_metro_age_pop, d_metropop, metro_ids, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, air_network, ch_travelers_r, ad_travelers_s, num_metro_zeros, num_child_zeros, num_adult_zeros, params_to_calc_C, manip_exp_params)
        # Note we are not using the incidence time series right now
        
        new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, outbreak_size_child, outbreak_size_adult = chain_binomial_one_simulation(d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, metro, d_metro_age_pop, d_metropop, metro_ids, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, air_network, ch_travelers_r, ad_travelers_s, num_metro_zeros, num_child_zeros, num_adult_zeros, params_to_calc_C, manip_exp_params)
        print outbreak_size_child
        print outbreak_size_adult
        outbreak_size = (outbreak_size_child + outbreak_size_adult)
        
        ## SB said don't worry about this for now ##
        # figure out if this is small outbreak or large epidemic
        #if outbreak_size > threshold * population_size: # if outbreak reached more than 10% of the population
            #large_epidemic_sizes.append(outbreak_size)

        # call it a large epidemic and save its size
        large_epidemic_sizes.append(outbreak_size)
        
        ##output - plot for each sim
        #plot_new_cases(metro_ids, time_end, new_cases_incidence_time_series_metro_child, new_cases_incidence_time_series_metro_adult, sim, alpha, ch_travelers_r, R0, gamma, beta)
        #plot_current_cases(metro_ids, time_end, incidence_time_series_metro_child, incidence_time_series_metro_adult, sim, alpha, ch_travelers_r, R0, gamma, beta)
        #write_csv_file(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult)
            
    # calculate average large epidemic size, and how frequent they were
    if large_epidemic_sizes:
        average_epidemic_size = np.mean(large_epidemic_sizes)/float(population_size)
        #probability_epidemic = len(large_epidemic_sizes)/float(num_sims)
    else:
        average_epidemic_size = 0
        #probability_epidemic = 0
    
    
    return average_epidemic_size    
    
###################################################
def read_edgelist_anne (filename):
# read metro area contact network from edgelist contained in file (filename)
    
    G = nx.Graph()
    file = open(filename, 'rU')
    reader = csv.reader(file)
    for row in reader:
        data = row[0].split('\t')
        G.add_edge(int(data[0]),int(data[1]), weight=float(data[2]))
        
    return G
    
####################################################
#def calc_avg_excess_degree(air_network):
#    deg_list = air_network.degree().values()
#    deg_list_minus_one = [d-1 for d in deg_list]
#    deg_times_deg_minus_one = [d*d1 for d,d1 in zip(deg_list, deg_list_minus_one)]
#    numerator = np.mean(deg_times_deg_minus_one)
#    denominator = np.mean(deg_list)
#        
#    return float(numerator)/denominator
#
####################################################
#def calculate_beta(R0, gamma, air_network):
#    
#    # calculate T value for given R0
#    #  where R0 = T * avg_excees_degree 
#    avg_excess_deg = calc_avg_excess_degree(air_network)
#    T = float(R0)/avg_excess_deg
#            
#    # calculate beta based on T = beta/(beta + gamma)
#    beta = T*gamma/(1-T)
#        
#    return beta

####################################################
def write_csv_file(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult):
    
    #time_series = range(0, time_end)
    csvfile = csv.writer(open('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/chain_binomial_output_nummetrozeros_%s_numchildzeros_%s_numadultzeros_%s.csv' % (num_metro_zeros, num_child_zeros, num_adult_zeros), 'wb'), delimiter = ',')
    csvfile.writerow(['metro_zero', 'time_step', 'metro_id', 'age', 'currently_infected', 'total_infected'])
    list_tuples = []
    for (met_zero, met_id, time_step) in incidence_time_series_metro_child:
        list_tuples.append((met_zero, met_id, time_step))
    sort_by_time_list_tuples = sorted(list_tuples, key=operator.itemgetter(2))
    sort_by_sim_list_tuples = sorted(sort_by_time_list_tuples, key=operator.itemgetter(0))
    for (met_zero, met_id, time_step) in sort_by_sim_list_tuples:
        csvfile.writerow([met_zero, time_step, met_id, 'C', (incidence_time_series_metro_child[(met_zero, met_id, time_step)]), (tot_incidence_time_series_child[(met_zero, met_id, time_step)])])
    for (sim, met_id, time_step) in sort_by_sim_list_tuples:
        csvfile.writerow([met_zero, time_step, met_id, 'A', (incidence_time_series_metro_adult[(met_zero, met_id, time_step)]), (tot_incidence_time_series_adult[(met_zero, met_id, time_step)])])

####################################################
def write_csv_file_no_travel(num_metro_zeros, num_child_zeros, num_adult_zeros, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, beta):
    
    #time_series = range(0, time_end)
    csvfile = csv.writer(open('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Unit_Tests/no_travel/chain_binomial_output_nummetrozeros_%s_numchildzeros_%s_numadultzeros_%s_beta_%s.csv' % (num_metro_zeros, num_child_zeros, num_adult_zeros, beta), 'wb'), delimiter = ',')
    csvfile.writerow(['metro_zero', 'time_step', 'metro_id', 'age', 'currently_infected', 'total_infected'])
    list_tuples = []
    for (met_zero, met_id, time_step) in incidence_time_series_metro_child:
        list_tuples.append((met_zero, met_id, time_step))
    sort_by_time_list_tuples = sorted(list_tuples, key=operator.itemgetter(2))
    sort_by_sim_list_tuples = sorted(sort_by_time_list_tuples, key=operator.itemgetter(0))
    for (met_zero, met_id, time_step) in sort_by_sim_list_tuples:
        csvfile.writerow([met_zero, time_step, met_id, 'C', (incidence_time_series_metro_child[(met_zero, met_id, time_step)]), (tot_incidence_time_series_child[(met_zero, met_id, time_step)])])
    for (sim, met_id, time_step) in sort_by_sim_list_tuples:
        csvfile.writerow([met_zero, time_step, met_id, 'A', (incidence_time_series_metro_adult[(met_zero, met_id, time_step)]), (tot_incidence_time_series_adult[(met_zero, met_id, time_step)])])

####################################################
def write_csv_file_test_beta(num_metro_zeros, num_child_zeros, num_adult_zeros, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, beta):
    
    #time_series = range(0, time_end)
    csvfile = csv.writer(open('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Unit_Tests/test_beta_value/chain_binomial_output_nummetrozeros_%s_numchildzeros_%s_numadultzeros_%s_beta_%s_newcontactmatrix.csv' % (num_metro_zeros, num_child_zeros, num_adult_zeros, beta), 'wb'), delimiter = ',')
    csvfile.writerow(['metro_zero', 'time_step', 'metro_id', 'age', 'currently_infected', 'total_infected'])
    list_tuples = []
    for (met_zero, met_id, time_step) in incidence_time_series_metro_child:
        list_tuples.append((met_zero, met_id, time_step))
    sort_by_time_list_tuples = sorted(list_tuples, key=operator.itemgetter(2))
    sort_by_sim_list_tuples = sorted(sort_by_time_list_tuples, key=operator.itemgetter(0))
    for (met_zero, met_id, time_step) in sort_by_sim_list_tuples:
        csvfile.writerow([met_zero, time_step, met_id, 'C', (incidence_time_series_metro_child[(met_zero, met_id, time_step)]), (tot_incidence_time_series_child[(met_zero, met_id, time_step)])])
    for (sim, met_id, time_step) in sort_by_sim_list_tuples:
        csvfile.writerow([met_zero, time_step, met_id, 'A', (incidence_time_series_metro_adult[(met_zero, met_id, time_step)]), (tot_incidence_time_series_adult[(met_zero, met_id, time_step)])])

###################################################
def write_csv_file_experiments(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, beta, intervention):
    
    #time_series = range(0, time_end)
    csvfile = csv.writer(open('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Experiments/%s/chain_binomial_output_nummetrozeros_%s_numchildzeros_%s_numadultzeros_%s_beta_%s_intervention_%s.csv' % (intervention, num_metro_zeros, num_child_zeros, num_adult_zeros, beta, intervention), 'wb'), delimiter = ',')
    csvfile.writerow(['metro_zero', 'time_step', 'metro_id', 'age', 'currently_infected', 'total_infected'])
    list_tuples = []
    for (met_zero, met_id, time_step) in incidence_time_series_metro_child:
        list_tuples.append((met_zero, met_id, time_step))
    sort_by_time_list_tuples = sorted(list_tuples, key=operator.itemgetter(2))
    sort_by_sim_list_tuples = sorted(sort_by_time_list_tuples, key=operator.itemgetter(0))
    for (met_zero, met_id, time_step) in sort_by_sim_list_tuples:
        csvfile.writerow([met_zero, time_step, met_id, 'C', (incidence_time_series_metro_child[(met_zero, met_id, time_step)]), (tot_incidence_time_series_child[(met_zero, met_id, time_step)])])
    for (sim, met_id, time_step) in sort_by_sim_list_tuples:
        csvfile.writerow([met_zero, time_step, met_id, 'A', (incidence_time_series_metro_adult[(met_zero, met_id, time_step)]), (tot_incidence_time_series_adult[(met_zero, met_id, time_step)])])


###################################################
def plot_new_cases (metro_ids, time_end, d_new_cases_child, d_new_cases_adult, metro_zero, alpha, ch_travelers_r, R0, gamma, beta, intervention):
# print out time series of new cases at each time step
# one line for each metro area
# one plot for child, one for adult
        
# child
    for met_id in metro_ids:
        time_series = range(0, time_end)
        graphxax = time_series
        graphyax = [d_new_cases_child[(met_id, time_step)] for time_step in time_series]
        plt.plot(graphxax, graphyax)
        
    #plt.xlim([0, 9])
    #plt.ylim([0, 1000000])
    plt.xlabel('Time Step')
    plt.ylabel('New Cases - Child')
    #plt.xticks(range(0, 10), wklab)
    #plt.show()
    #plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/New_Cases/Child/chain_binomial_new_cases_child_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Experiments/%s/Plots/Child/chain_binomial_new_cases_child_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (intervention, alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
    plt.close()

#might need separate for adult
#adult
    for met_id in metro_ids:
        time_series = range(0, time_end)
        graphxax = time_series
        graphyax = [d_new_cases_adult[(met_id, time_step)] for time_step in time_series]
        plt.plot(graphxax, graphyax)
        
    #plt.xlim([0, 40])
    #plt.ylim([0, 1000000])
    plt.xlabel('Time Step')
    plt.ylabel('New Cases - Adult')
    #plt.xticks(range(0, 10), wklab)
    #plt.show()
    #plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/New_Cases/Adult/chain_binomial_new_cases_adult_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Experiments/%s/Plots/Adult/chain_binomial_new_cases_child_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (intervention, alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
    plt.close()
         

###################################################
def plot_current_cases(metro_ids, time_end, d_metro_infected_child, d_metro_infected_adult, metro_zero, alpha, ch_travelers_r, R0, gamma, beta):
# time series for currently infected cases
# one plot for child one for adult
# one line for each metro area

    
# child
    #for met_id in metro_ids:
    #    time_series = range(0, time_end)
    #    graphxax = time_series
    #    graphyax = [d_metro_infected_child[(metro_zero, met_id, time_step)] for time_step in time_series]
    #    plt.plot(graphxax, graphyax)
    #    
    ##plt.xlim([0, 9])
    ##plt.ylim([0, 1000000])
    #plt.xlabel('Time Step')
    #plt.ylabel('Current Cases - Child')
    ##plt.xticks(range(0, 10), wklab)
    ##plt.show()
    #plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Current_Cases/Child/chain_binomial_current_cases_child_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
    #plt.close()

#adult
    for met_id in metro_ids:
        time_series = range(0, time_end)
        graphxax = time_series
        graphyax = [d_metro_infected_adult[(metro_zero, met_id, time_step)] for time_step in time_series]
        plt.plot(graphxax, graphyax)
        
    #plt.xlim([0, 40])
    #plt.ylim([0, 1000])
    plt.xlabel('Time Step')
    plt.ylabel('Current Cases - Adult')
    #plt.xticks(range(0, 10), wklab)
    #plt.show()
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Current_Cases/Adult/chain_binomial_current_cases_adult_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
    plt.close()

###################################################
def plot_new_cases_unit_tests (unit_test, metro_ids, time_end, d_new_cases_child, d_new_cases_adult, metro_zero, alpha, ch_travelers_r, R0, gamma, beta):
# print out time series of new cases at each time step
# one line for each metro area
# one plot for child, one for adult
        
# child
    for met_id in metro_ids:
        time_series = range(0, time_end)
        graphxax = time_series
        graphyax = [d_new_cases_child[(met_id, time_step)] for time_step in time_series]
        plt.plot(graphxax, graphyax)
        
    #plt.xlim([0, 9])
    #plt.ylim([0, 1000000])
    plt.xlabel('Time Step')
    plt.ylabel('New Cases - Child')
    #plt.xticks(range(0, 10), wklab)
    #plt.show()
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Unit_Tests/%s/Child/chain_binomial_new_cases_child_unit_test_%s_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (unit_test, unit_test, alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
    plt.close()

#might need separate for adult
#adult
    for met_id in metro_ids:
        time_series = range(0, time_end)
        graphxax = time_series
        graphyax = [d_new_cases_adult[(met_id, time_step)] for time_step in time_series]
        plt.plot(graphxax, graphyax)
        
    #plt.xlim([0, 9])
    #plt.ylim([0, 1000000])
    plt.xlabel('Time Step')
    plt.ylabel('New Cases - Adult')
    #plt.xticks(range(0, 10), wklab)
    #plt.show()
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Unit_Tests/%s/Adult/chain_binomial_new_cases_adult_unit_test_%s_alpha_%1.2f_r_%s_R0_%s_gamma_%s_beta_%s_metrozero_%s.png' % (unit_test, unit_test, alpha, ch_travelers_r, R0, gamma, beta, metro_zero))
    plt.close()   
     
###################################################
def beta_test (betas, gamma, R0, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params):
    
    d_epi_size, d_epi_size_std = {}, {}
    for beta in betas:
        average_epidemic_size, std_dev, _, _, _, _ = chain_binomial_monte_carlo(R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params)
        #add epi_size to dict
	#divide by pop size to make it a prop of pop
	#population_size = sum([d_metropop[x] for x in metro_ids])
        #d_epi_size[(beta)] = (average_epidemic_size/population_size)
        d_epi_size[(beta)] = average_epidemic_size
	d_epi_size_std[(beta)] = std_dev
        
    graphxax = betas
    graphyax = [d_epi_size[(beta)] for beta in betas]
    #print graphyax
    yerr = [d_epi_size_std[(beta)] for beta in betas]
    #print yerr
    plt.errorbar(graphxax, graphyax, yerr=yerr)
    
    plt.xlabel('Beta Value')
    plt.ylabel('Epidemic Size') 
    #plt.show()
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Beta_v_Epi_Size/chain_binomial_beta_v_epidemic_size_alpha_%1.2f_r_%s_R0_%s_gamma_%s.png' % (alpha, ch_travelers_r, R0, gamma))
    plt.close()
   
###################################################  
def gamma_test (beta, gammas, R0, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params):
    
    d_epi_size, d_epi_size_std = {}, {}
    for gamma in gammas:
        average_epidemic_size, std_dev, _, _, _, _ = chain_binomial_monte_carlo(R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, manip_exp_params)
        #add epi_size to dict
        d_epi_size[(gamma)] = average_epidemic_size
        d_epi_size_std[(gamma)] = std_dev
        
    graphxax = gammas
    graphyax = [d_epi_size[(gamma)] for gamma in gammas]
    yerr = [d_epi_size_std[(gamma)] for gamma in gammas]
    #print yerr
    plt.errorbar(graphxax, graphyax, yerr=yerr)
    
    plt.xlabel('Gamma Value')
    plt.ylabel('Epidemic Size') 
    #plt.show()
    plt.savefig('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/Deterministic_Model/Diagnostic_Plots/Gamma_v_Epi_Size/chain_binomial_gamma_v_epidemic_size_alpha_%1.2f_r_%s_R0_%s_beta_%s.png' % (alpha, ch_travelers_r, R0, beta))
    plt.close()   

###################################################  
def orig_params_to_exp(params_to_calc_C):
    
    p_c, p_a, q_c, q_a, alpha = params_to_calc_C	
    C = pop_func.calc_contact_matrix_pqa(p_c, p_a, q_c, q_a, alpha)
    beta = 0.037
    
    return C, beta
    
###################################################    
def def_exp_params(experiment, intervention, timing, length_before, length_after, beta_exp, C):
    
    exp_param_list, exp_inter_list = [], [] #create lists for parameters, and list for interventions
    exp_param_list.append(experiment) # append experiment (T or F) to parameter list
    
    if experiment == 'yes':
        exp_inter_list.append(intervention)
        
    time_exp_start, time_exp_end = exp_func.set_time_start_and_end(timing, length_before, length_after)
    exp_param_list.append(time_exp_start)
    exp_param_list.append(time_exp_end)

    exp_param_list.append(beta_exp)
    
    if intervention == 'red_C_cc':
        C_cc_red = exp_func.reduce_C_cc(C)
        exp_param_list.append(C_cc_red)
    
    return exp_param_list, exp_inter_list  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
###################################################
if __name__ == "__main__":
    
    # READ METRO NETWORK FROM FILE
    filename_metropop = 'Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/metedges.txt'
    d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
    filename_air_network = 'Dropbox/Anne_Bansal_lab/Python_Scripts/Modeling_Project/air_traffic_data/air_traffic_edgelist.txt'
    air_network = read_edgelist_anne(filename_air_network)
    # READ US population data
    us_popdata = csv.reader(open('Dropbox/Anne_Bansal_lab/SDI_Data/totalpop_age.csv', 'r'),delimiter = ',')
    dict_popdata, ages, years = pop_func.import_popdata(us_popdata, 0, 1, 2)
    dict_childpop, dict_adultpop = pop_func.pop_child_adult (dict_popdata, years)
    # READ Germany contact data
    #filename_germ_contact_data = 'Dropbox/Anne_Bansal_lab/Contact_Data/polymod_germany_contact_matrix_Mossong_2008.csv'
    filename_germ_within_group_contact_data = 'Dropbox/Anne_Bansal_lab/Contact_Data/within_group_polymod_germany_contact_matrix_Mossong_2008.csv'
    filename_germ_all_contact_data = 'Dropbox/Anne_Bansal_lab/Contact_Data/all_ages_polymod_germany_contact_matrix_Mossong_2008.csv'
    # READ Germany population data
    filename_germ_pop_data = 'Dropbox/Anne_Bansal_lab/UNdata_Export_2008_Germany_Population.csv'
        
    # DEFINE POPULATION PARAMETERS
    year = 2010
    alpha = pop_func.calc_alpha(year, dict_childpop, dict_adultpop)
    d_metro_age_pop = pop_func.calc_metro_age_pop(filename_metropop, alpha)
    ch_travelers_r = 0.0 # fraction of children who travel
    ad_travelers_s = 1
    
    # CONTACT MATRIX
    q_c, q_a, p_c, p_a, _, _ = pop_func.calc_p(filename_germ_within_group_contact_data, filename_germ_pop_data, filename_germ_all_contact_data)
    params_to_calc_C = []
    params_to_calc_C.append(p_c)
    params_to_calc_C.append(p_a)
    params_to_calc_C.append(q_c)
    params_to_calc_C.append(q_a)
    params_to_calc_C.append(alpha)
    C = pop_func.calc_contact_matrix_pqa(p_c, p_a, q_c, q_a, alpha)
    #print C
                            
    # DEFINE DISEASE PARAMETERS
    R0 = 1.2
    gamma = 0.5 # recovery rate based on (1/gamma) day infectious period
    test_gammas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    #beta = 0.037 (calc from R=1.2)
    #beta = 0.08
    test_betas = [0.005, 0.015, 0.025, 0.035, 0.045, 0.055, 0.065]
    #beta = calculate_beta(R0, gamma, air_network)
    #beta = 0
    beta = 0.037 #gamma 0.5
    #jbeta = 0.005
    num_metro_zeros = 1 # set how many metros to select patients from to start with
    num_child_zeros = 1
    num_adult_zeros = 0
    time_end = 150
    #time_end = 7 # to test if contact matrix changes
    
    # DEFINE TRAVEL PARAMETERS
    theta_susc = 1
    theta_infc = 1
    theta_recv = 1
    
    # DEFINE EXP PARAMETERS
    #exp_param_list, exp_inter_list = [], []
    manip_exp_params = []
    experiment = 'yes' #turn exp on or off
    manip_exp_params.append(experiment)
    #exp_param_list.append(experiment)
    #time_list, exp_list = exp_func.list_exp_names()
    #intervention = exp_list[0] # identify which experiment from the list
    intervention = 'red_C_cc'
    manip_exp_params.append(intervention)
    #if experiment == 'yes':
    #    exp_inter_list.append(intervention)
    #timing = time_list[1] # 1 = time_steps
    timing = 'real_time'
    length_before = 7 # days
    length_after = 7 # days
    manip_exp_params.append(timing)
    manip_exp_params.append(length_before)
    manip_exp_params.append(length_after)
    #time_exp_start, time_exp_end = exp_func.set_time_start_and_end(timing, length_before, length_after)
    #exp_param_list.append(time_exp_start)
    #exp_param_list.append(time_exp_end)

    
    beta_exp = 0.02
    manip_exp_params.append(beta_exp)
    #exp_param_list.append(beta_exp)

    #test assigning parameters

    #time_exp_end = exp_func.set_time_end(timing)
    
    
    #
    #if intervention == 'red_C_cc':
    #    C_cc_red = exp_func.reduce_C_cc(C)
    #    exp_param_list.append(C_cc_red)

    #exp_param_list = [experiment, time_exp_start, time_exp_end, beta_exp, C_cc_red]
    
    # RUN EPIDEMIC SIMULATIONS
    average_epidemic_size, std_dev, adult_epi_size, child_epi_size, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult = chain_binomial_monte_carlo_plots_csv(R0, beta, gamma, alpha, theta_susc, theta_infc, theta_recv, time_end, num_metro_zeros, num_child_zeros, num_adult_zeros, d_metropop, metro_ids, filename_metropop, air_network, ch_travelers_r, ad_travelers_s, params_to_calc_C, intervention, manip_exp_params)

    #print incidence_time_series_metro_child
    
    # UNIT TESTS
    #plot beta values v. avg epi size
    #beta_test(test_betas, gamma, num_metro_zeros, num_child_zeros, num_adult_zeros)
    #gamma_test(beta, test_gammas, num_metro_zeros, num_child_zeros, num_adult_zeros)
    
    # OUTPUT RESULTS
    print "\nAverage Large Epidemic Size = ", round(100*average_epidemic_size,2), '%.\n'
    print "\nAverage Adult Epidemic Size = ", round(100*adult_epi_size,2), '%.\n'
    print "\nAverage Child Epidemic Size = ", round(100*child_epi_size,2), '%.\n'

    #use a chain binomial monte carlo function that has csv in the title, not because it writes the csv but bc it outputs the data necessary for this function to write csv
    write_csv_file_experiments(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, beta, intervention)
    #write_csv_file_test_beta(num_metro_zeros, num_child_zeros, num_adult_zeros, incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, beta)
    #write_csv_file(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult)
