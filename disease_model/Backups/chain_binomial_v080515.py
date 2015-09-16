import networkx as nx # to work with the population contact network
import random as rnd
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

### local modules ###
sys.path.append('/home/anne/Dropbox/Anne_Bansal_Lab')

### functions ###
import population_parameters as pop_func
    

###################################################    
def infected_contact_child(C, Infc_A, Infc_C):
# returns infected degree based on contact probabilities for adults

    C_cc = C.item((0, 0))
    C_ca = C.item((0, 1))
    
    return ((C_cc*Infc_C)+(C_ca*Infc_A))
        
###################################################    
def infected_contact_adult(C, Infc_A, Infc_C):
# returns infected degree based on contact probabilities for children

    C_ac = C.item((1, 0))
    C_aa = C.item((1, 1))

    return ((C_ac*Infc_C)+(C_aa*Infc_A))
    
###################################################    
def susc_infc_event_child(infected_contact_child, beta):
# returns true or false for infection event with probability (1-exp(-beta*infected_degree))

    return (rnd.random() < (1-np.exp(-beta*infected_contact_child)))
    
###################################################    
def susc_infc_event_adult(infected_contact_adult, beta):
# returns true or false for infection event with probability (1-exp(-beta*infected_degree))

    return (rnd.random() < (1-np.exp(-beta*infected_contact_adult)))
    
###################################################
def susc_infc_binomial_prob_child (met_id, d_Infc, C, beta):
# returns probability of infection for binomial function    
    # Si --> Ii
    # determine how many susceptibles get infected in each metro area
  
    Infc_C = d_Infc[(met_id, 'C')] #number children infected in that metro area
    Infc_A = d_Infc[(met_id, 'A')] #number adults infected in that metro area
    infected_degree_of_s_C = infected_contact_child(C, Infc_A, Infc_C) #number of infected contacts
    
    return (susc_infc_event_child(infected_degree_of_s_C, beta))
    
    #if susc_infc_event_child(infected_degree_of_s_C, beta):
    #    new_cases = new_cases + 1
        
###################################################
def susc_infc_binomial_prob_adult (met_id, d_Infc, C, beta):
# returns probability of infection for binomial function    

    # Si --> Ii
    # determine how many susceptibles get infected in each metro area

    Infc_C = d_Infc[(met_id, 'C')]
    Infc_A = d_Infc[(met_id, 'A')]  
    infected_degree_of_s_A = infected_contact_adult(C, Infc_A, Infc_C)
    
    return (susc_infc_event_adult(infected_degree_of_s_A, beta))
 
              
###################################################    
def infc_recv_event(gamma):
# returns true or false for recovery event with probability (gamma)
#   assumes exponential distribution for infectious period
# same for children and adults

    return (rnd.random() < gamma)

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
def update_SI(met_id, d_Susc, d_Infc, num_new_infc_child, num_new_infc_adult):
# updates dictionary values for d_Susc, d_Infc
    
    #children
    d_Susc[(met_id, 'C')] -= num_new_infc_child
    d_Infc[(met_id, 'C')] += num_new_infc_child
    
    #adults
    d_Susc[(met_id, 'A')] -= num_new_infc_adult
    d_Infc[(met_id, 'A')] += num_new_infc_adult
    
###################################################
def update_IR(met_id, d_Infc, d_Recv, new_recov_child, new_recov_adult, d_metro_infected_child, d_metro_infected_adult, time_step):
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
    d_metro_infected_child[(met_id, time_step)] = num_infc_child
    d_metro_infected_adult[(met_id, time_step)] = num_infc_adult
    #d_metro_tot_infected_child[(met_id)] = 
    
###################################################
def travel_btwn_metros(air_network, d_Susc, d_Infc, d_Recv, d_prob_travel_C, d_prob_travel_A, theta_susc, theta_infc, theta_recv):
    
    edges = air_network.edges() 
    for (i, j) in edges:
        # Si <-> Sj
        
        # children
        ch_travel_i_j = ((d_prob_travel_C[(i, j)]) * (d_Susc[(i, 'C')]) * (theta_susc)) # number of susceptible children traveling from i to j
        ch_travel_j_i = ((d_prob_travel_C[(j, i)]) * (d_Susc[(j, 'C')]) * (theta_susc)) # number of susceptible children traveling from j to i
        net_travel = ((ch_travel_j_i) - (ch_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)
        # update pop sizes
        d_Susc[(i, 'C')] = ((d_Susc[(i, 'C')]) + net_travel)
        d_Susc[(j, 'C')] = ((d_Susc[(j, 'C')]) - net_travel)
	        
        #adults
        ad_travel_i_j = ((d_prob_travel_A[(i, j)]) * (d_Susc[(i, 'A')]) * (theta_susc)) # number of susceptible adults traveling from i to j
        ad_travel_j_i = ((d_prob_travel_A[(j, i)]) * (d_Susc[(j, 'A')]) * (theta_susc)) # number of susceptible adults traveling from j to i
        net_travel = ((ad_travel_j_i) - (ad_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)
        # update pop sizes	        
        d_Susc[(i, 'A')] = ((d_Susc[(i, 'A')]) + net_travel)
        d_Susc[(j, 'A')] = ((d_Susc[(j, 'A')]) - net_travel)
               
               
        # Ii <-> Ij
    
        # children
        ch_travel_i_j = ((d_prob_travel_C[(i, j)]) * (d_Infc[(i, 'C')]) * (theta_infc)) # number of infected children traveling from i to j
        ch_travel_j_i = ((d_prob_travel_C[(j, i)]) * (d_Infc[(j, 'C')]) * (theta_infc)) # number of infected children traveling from j to i
        net_travel = ((ch_travel_j_i) - (ch_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)
        # update pop sizes
        d_Infc[(i, 'C')] = ((d_Infc[(i, 'C')]) + net_travel)
        d_Infc[(j, 'C')] = ((d_Infc[(j, 'C')]) - net_travel)
	        
        #adults
        ad_travel_i_j = ((d_prob_travel_A[(i, j)]) * (d_Infc[(i, 'A')]) * (theta_infc)) # number of infected adults traveling from i to j
        ad_travel_j_i = ((d_prob_travel_A[(j, i)]) * (d_Infc[(j, 'A')]) * (theta_infc)) # number of infected adults traveling from j to i
        net_travel = ((ad_travel_j_i) - (ad_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)
        # update pop sizes	        
        d_Infc[(i, 'A')] = ((d_Infc[(i, 'A')]) + net_travel)
        d_Infc[(j, 'A')] = ((d_Infc[(j, 'A')]) - net_travel)


        # Ri <-> Rj
    
        # children
        ch_travel_i_j = ((d_prob_travel_C[(i, j)]) * (d_Recv[(i, 'C')]) * (theta_recv)) # number of recovered children traveling from i to j
        ch_travel_j_i = ((d_prob_travel_C[(j, i)]) * (d_Recv[(j, 'C')]) * (theta_recv)) # number of recovered children traveling from j to i
        net_travel = ((ch_travel_j_i) - (ch_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)
        # update pop sizes
        d_Recv[(i, 'C')] = ((d_Recv[(i, 'C')]) + net_travel)
        d_Recv[(j, 'C')] = ((d_Recv[(j, 'C')]) - net_travel)
	        
        #adults
        ad_travel_i_j = ((d_prob_travel_A[(i, j)]) * (d_Recv[(i, 'A')]) * (theta_recv)) # number of recovered adults traveling from i to j
        ad_travel_j_i = ((d_prob_travel_A[(j, i)]) * (d_Recv[(j, 'A')]) * (theta_recv)) # number of recovered adults traveling from j to i
        net_travel = ((ad_travel_j_i) - (ad_travel_i_j)) # difference - if (+), j pop decreases, and i pop increases, vice versa for (-)
        # update pop sizes	        
        d_Recv[(i, 'A')] = ((d_Recv[(i, 'A')]) + net_travel)
        d_Recv[(j, 'A')] = ((d_Recv[(j, 'A')]) - net_travel)
        
        
###################################################
def chain_binomial_one_simulation(d_metro_age_pop, metro_ids, beta, gamma, air_network, num_metro_zeros, num_child_zeros, num_adult_zeros, C):
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
    d_metro_infected_child, d_metro_infected_adult = {}, {}
    #record number total infected at each metro area
    d_metro_tot_infected_child, d_metro_tot_infected_adult = {}, {}
    # record only new cases at each time step and each metro
    d_new_cases_child, d_new_cases_adult = {}, {}
    # record new cases for each time step (agg over metro areas)
    d_tot_new_cases_child, d_tot_new_cases_adult = {}, {}
    
    # infect patient_zeros and upated Susc and Infc lists
    metro_zeros = rnd.sample(metro_ids, num_metro_zeros)
    #second patient_zero selection for indv - select fixed number of patient_zeros - #children, #adults
    for met_id in metro_zeros:
        update_SI(met_id, d_Susc, d_Infc, num_child_zeros, num_adult_zeros)
        num_infc_child = d_Infc[(met_id, 'C')]
        num_infc_adult = d_Infc[(met_id, 'A')]
        d_metro_infected_child[(met_id, time_step)] = num_infc_child
        d_metro_infected_adult[(met_id, time_step)] = num_infc_adult
        d_metro_tot_infected_child[(met_id, time_step)] = num_infc_child
        d_metro_tot_infected_adult[(met_id, time_step)] = num_infc_adult
        d_new_cases_child[(met_id, time_step)] = num_infc_child
        d_new_cases_adult[(met_id, time_step)] = num_infc_adult
       
    metros_not_zeros = [metro for metro in metro_ids if metro not in metro_zeros] 
    for met_id in metros_not_zeros:
        d_metro_infected_child[(met_id, time_step)] = 0
        d_metro_infected_adult[(met_id, time_step)] = 0
        d_metro_tot_infected_child[(met_id, time_step)] = 0
        d_metro_tot_infected_adult[(met_id, time_step)] = 0
        d_new_cases_child[(met_id, time_step)] = 0
        d_new_cases_adult[(met_id, time_step)] = 0
    
    num_infc_child = sum([d_Infc[(met_id, 'C')] for met_id in metro_ids])
    num_infc_adult = sum([d_Infc[(met_id, 'A')] for met_id in metro_ids])
    d_nat_infected_child[(time_step)] = num_infc_child
    d_nat_infected_adult[(time_step)] = num_infc_adult
    
    # while there are infected individuals
    # go to next time step
    num_infected = ((d_nat_infected_child[time_step]) + (d_nat_infected_adult[time_step])) 
    while num_infected > 0 and time_step < time_end:
        
        time_step += 1 #update clock
    
        # TRAVEL #
    
        # create two dictionaries with probabilities of travel for each age group, keys being tuples of cities: (i, j) and (j, i)
        d_prob_travel_C, d_prob_travel_A = pop_func.calc_prob_travel(air_network, alpha, ch_travelers_r, d_metropop)
        
        #update population sizes for S, I, R for each metro
        travel_btwn_metros(air_network, d_Susc, d_Infc, d_Recv, d_prob_travel_C, d_prob_travel_A, theta_susc, theta_infc, theta_recv)
        
        # DISEASE #
    
        for met_id in metro_ids:
            # Si --> Ii
            # determine how many susceptibles get infected in each metro area
            
            # child
            Susc_C = d_Susc[(met_id, 'C')] # number of susc children in metro area
            prob = susc_infc_binomial_prob_child(met_id, d_Infc, C, beta) # calc probability of infection
            # determine how many are infected (coin flip 'Susc_C' number of times, with probability 'prob' each flip will result in infected)
            new_cases_child = np.random.binomial(Susc_C, prob) # determine how many are infected
            d_new_cases_child[(met_id, time_step)] = new_cases_child
            previous_time_steps = range(0, time_step)
            previous_cases = sum([d_metro_tot_infected_child[(met_id, time)] for time in previous_time_steps])
            d_metro_tot_infected_child[(met_id, time_step)] = previous_cases + new_cases_child
                         
            #adult
            Susc_A = d_Susc[(met_id, 'A')] # number of susc adults in metro area
            prob = susc_infc_binomial_prob_adult(met_id, d_Infc, C, beta) # calc probability of infection
            new_cases_adult = np.random.binomial(Susc_A, prob) 
            d_new_cases_adult[(met_id, time_step)] = new_cases_adult
            previous_time_steps = range(0, time_step)
            previous_cases = sum([d_metro_tot_infected_adult[(met_id, time)] for time in previous_time_steps])
            d_metro_tot_infected_adult[(met_id, time_step)] = previous_cases + new_cases_adult
                    
            #subtract from Si, add to Ii
            update_SI(met_id, d_Susc, d_Infc, new_cases_child, new_cases_adult)

            # Ii --> Ri
            # determine how many child / adult infected will recover in each metro area
                
            #child
            Infc_C = d_Infc[(met_id, 'C')] # number of infected children in metro area
            prob = gamma # probability of recovery = gamma
            new_recov_child = np.random.binomial(Infc_C, prob)
            
            #adult
            Infc_A = d_Infc[(met_id, 'A')]
            prob = gamma
            new_recov_adult = np.random.binomial(Infc_A, prob)
            
            # subtract from Ii, add to Ri            
            update_IR(met_id, d_Infc, d_Recv, new_recov_child, new_recov_adult, d_metro_infected_child, d_metro_infected_adult, time_step)
        
        #record how many total infected across metro_ids at this time step
        num_infc_child = sum([d_Infc[(met_id, 'C')] for met_id in metro_ids])
        num_infc_adult = sum([d_Infc[(met_id, 'A')] for met_id in metro_ids])
        d_nat_infected_child[(time_step)] = num_infc_child
        d_nat_infected_adult[(time_step)] = num_infc_adult
        d_tot_new_cases_adult[(time_step)] = sum([d_new_cases_adult[(met_id, time_step)] for met_id in metro_ids])
        d_tot_new_cases_child[(time_step)] = sum([d_new_cases_child[(met_id, time_step)] for met_id in metro_ids])        
        print d_tot_new_cases_child[(time_step)]
        # go back thru while loop
        # next time step
        # travel again
        # S --> I --> R
        
    # Note num_newly_infected is the incidence time series      
    return d_metro_infected_child, d_metro_infected_adult, d_metro_tot_infected_child, d_metro_tot_infected_adult, sum(d_tot_new_cases_child.values()), sum(d_tot_new_cases_adult.values()) # return total number infected in outbreak

###################################################
def chain_binomial_monte_carlo(beta, gamma, num_sims, num_metro_zeros, num_child_zeros, num_adult_zeros):
# run many (num_sims) instances of chain binomial simulation and return average epidemic size    

    #d_metropop, metro_ids = pop_func.import_metropop(filename_metropop, 2, 3)
    population_size = sum([d_metropop[x] for x in metro_ids])
    threshold = 0.10 # 10% of population size is our threshold for a large epidemic

    d_metro_age_pop = pop_func.calc_metro_age_pop(filename_metropop, alpha)

    large_epidemic_sizes= [] # will keep list of outbreak sizes that are large epidemics
    for sim in range(1,num_sims):
        #incidence_time_series, outbreak_size = chain_binomial_one_simulation(d_metro_age_pop, metro_ids, beta, gamma, air_network, num_metro_zeros, num_child_zeros, num_adult_zeros, C)
        # Note we are not using the incidence time series right now
        
        incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, outbreak_size_child, outbreak_size_adult = chain_binomial_one_simulation(d_metro_age_pop, metro_ids, beta, gamma, air_network, num_metro_zeros, num_child_zeros, num_adult_zeros, C)

        outbreak_size = (outbreak_size_child + outbreak_size_adult)
        
        ## SB said don't worry about this for now ##
        # figure out if this is small outbreak or large epidemic
        #if outbreak_size > threshold * population_size: # if outbreak reached more than 10% of the population
            #large_epidemic_sizes.append(outbreak_size)

        # call it a large epidemic and save its size
        large_epidemic_sizes.append(outbreak_size)
            
    # calculate average large epidemic size, and how frequent they were
    if large_epidemic_sizes:
        average_epidemic_size = np.mean(large_epidemic_sizes)/float(population_size)
        probability_epidemic = len(large_epidemic_sizes)/float(num_sims)
    else:
        average_epidemic_size = 0
        probability_epidemic = 0
    
    
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

###################################################
def write_csv_file(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult):
    
    #time_series = range(0, time_end)
    csvfile = csv.writer(open('/home/anne/Dropbox/Anne_Bansal_lab/Modeling_Project_Outputs/chain_binomial_output_nummetrozeros_%_numchildzeros_%_numadultzeros_%_numsims_%.csv' % (num_metro_zeros, num_child_zeros, num_adult_zeros, num_sims), 'wb'), delimiter = ',')
    csvfile.writerow(['time_step', 'metro_id', 'age', 'currently_infected', 'total_infected'])
    for (met_id, time_step) in incidence_time_series_metro_child:
        csvfile.writerow([time_step, met_id, 'C', (incidence_time_series_metro_child[(met_id, time_step)]), (tot_incidence_time_series_child[(met_id, time_step)])])
    for (met_id, time_step) in incidence_time_series_metro_adult:
        csvfile.writerow([time_step, met_id, 'A', (incidence_time_series_metro_adult[(met_id, time_step)]), (tot_incidence_time_series_adult[(met_id, time_step)])])


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
    filename_germ_contact_data = 'Dropbox/Anne_Bansal_lab/Contact_Data/polymod_germany_contact_matrix_Mossong_2008.csv'
    filename_germ_pop_data = 'Dropbox/Anne_Bansal_lab/UNdata_Export_2008_Germany_Population.csv'
        
    # DEFINE POPULATION PARAMETERS
    year = 2010
    alpha = pop_func.calc_alpha(year, dict_childpop, dict_adultpop)
    ch_travelers_r = 0.0 # fraction of children who travel
    
    # CONTACT MATRIX
    C = pop_func.calc_contact_matrix(filename_germ_contact_data, filename_germ_pop_data, alpha)
                            
    # DEFINE DISEASE PARAMETERS
    R0 = 1.2
    gamma = 0.2 # recovery rate based on (1/gamma) day infectious period
    #beta = calculate_beta(R0, gamma, air_network)
    beta = 0.015
    #beta = 0.005
    num_metro_zeros = 1 # set how many metros to select patients from to start with
    num_child_zeros = 1
    num_adult_zeros = 0
    time_end = 300
    
    # DEFINE TRAVEL PARAMETERS
    theta_susc = 1
    theta_infc = 1
    theta_recv = 1
    
    # RUN EPIDEMIC SIMULATIONS
    #num_sims = 250 # if debugging, reduce this number to something small like 10
    num_sims = 10
    average_epidemic_size = chain_binomial_monte_carlo(beta, gamma, num_sims, num_metro_zeros, num_child_zeros, num_adult_zeros)
    
    # OUTPUT RESULTS
    print "\nAverage Large Epidemic Size = ", round(100*average_epidemic_size,2), '%.\n'
    incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult, outbreak_size_child, outbreak_size_adult = chain_binomial_one_simulation(d_metro_age_pop, metro_ids, beta, gamma, air_network, num_metro_zeros, num_child_zeros, num_adult_zeros, C)	
    write_csv_file(incidence_time_series_metro_child, incidence_time_series_metro_adult, tot_incidence_time_series_child, tot_incidence_time_series_adult)