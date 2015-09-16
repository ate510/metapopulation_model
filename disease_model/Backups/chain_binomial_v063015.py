import networkx as nx # to work with the population contact network
import random as rnd
import numpy as np
import matplotlib.pyplot as plt

###################################################
def infected_degree(node, network, infected_list):
# returns number of neighbors of (node) that are infected

    neighbors = network.neighbors(node)
    
    infected_neighbors = list(set(neighbors) & set(infected_list))
    
    return len(infected_neighbors)
    
###################################################    
def susc_infc_event(infected_degree, beta):
# returns true or false for infection event with probability (1-exp(-beta*infected_degree))

    return (rnd.random() < (1-np.exp(-beta*infected_degree)))
    
###################################################    
def infc_recv_event(gamma):
# returns true or false for recovery event with probability (gamma)
#   assumes exponential distribution for infectious period

    return (rnd.random() < gamma)
    
###################################################
def update_SI(Susc, Infc, num_infected, time_step, new_cases):
    
    for case in new_cases:
        Susc.remove(case)
        Infc.append(case)
        
    num_infected[time_step] = len(new_cases)
    
###################################################
def update_IR(Infc, Recv, new_recv):
    
    for case in new_recv:
        Infc.remove(case)
        Recv.append(case)
        

###################################################
def chain_binomial_one_simulation(beta, gamma, contact_network, num_patient_zeros):
# given a per contact, per time step tranmssion probability (beta) and
#  a per time step recovery probability (gamma), use the SIR chain binomial
#  model to simulate ONE outbreak on the population (provided by contact_network)
#  and return the total number of infected individuals (num_infected)

    Susc = contact_network.nodes() # List of susceptible nodes (everyone in the population)
    Infc = [] # List of infecteds/infectious (empty for now)
    Recv = [] # List of recovereds (empty for now)
    
    time_step = 0 # clock counter keeping track of current time
    num_newly_infected = {} # this dictionary will keep track of how many infected in current time step
    
    # infect patient_zeros and upated Susc and Infc lists
    patient_zeros = rnd.sample(contact_network.nodes(), num_patient_zeros)
    update_SI(Susc, Infc, num_newly_infected, time_step, patient_zeros)
    
    while Infc: # while there are infectious individuals in the list
                
        time_step += 1 # update the clock
        new_cases = []
        new_recov = []

        # for every susceptible indivdiual in population check if they get infected in this time step
        for s in Susc:
            infected_degree_of_s = infected_degree(s, contact_network, Infc)
            if susc_infc_event(infected_degree_of_s, beta):
                new_cases.append(s)
            
        # for every infected indivdiuals in population check if they get infected in this time step
        for i in Infc:
            if infc_recv_event(gamma):
                new_recov.append(i)
            
        update_SI(Susc, Infc, num_newly_infected, time_step, new_cases)
        update_IR(Infc, Recv, new_recov)
        
    # Note num_newly_infected is the incidence time series      
    return num_newly_infected, sum(num_newly_infected.values()) # return total number infected in outbreak

###################################################
def chain_binomial_monte_carlo(beta, gamma, contact_network, num_sims, num_patient_zeros):
# run many (num_sims) instances of chain binomial simulation and return average epidemic size    

    population_size = contact_network.number_of_nodes()
    threshold = 0.10 # 10% of population size is our threshold for a large epidemic

    large_epidemic_sizes= [] # will keep list of outbreak sizes that are large epidemics
    for sim in range(1,num_sims):
        incidence_time_series, outbreak_size = chain_binomial_one_simulation(beta, gamma, contact_network, num_patient_zeros)
        # Note we are not using the incidence time series right now
        
        # figure out if this is small outbreak or large epidemic
        if outbreak_size > threshold * population_size: # if outbreak reached more than 10% of the population
            
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
def read_contact_network(filename):
# read population contact network from edgelist contained in file (filename)

    G = nx.Graph()
    G = nx.read_edgelist(filename, delimiter = '\t', data=True)
    
    return G
    
###################################################
def calc_avg_excess_degree(contact_network):
    deg_list = contact_network.degree().values()
    deg_list_minus_one = [d-1 for d in deg_list]
    deg_times_deg_minus_one = [d*d1 for d,d1 in zip(deg_list, deg_list_minus_one)]
    numerator = np.mean(deg_times_deg_minus_one)
    denominator = np.mean(deg_list)
        
    return float(numerator)/denominator

###################################################
def calculate_beta(R0, gamma, contact_network):
    
    # calculate T value for given R0
    #  where R0 = T * avg_excees_degree 
    avg_excess_deg = calc_avg_excess_degree(contact_network)
    T = float(R0)/avg_excess_deg
            
    # calculate beta based on T = beta/(beta + gamma)
    beta = T*gamma/(1-T)
        
    return beta
    
###################################################
if __name__ == "__main__":
    
    # READ POPULATION NETWORK FROM FILE
    filename_contact_network = '../air_traffic_data/air_traffic_edgelist.txt'
    contact_network = read_contact_network(filename_contact_network)
            
    # DEFINE DISEASE PARAMETERS
    R0 = 1.2
    gamma = 0.2 # recovery rate based on (1/gamma) day infectious period
    beta = calculate_beta(R0, gamma, contact_network)
    num_patient_zeros = 1 # set how many patient zeros to start with
    
    # RUN EPIDEMIC SIMULATIONS
    num_sims = 250 # if debugging, reduce this number to something small like 10
    average_epidemic_size = chain_binomial_monte_carlo(beta, gamma, contact_network, num_sims, num_patient_zeros)
    
    # OUTPUT RESULTS
    print "\nAverage Large Epidemic Size = ", round(100*average_epidemic_size,2), '%.\n'