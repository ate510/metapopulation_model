C = 0
A = 0

num_child_zeros = 1
num_adult_zeros = 1

dict_metro_infc = {}
dict_metro_infc[(1, 'A')] = 0
dict_metro_infc[(1, 'C')] = 0
dict_metro_infc[(2, 'A')] = 0
dict_metro_infc[(2, 'C')] = 0

dict_metro_susc = {}
dict_metro_susc[(1, 'A')] = 1
dict_metro_susc[(1, 'C')] = 1
dict_metro_susc[(2, 'A')] = 1
dict_metro_susc[(2, 'C')] = 1

metro_ids = [1, 2]

for met_id in metro_ids:
    #children
    dict_metro_susc[(met_id, 'C')] = ((dict_metro_susc[(met_id, 'C')]) - num_child_zeros)
    dict_metro_infc[(met_id, 'C')] = ((dict_metro_infc[(met_id, 'C')]) + num_child_zeros)
    C = (C + (dict_metro_infc[(met_id, 'C')]))
    #adults
    dict_metro_susc[(met_id, 'A')] = ((dict_metro_susc[(met_id, 'A')]) - num_adult_zeros)
    dict_metro_infc[(met_id, 'A')] = ((dict_metro_infc[(met_id, 'A')]) + num_adult_zeros)
    A = (A + (dict_metro_infc[(met_id, 'A')]))
    