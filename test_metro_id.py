metro_ids = [1, 2, 3]
x = [100, 150, 75]
metro_pops = dict(zip(metro_ids, x))

assign_metro = {}
last_index = ((len(metro_ids)))
for i in range(0, last_index): # range gives indexes
    start_id = sum(x[0:i])
    end_id = sum(x[0:i+1]) - 1
    for j in range(start_id, (end_id + 1)):
        pop = ((end_id - start_id) + 1) # might be an issue to use pop to grab metro_id if multiple have the same pop
        metro_id = [key for key in metro_pops if metro_pops[key] == pop]
        assign_metro[j] = metro_id[0]
        #metro_id_index = (i - 1)
        #metro_id = metro_ids[metro_id_index] 
        #assign_metro[j] = metro_id
        
print [key for key in assign_metro if assign_metro[key] == 3]

    
