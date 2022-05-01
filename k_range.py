from .methods import sum_method
import numpy as np
from .ratio import ratio_matrix
import scipy as scp



#Function to understand the range of suitable k for sum method, then it has to verify also the suitable 
#pc values, after that it we have to compute the pc for ehich it is most parallel


def method_range(final_k, densities, critical_pcs, residues, dataframe, cutoffs, method=sum_method):
    final_k_list = []
    for k_init, residue, pc, cutoff in zip(final_k, residues, critical_pcs, cutoffs):
        #The result of the evaluation with the method
        difference = 0
        final_k_list= []
        #The central k 
        final_k_list.append(k_init)
        initial_k = k_init + 1
        dens = np.abs(np.array(densities) - pc)
        while (initial_k < 280):
            rat_array = ratio_matrix([initial_k], densities, dataframe)
            results = scp.stats.linregress(np.log(dens), np.log(rat_array[0,:]) )
            line = results.intercept + results.slope*np.log(dens)
                

            residues_new = np.log(rat_array[0,:])[0:] - (line)
            new_method = method(residues_new)
            difference = abs(new_method - residue)
            #print(initial_k, new_method, difference, residue)
            if abs(difference) < cutoff:
                final_k_list.append(initial_k)
            

            #updating the k
            initial_k += 1
        
        initial_k = k_init - 1

        #reset differnece
        difference = 0
        while  (initial_k > -1):
            rat_array = ratio_matrix([initial_k], densities, dataframe)
            results = scp.stats.linregress(np.log(dens), np.log(rat_array[0,:]) )
            
            
            
            line = results.intercept + results.slope*np.log(dens)
            residues_new = np.log(rat_array[0,:])[0:] - (line)
            new_method = method(residues_new)
            difference = abs(new_method - residue)
            if abs(difference) < cutoff:
                final_k_list.append(initial_k)

            #updating the k
            initial_k -= 1

    

    return(final_k_list) 
            



def info_method_range():
    print('''
    Function to understand the range of suitable k for the custom method for each critica p provided 
    in input

    Takes in input:
    1) final_k: the array of ks to test
    2) densities: the array of densities
    3) critical_pcs: the array of critical pcs for each k
    4) residues: an array of values, one for every k, that represent the method value
    5) dataframe: the dataframe from which to retrieve the size (the dataframe has to have the format
    specified in info_dataframe)
    6) cutoffs: if the value returned by methods is distant more the cutoffs from the reference value
    it gets discarded
    7) methods=sum_method: The custom method acting on the residues to evaluate linearity
    
    
    
    
    
    
    ''')