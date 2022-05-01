import numpy as np
from .ratio import ratio_matrix
import scipy as scp
import pandas as pd
from .methods import min_method, sum_method
from .ratio import ratio_matrix


def best_pc(densities, dataframe, chunk_length=5, k_space=1, pc_space=0.5):
    #densities has to be a numpy array
    pc_max = np.max(densities)
    
    #The pcs to test
    pc_test = np.arange(0, pc_max, pc_space)
    k_array = np.arange(0, 101, k_space)
    #The length of each chunk ok k_values
    chunk_length = chunk_length
    
    
    #The dictionary with all the std of the phi of each pc
    std_dict = {}
    
    #The dictionary with all the means of the phi of each pc
    mean_dict = {}
    
    #Dictionary with the final k for every pc
    final_k_dict = {}
    
    #Dictionary for the mean r2 for every pc
    final_r2 = {}
    
    #Creating a sub-dictionary with all the possible pc_positions
    for pc_pos in np.arange(len(densities)):
        std_dict[pc_pos] = {}
        mean_dict[pc_pos] = {}
    
    for pc in pc_test:
        
        #Getting how many densities are below pc
        pc_position = 0
        for density in densities:
            if density < pc:
                pc_position += 1
        
        
        #The densities to try
        dens = densities - pc
        
        #The dictionary for the values of the means of r^2
        r2_means = {}
        
        for i in np.arange(int((np.max(k_array))/chunk_length)):
            #The ks to be sent in input to the ratio_matrix function
            k_ratio = k_array[i:chunk_length + i]
            
            #The matrix of ratios
            #here we pass as densities denstoes and nt dens because the ratios matrix
            #is independent from the densities and we need the original densities to get the delta p, delta s...
            
            rat_matrix = ratio_matrix(k_ratio, densities, dataframe)
            
            #alterntive depending on the import
            #rat_matrix = ap.ratio.ratio_matrix(k_ratio, densities, dataframe)
            
            #The list of the r^2 values for each chunk
            r2_list = []
            
            #computing the r^2 values
            for row in np.arange(rat_matrix.shape[0]):
                #We start from pc postion in y to have the same length between x and y
                results = scp.stats.linregress(x = np.log(dens), y=np.log(rat_matrix[row, :]))
                r2_list.append(results.rvalue)
                
            r2_means[i] = np.mean(r2_list)
            
        #getting the chunk with the ks that have the highest mean of k_values    
        max_chunk = max(r2_means, key=r2_means.get)
        max_r2 = r2_means[max_chunk]
        
        final_r2[pc] = max_r2
        
        #The final k_value chosen for this pc
        final_k = k_array[max_chunk:  max_chunk + chunk_length]
        
        final_k_dict[pc] = final_k
        
        #The ratio matrix for the values found for k used to estimate phi
        ratio_phi = ratio_matrix(final_k, densities, dataframe)
        
        #Alternative based on the input
        #ratio_phi = ap.ratio.ratio_matrix(final_k, densities, dataframe)
        
        #The list of the phi for a given pc
        phi_list = []
        
        
        for row in np.arange(ratio_phi.shape[0]):    
            #We start from pc postion in y to have the same length between x and y
            results_phi = scp.stats.linregress(x = np.log(dens), y=np.log(ratio_phi[row, :]))
            phi_list.append(results_phi.slope)
            
            
        std_dict[pc_position][pc] = np.std(phi_list)
        mean_dict[pc_position][pc] = np.mean(phi_list)
        
    
   #The dictionaty with the phi of the minimum std  
    min_std_pc = {}
    
    
    min_mean_pc = {}
    
    
    for position in np.arange(len(densities)):
        min_std_pc[position] = {}
        
        
        min_pc = min(std_dict[position], key=std_dict[position].get)
        
        
        min_std = np.min(list(std_dict[position].values()))
        
        k_min_std = final_k_dict[min_pc]
        
        max_r2_pc = final_r2[min_pc]
        
        min_std_pc[position] = [min_pc, min_std, k_min_std, max_r2_pc]
        min_mean_pc[position] = mean_dict[position][min_pc]
    
    
    return_1 = pd.DataFrame.from_dict(min_std_pc, orient='index', columns=['pc', 'std', 'k', 'max r2'])
    return_2 = pd.DataFrame.from_dict(min_mean_pc, orient='index', columns = ['phi'])
    
    return return_1, return_2


def info_best_pc():
    print('''
    DEPRECATED
    Takes in input,
    1) densities: an array of densities,
    2) dataframe: a dataframe from which to read the sizes (the dataframe has to have the format 
    3) chunk_length = 5: the length of each chunk k that gets tested
    4) k_space = 1, the step used to move through the ks to test
    5) pc_space = 0.5, the step used to move through the pcs to test
    defined in the info_dataframe function of this method)
    This function searches for the range of k for which the (k+1)-th moment/k-th moment plotted
    as a function of the densities is straigther (in log-log), then it searches the pc that makes the lines
    more parallel (plotting abs(p - pc)) in the x axis.
    It is the version with pc positions still present and thus is deprecated

    , k_space=1, pc_space=0.5
    ''')



def linearity_pc(densities, dataframe, k_space=1, pc_space=0.5, method=min_method):
    #densities has to be a numpy array
    #we exclude the last one because the analysis on the last one does not make sense
    pc_max = np.max(densities[:-1])
    
    #The pcs to test
    pc_test = np.arange(0, pc_max, pc_space)

    #Here we compose the k_array defined as the concatenation of 2 k_tmp
    k_tmp_1 = np.linspace(0, 1, 11)[:-1]
    k_tmp_2 = np.arange(1, 250, k_space) 
    k_array = np.concatenate((k_tmp_1, k_tmp_2), axis=0)
    
    #met_dict stands for method_dict: the result of the chosen method for a given k will be stored here
    met_dict = {}
    
    #The pcs giving the best result of the method 
    best_pc_dict = {}

    #
    best_met_dict = {}
    
    

    for k in k_array:
        
        
        #for pc_pos in np.arange(len(densities[:-1])):
        met_dict[k] = {}
        best_pc_dict[k] = {}
        best_met_dict[k] = {}

        for pc in pc_test:

            #pc_position = 0
            #for density in densities:
            #    if density < pc:
            #        pc_position += 1
            
            
            dens = np.abs(np.array(densities) - pc)
            
            rat_array = ratio_matrix([k], densities, dataframe)
            
            
            for row in np.arange(rat_array.shape[0]): 
                
                #print(np.log(rat_array[row,:])[pc_position:])
            #We start from pc postion in y to have the same length between x and y
                results = scp.stats.linregress(np.log(dens), np.log(rat_array[row,:]) )
                
                line = results.intercept + results.slope*np.log(dens)
                
                residues = np.log(rat_array[row,:]) - (line)  
                
                met_dict[k][pc] = method(residues)
                
        
        #for pc_pos in np.arange(len(densities[:-1])):
            #best_dens_dict[k][pc_pos] = np.min(met_dict[(k, pc_pos)].values())[0]
                
        best_pc_dict[k] = min(met_dict[k], key=met_dict[k].get)
            
        
        
        best_met_dict[k] = np.min(list((met_dict[k].values())))
            
    
    
    res = pd.DataFrame.from_dict(best_met_dict, orient='Index', columns=['residues'])
    pc_ = pd.DataFrame.from_dict(best_pc_dict, orient='Index', columns=['pc'])
    
    res['k'] = res.index
    
    
    

    pc_['k'] = res.index
    

    res.reset_index()
    #tmp_res = res.groupby(by='pc_pos')['residues'].idxmin()
    
    #final_res = res.loc[tmp_res]
    
    #final_df = pd.merge(final_res, pc_, on='(k, pc_pos)')

    #return1 = pd.DataFrame.from_dict(best_set, orient='index', columns=['best set'])
    #return(best_met_dict, best_pc_dict)
    return(res, pc_)   

              

def info_linearity_pc():
    print('''
    This function is is similar to best_pc and has the same goal, it does it though with a custom method
    acting on the residues and returns results with a better visualization and above all it doesn't act
    with pc_position
    Takes in input,
    1) densities: an array of densities,
    2) dataframe: a dataframe from which to read the sizes (the dataframe has to have the format 
    3) chunk_length = 5: the length of each chunk k that gets tested
    4) pc_space = 0.5, the step used to move through the pcs to test
    5) method = methods.min_method: the method used for the evaluation, acting on the residues
    ''')




#Function to determine the range of the pc, eventually to insert in the above function

def method_range_pc(final_k_dict, densities, critical_pcs, residues, dataframe,  cutoffs, method=sum_method):
    final_pc_list_max = []
    final_pc_list_min = []
   
    #for pc_ in  pc_positions:
     #   max_k_dict[pc_] = np.max(final_k_dict[pc_])
     #   min_k_dict[pc_] = np.min(final_k_dict[pc_])


    for residue, critical_pc, cutoff in zip(residues, critical_pcs, cutoffs):
        
        difference = 0
        #The result of the evaluation with the method
        #final_pc_list_max[pc_pos] = []
        #final_pc_list_min[pc_pos] = []
        
        #The central pc 
        initial_pc = critical_pc

        max_k = np.max(final_k_dict)
        #print(max_k)
        while  (initial_pc < np.max(densities)):
            dens = np.abs(np.array(densities) - initial_pc)
            rat_array = ratio_matrix([max_k], densities, dataframe)
            results = scp.stats.linregress(np.log(dens), np.log(rat_array[0,:])[:] )
            line = results.intercept + results.slope*np.log(dens)
                

            residues_new = np.log(rat_array[0,:])[0:] - (line)
            new_method = method(residues_new)
            difference = abs(new_method - residue)
            #print(pc, difference)
            #print(initial_k, new_method, difference, residue)
            if abs(difference) < cutoff:
                final_pc_list_max.append(initial_pc)
            

            #updating the k
            initial_pc += 0.1
        
        initial_pc = critical_pc - 0.1

        #reset differnece
        difference = 0
        while  (initial_pc > 0):
            dens = np.abs(np.array(densities) - initial_pc)
            rat_array = ratio_matrix([max_k], densities, dataframe)
            results = scp.stats.linregress(np.log(dens), np.log(rat_array[0,:])[:] )
            line = results.intercept + results.slope*np.log(dens)
                

            residues_new = np.log(rat_array[0,:])[0:] - (line)
            new_method = method(residues_new)
            difference = abs(new_method - residue)
            #print(initial_k, new_method, difference, residue)
            if (abs(difference) < cutoff):
                final_pc_list_max.append(initial_pc)
            

            #updating the k
            initial_pc -= 0.1
        

        initial_pc = critical_pc
        
        min_k = np.min(final_k_dict)
        
        while  (initial_pc < 33.33):
            dens = np.abs(np.array(densities) - initial_pc)
            rat_array = ratio_matrix([min_k], densities, dataframe)
            results = scp.stats.linregress(np.log(dens), np.log(rat_array[0,:])[:] )
            line = results.intercept + results.slope*np.log(dens)
                

            residues_new = np.log(rat_array[0,:])[0:] - (line)
            new_method = method(residues_new)
            difference = abs(new_method - residue)
            #print(initial_k, new_method, difference, residue)
            if (abs(difference) < cutoff):    
                final_pc_list_min.append(initial_pc)
            

            #updating the k
            initial_pc += 0.1
        
        initial_pc = critical_pc - 0.1

        #reset differnece
        difference = 0
        while (initial_pc > 0):
            dens = np.abs(np.array(densities) - initial_pc)
            rat_array = ratio_matrix([min_k], densities, dataframe)
            results = scp.stats.linregress(np.log(dens), np.log(rat_array[0,:]) )
            line = results.intercept + results.slope*np.log(dens)
                

            residues_new = np.log(rat_array[0,:])[0:] - (line)
            new_method = method(residues_new)
            difference = abs(new_method - residue)
            #print(initial_k, new_method, difference, residue)
            if (abs(difference) < cutoff):
                final_pc_list_min.append(initial_pc)
            

            #updating the k
            initial_pc -= 0.1
            
    
  


    return(final_pc_list_max, final_pc_list_min) 
            
            
def method_range_pc():
    print('''
    This function, given a range ok k givan as an input, returns the pcs for which
    the (k+1)-th/k_th moment of the size plotted against abs(p-pc) (in log-log).
    We accept pcs that have the return of the custom method not being too far from a certain value

    Takes in input,
    1) final_k_dict: despite the name, an array ok ks to test
    2) densities: an array of densities
    3) criritical_pcs: The pcs from were we start the search for each k
    4) residues: an array of values, one for every k, that represent the starting method value
    for that k, used as a reference value for the search of the othe pcs
    5) dataframe: The dataframe from were to read the size, has the format specified in info_dataframe
    6) cutoffs: if the value returned by methods is distant more the cutoffs from the reference value
    it gets discarded
    7) methods=sum_method: The custom method acting on the residues to evaluate linearity
    ''')

#Now we want to find the pcs were the lines are more parallel    

def searching_phi(range_pc_max, range_k, densities, dataframe):

    slope_std_list = {}
    slope_phi_list = {}

    min_std = []
    #min_phi_std = []
    
    #for pc_pos in pc_positions:
    pc_to_try = range_pc_max.copy()
    k_to_span = range_k.copy()

    slope_std_list = {}
    slope_phi_list= {}

        
    for pc in pc_to_try:
        slope_std_list[pc] = []
        slope_phi_list[pc] = []

        k_phi = []
        
        for k in k_to_span:
            dens = np.abs(np.array(densities) - pc)
            rat_array = ratio_matrix([k], densities, dataframe)
            results = scp.stats.linregress(np.log(dens), np.log(rat_array[0,:])[:] )
            k_phi.append(results.slope)

        slope_std_list[pc].append(np.std(k_phi))
        slope_phi_list[pc].append(np.mean(k_phi))

    
        
    min_std = np.min(list(slope_std_list.values()))
    #pc that corresponds to the minimum std
    min_phi_std = min(slope_std_list, key=slope_std_list.get)
    #print(slope_std_list[min_phi_std])
    
    
    final_phi = slope_phi_list[min_phi_std]
    #
    #std_values_min = slope_std_list[min_phi_std]

    

    
    #min_std_dict = pd.DataFrame.from_dict(min_std, orient='Index', columns=['pc_s'])
    #min_std_dict['Index'] = min_std_dict.index

    #min_phi_std_dict = pd.DataFrame.from_dict(min_phi_std, orient='Index', columns=['phi_s'])
    #min_phi_std_dict['Index'] = min_std_dict.index

    #std_values_min_dict = pd.DataFrame.from_dict(std_values_min, orient='Index', columns=['std_values'])
    #std_values_min_dict['Index'] = std_values_min_dict.index

    #displacement_dict = pd.DataFrame.from_dict(displacement, orient='Index', columns=['Displacement'])
    #displacement_dict['Index'] = displacement_dict.index


    #final_df = min_std_dict.copy()

    #final_df = final_df.merge(min_phi_std_dict, how='inner', on='Index')

    #final_df = final_df.merge(std_values_min_dict, how='inner', on='Index')

    #final_df = final_df.merge(displacement_dict, how='inner', on='Index')


        
    return(final_phi, min_phi_std, min_std)




def info_searching_phi():
    print('''
    This function, for each k and also having the range of pcs, the pc for which the 
    (k+1)-th moment/k-th moment plotted as a function of the densities is straigther (in log-log),.
    It returns the slope found, the std of the optimal set of lines and the best pc.
    
    Takes in input:
    1) range_pc_max: the array of suitable pcs
    2) range_k: The range of suitable ks
    3) densities: an array of densities
    4) dataframe: the dataframe from which to retrieve the size (the dataframe has to have the format
    specified in info_dataframe)
    ''')


