
import numpy as np
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter



def func(k, dens, densities_plot, dataframe):
    
    def ratio(k_list, dens_list, df):
        #Generating the matrix to fill:
        #it will have a row for every value of k and a column
        #for each density

        ratios = np.empty((len(k_list), len(dens_list)))
        for i, k in enumerate(k_list):
            
            #Here we take from our dataframe the columns of interest
            #for each density
            for j, density in enumerate(dens_list):
                for column in df.columns:
                    if str(density) == column[0:column.find(' ')]:
                        if 'delta P' in column:
                            delta_p = df[column]
                        elif 'delta s' in column:
                            delta_s = df[column]
                        elif 'L-DNA' in column:
                            size = df[column]
                        else:
                            continue

                #Here first we treat the nans, then we compute
                #the k and the k+1 moment returning their fraction,
                
                def moment_k(k_value, dp, ds, size):
                    y_1 = size**k_value*(-1*(dp/ds))
                    y_2 = size**(k_value+1)*(-1*(dp/ds))
                    
                    #Here we mask the sizes with the indexes
                    #of the ys because we are deleting
                    #the corresponding records that are corrupted in y
                    #(and not in size but we can't use them)
                    size_1 = size.copy()
                    size_1 = size_1[~np.isinf(y_1)]
                    size_2 = size.copy()
                    size_2 = size_2[~np.isinf(y_2)]
                    size_1 = size_1[~np.isnan(y_1)]
                    size_2 = size_2[~np.isnan(y_2)]
                    y_1 = y_1[~np.isnan(y_1)]
                    y_2 = y_2[~np.isnan(y_2)]
                    y_1 = y_1[~np.isinf(y_1)]
                    y_2 = y_2[~np.isinf(y_2)]
                    
                    

                    moment_1 = np.trapz(y=y_1, x=size_1)
                    moment_2 = np.trapz(y=y_2, x=size_2)

                    
                    return moment_2/moment_1
                ratios[i,j] = moment_k(k, delta_p, delta_s, size)
        return ratios

    #Here we generate the matrix of ratios to give as an input to the
    #plotting function
    tmp = ratio(k, dens, dataframe)

    #Here below is the code for the plotting
    def plotting(ratios, densities, k_list):
        fig, ax = plt.subplots(1,1 , figsize=(15, 10))
        for i, column in enumerate(ratios):
            ax.scatter(densities, ratios[i,:] , marker='o'  , s=200,  edgecolor='black', alpha=0.7, label=f'k = {round(k_list[i],3)}' )

            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

            ax.xaxis.set_tick_params(which='minor', size=5)
            ax.xaxis.set_tick_params(which='major', labelsize=14 , size=10)
            ax.yaxis.set_tick_params(which='minor', size=5)
            ax.yaxis.set_tick_params(labelsize=14, size=10)

            ax.xaxis.get_offset_text().set_fontsize(16)
            ax.yaxis.get_offset_text().set_fontsize(16)

            ax.set_title("Moments ratio - DEAD-box helicase - LOG_LOG", fontsize=20)
            ax.set_ylabel('log(K+1 moment / k moment)', fontsize=16) 

            labels=[str(i) for i in densities]
            ax.set_xticks(densities)
            ax.set_xticklabels(labels=labels)
            ax.legend(shadow=True, framealpha=1, facecolor='aliceblue', edgecolor='black',  prop={'weight':'bold','size':16})
        return fig, ax
    return ratio(k, dens, dataframe).view(), plotting(tmp, densities_plot, k)




def info_func():
    print('''
    func expects in total 4 input: 

    k: the set of k that has to be a list or array-like

    dens: array-like of densities 

    densities_plot: the densities of the x_axis of the plot

    dataframe: the dataframe containing the data, the densities have to be at the start of the name of the columns
    and then a space has to follow, the delta P, delta s and size columns have to contain, respectively,
    delta P, delta s, L-DNA




    ''')






def ratio_matrix(k_list, dens_list, df):
#Generating the matrix to fill:
#it will have a row for every value of k and a column
#for each density

    ratios = np.empty((len(k_list), len(dens_list)))
    for i, k in enumerate(k_list):
        
        #Here we take from our dataframe the columns of interest
        #for each density
        for j, density in enumerate(dens_list):
            for column in df.columns:
                if str(density) == column[0:column.find(' ')]:
                    if 'delta P' in column:
                        delta_p = df[column]
                    elif 'delta s' in column:
                        delta_s = df[column]
                    elif 'L-DNA' in column:
                        size = df[column]
                    else:
                        continue

            #Here first we treat the nans, then we compute
            #the k and the k+1 moment returning their fraction,
            
            def moment_k(k_value, dp, ds, size):
                y_1 = size**k_value*(-1*(dp/ds))
                y_2 = size**(k_value+1)*(-1*(dp/ds))
                
                #Here we mask the sizes with the indexes
                #of the ys because we are deleting
                #the corresponding records that are corrupted in y
                #(and not in size but we can't use them)
                size_1 = size.copy()
                size_1 = size_1[~np.isinf(y_1)]
                size_2 = size.copy()
                size_2 = size_2[~np.isinf(y_2)]
                size_1 = size_1[~np.isnan(y_1)]
                size_2 = size_2[~np.isnan(y_2)]
                y_1 = y_1[~np.isnan(y_1)]
                y_2 = y_2[~np.isnan(y_2)]
                y_1 = y_1[~np.isinf(y_1)]
                y_2 = y_2[~np.isinf(y_2)]
                
                

                moment_1 = np.trapz(y=y_1, x=size_1)
                moment_2 = np.trapz(y=y_2, x=size_2)

                
                return moment_2/moment_1
            ratios[i,j] = moment_k(k, delta_p, delta_s, size)
    return ratios



def info_ratio_matrix():
    print('''
    This function returns the ratio matrix, the same used for potting in func
    
    Takes in input:
    1) k_list: an array of ks
    2) dens_list: an array of densities
    3) df: the dataframe with the same format described in info_df
    
    
    
    
    
    
    ''')




def moment(k_list, dens_list, df):
#Generating the matrix to fill:
#it will have a row for every value of k and a column
#for each density

    
    ratios = np.empty((len(k_list), len(dens_list)))
    for i, k in enumerate(k_list):
        
        #Here we take from our dataframe the columns of interest
        #for each density
        for j, density in enumerate(dens_list):
            for column in df.columns:
                if str(density) == column[0:column.find(' ')]:
                    if 'delta P' in column:
                        delta_p = df[column]
                    elif 'delta s' in column:
                        delta_s = df[column]
                    elif 'L-DNA' in column:
                        size = df[column]
                    else:
                        continue

            #Here first we treat the nans, then we compute
            #the k and the k+1 moment returning their fraction,
            
            def moment_k(k_value, dp, ds, size):
                y_1 = size**k_value*(-1*(dp/ds))
                y_2 = size**(k_value+1)*(-1*(dp/ds))
                
                #Here we mask the sizes with the indexes
                #of the ys because we are deleting
                #the corresponding records that are corrupted in y
                #(and not in size but we can't use them)
                size_1 = size.copy()
                size_1 = size_1[~np.isinf(y_1)]
                size_2 = size.copy()
                size_2 = size_2[~np.isinf(y_2)]
                size_1 = size_1[~np.isnan(y_1)]
                size_2 = size_2[~np.isnan(y_2)]
                y_1 = y_1[~np.isnan(y_1)]
                y_2 = y_2[~np.isnan(y_2)]
                y_1 = y_1[~np.isinf(y_1)]
                y_2 = y_2[~np.isinf(y_2)]
                
                

                moment_1 = np.trapz(y=y_1, x=size_1)
                

                
                return moment_1
            ratios[i,j] = moment_k(k, delta_p, delta_s, size)

    return ratios


def info_moment():
    print('''
    Returns the k-th momenths of an array of ks

    Takes in input:
    1) k_list: The array of ks to compite the moments on
    2) dens_list: The array of densities
    3) df: the dataframe with the same format described in info_df
    ''')