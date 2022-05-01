import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib

def cumulative(array, title='Cumulative', xlabel='x axis', ylabel='y axis'):
    fig, ax = plt.subplots(1,1,figsize=(15,10))
    array = np.array(array)
    array = np.sort(array)
    array = array[~np.isnan(array)]
    cumul = 1 - np.arange(0, len(array))/(len(array))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.minorticks_on()
    ax.scatter(array, cumul, title, xlabel, ylabel, s=150, edgecolor='black', alpha=0.7)
    return fig.tight_layout(), ax



def info_cumulative():
    print('''
    The cumulativr function plots the survivability in a log-log scale.
    It accepts as an input the array to plot the cumulative over.

    
    It also accepts the following key-word arguments:

    scatter: a dictionary with the usual keyword arguments for a scatterplot


    Other keyword arguments to be added.
        
    
    
    
    
    
    ''')




def cumulative_data(array):
    array = np.array(array)
    array = np.sort(array)
    array = array[~np.isnan(array)]
    cumul = 1 - np.arange(0, len(array))/(len(array))
    return cumul


def info_cumulative_data():
    print('''
    The function cumulative data returnd the survivability ( p(X>x) ) of a given array provided in input
    
    
    
    
    
    
    
    ''')





def cumulative_matrix(array, densities, title='Cumulative', xlabel='x axis', ylabel='y axis', log_log = True, fontsize_title = 20, fontsize_x = 15, fontsize_y = 15, offset_x = 14, offset_y = 14, s=250, legend_size=14, unit_of_measure='$\\mu M$', k_round=2):
    fig, ax = plt.subplots(1,1,figsize=(15,10))
    for column, density in zip(range(array.shape[1]), densities):
        arr = array[:,column]
        arr = np.array(arr)
        arr = np.sort(arr)
        arr = arr[~np.isnan(arr)]
        cumul = 1 - np.arange(0, len(arr))/(len(arr))
        if log_log == True:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.yaxis.offsetText.set_fontsize(offset_y)
            ax.xaxis.offsetText.set_fontsize(offset_x)
        ax.minorticks_on()
        leg_num = round(float(density), k_round)
        ax.scatter(arr, cumul, s=s, edgecolor='black', alpha=0.5, label=f'{leg_num} {unit_of_measure}')
        ax.set_title(title, fontsize=fontsize_title)
        ax.set_xlabel(xlabel, fontsize=fontsize_x)
        ax.set_ylabel(ylabel, fontsize=fontsize_y)
        ax.xaxis.set_tick_params(labelsize=15, size=10)
        ax.xaxis.set_tick_params(which='minor', size=4)
        ax.yaxis.set_tick_params(labelsize=15, size=10)
        ax.yaxis.set_tick_params(which='minor', size=4)
        ax.legend(shadow = True, framealpha=1, facecolor='aliceblue', edgecolor='black',  prop={'weight':'bold','size':legend_size}, loc='best')
    return fig.tight_layout(), ax



#To include the fact that densities has to have the same number of columns
#as array


def cumulative_collapse(array, alpha, phi, densities, title='Cumulative', xlabel='x axis', ylabel='y axis', log_log = False, fontsize_title = 20, fontsize_x = 15, fontsize_y = 15, offset_x = 14, offset_y = 14, s=250, legend_size=14, unit_of_measure='$\\mu M$', k_round=2):
    fig, ax = plt.subplots(1,1,figsize=(15,10))
    for density, column in zip(densities, range(array.shape[1])):
        arr = array[:,column]
        arr = np.array(arr)
        arr = np.sort(arr)
        arr = arr[~np.isnan(arr)]
        cumul = 1 - np.arange(0, len(arr))/(len(arr))
        if log_log == True:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.yaxis.offsetText.set_fontsize(offset_y)
            ax.xaxis.offsetText.set_fontsize(offset_x)
        ax.minorticks_on()
        leg_num = round(float(density), k_round)
        ax.scatter(arr/(densities[column]**phi), cumul*(arr**alpha), s=s, edgecolor='black', alpha=0.7, label=f'{leg_num} {unit_of_measure}')
        ax.set_title(title, fontsize=fontsize_title)
        ax.set_xlabel(xlabel, fontsize=fontsize_x)
        ax.set_ylabel(ylabel, fontsize = fontsize_y)
        ax.xaxis.set_tick_params(labelsize=15, size=10)
        ax.xaxis.set_tick_params(which='minor', size=4)
        ax.yaxis.set_tick_params(labelsize=15, size=10)
        ax.yaxis.set_tick_params(which='minor', size=4)
        ax.legend(shadow=True, framealpha=1, facecolor='aliceblue', edgecolor='black',  prop={'weight':'bold','size':legend_size}, loc='best')
    return fig.tight_layout(), ax


