#Boxplot code
from turtle import title
import numpy as np
import matplotlib.pyplot as plt


def boxplot(y, colors, figsize=(15,10), linewidth=2, color_median='black', linewidth_median = 3, size_props = 12, title = 'boxplot', font_title = 20, xlabel='x_axis', x_font=15, ylabel='y_axis', y_font=15, labels=False, list_labels=[] ):
    
    fig, ax = plt.subplots(1,1,figsize=figsize)
    
    bp = ax.boxplot(y,boxprops=dict(linewidth=linewidth), patch_artist=True, medianprops=dict(color=color_median, linewidth=linewidth_median), flierprops=dict(markersize=size_props))
    
    
    ax.set_title(title, fontsize=font_title)
    

    ax.set_xlabel(xlabel, fontsize=x_font)

    ax.set_ylabel(ylabel, fontsize=y_font)

   
    

    #Still to insert the case were y is a matrix
    
    for array, x_point, color in zip(y, range(1, len(y)+1), colors) :
        tmp = np.ones([len(array),])*x_point + 0.35
        ax.scatter(tmp, array, s=77, alpha=0.5, edgecolor='black', color=color)

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
   

    ax.minorticks_on()
    ax.xaxis.set_tick_params(which='minor', bottom=False , size=10)
    ax.yaxis.set_tick_params(labelsize=14, size=10)
    ax.set_xticklabels(labels)

    if labels:
        ax.set_xticklabels(list_labels)

    
    return fig, ax




def info_boxplot():
    print('''
    The boxplot function takes only one mandatory argument: y that is either a matrix-like with the series to plot as columns or an array-like of array-likes.
    
    Then it accepts the following key-word arguments:
    1) y: a list of the sequences to plot
    2) colors: a list with the colors of the boxes,
    

    3) title, xlabel and ylabel: self explanatory, we have also as far as the fontsize font_title, 
        x_font, y_font.

    4) color_median: the color of the median
    5) linewidth_median: self explanatory
    6) size_props: size of the outliers
    
    ''')
