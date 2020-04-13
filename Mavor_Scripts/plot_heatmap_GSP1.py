###############################
# Usage - "python csv_to_matrix_HSP90.py file_name.csv"
###############################

###############################
# Imports 
###############################

import sys
import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from numpy import ma
from  matplotlib import cbook
import cPickle as pic

###############################
# Parameters - set plotting parameters here
###############################
# WT sequence using 1 letter code and "*" for STOP
WT_seq = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL*'
# Plotting Params - set start and window lenght for plotting
plot_start = 1
plot_stop = 220
#plot_stop = len(WT_seq) # uncomment to plot to the end of the WT protein
plot_window_lenght = 55
# Color Params - set floor and ceiling for coloring the slopes
# floor = -0.5
# ceiling = 0.2
# Uncomment the paramaters bellow to use the full range of slopes
floor = None
ceiling = None
# Output file string - set a string to be appended to all the output file names
output_string = sys.argv[1][:-4]


###############################
# Classes & Functions
###############################
class MidPointNorm(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result

def plot(start,stop,floor, ceiling): # function for plotting heatmaps
    fig=plt.figure()
    plt.title(output_string)
    ax = plt.subplot(111)
    ylabels = range(21)
    for index in ylabels:
        ylabels[index]= transform_inv[index]
    ax.set_xticks(np.arange(0,plot_window_lenght, 5))
    ax.tick_params(axis = 'x', which= 'major', length=2, top=0)   
    ax.set_xticklabels(np.arange(start,stop, 5), size = 11)
    ax.set_yticks(np.array(range(len(ylabels))))
    ax.set_yticklabels(ylabels, fontsize = 5, fontweight='bold')
    cmap = plt.get_cmap('RdBu')
    cmap.set_under(cmap(0))
    cmap.set_over(cmap(1.0))
    cmap.set_bad(color = 'grey', alpha = 0.5)
    for i in np.arange(0.5,plot_window_lenght + 0.5,1):    
        ax.axvline(i, color='k', lw= 0.8)
    for i in np.arange(0.5,21.5,1):
        ax.axhline(i, color = 'k', lw = 0.8)
    norm = MidPointNorm(midpoint=0, vmin = floor, vmax = ceiling)
    colors = cmap(norm(matrix))
    for index, letter in enumerate(WT_seq):
        pos = index
        letter_val = transform[letter]
        colors[letter_val][pos] = (0.15,0.9,0.15,0.65)
    for index, letter in enumerate(WT_seq[start-1:stop]):
        pos = index
        letter_val = transform[letter]
        circ = plt.Circle((pos,letter_val), radius = 0.15, color = 'k')
        ax.add_patch(circ)
    colors = colors[:,start-1:stop]
    img = plt.imshow(colors, interpolation = 'nearest', cmap = cmap, norm = norm)
    #plt.colorbar(extend = 'both', orientation ='horizontal', label = 'Fitness')
    print output_string
    plt.savefig('%s_%s_%s.png'%(output_string, start, stop), dpi=300, bbox='tight') # comment me if you would not like to save the plots
    #plt.show() # uncomment me if you would like to see the plots
###############################
# Dictionaries
###############################
transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}
transform_inv = {0: 'STOP', 1: 'W', 2: 'F', 3: 'Y', 4: 'L', 5: 'I', 6: 'M', 7: 'V', 8: 'C', 9: 'A', 10: 'G', 11: 'P', 12: 'S', 13: 'T', 14: 'N', 15: 'Q', 16: 'H', 17: 'R', 18: 'K', 19: 'D', 20: 'E'}

###############################
# Main
###############################
# f = open(sys.argv[1], 'rb') # input data
# f.readline() # skip header
# f.readline() # skip header
# slope_dict = {} # creat dictionary for slopes

# for line in f:
#     l_list = line.strip().split(',')
#     try:
#         pos = int(l_list[0]) # position
#         mut = l_list[1] # mutation
#         slope = float(l_list[2]) # slope
#         if slope == 999.9: # these are skipped and will become nan in the matrix
#             pass
#         else:
#             slope_dict[(pos,mut)] = slope # populate dictionary
#     except ValueError: # valueError is thrown when the rows containing mutation --> slope association ends. Input stops here
#         pass

# matrix = np.empty((21,len(WT_seq))) # creattes an empty matrix for heatmap 21 X sequence lenght
# matrix.fill(np.nan) # fills the matrix with nan

# for mutant in slope_dict: # replaces nan with a slope for each observation
#     pos = mutant[0]
#     mut = mutant[1]
#     matrix[transform[mut]][pos-1] = slope_dict[mutant]
matrix = pic.load(open(sys.argv[1],'rb'))
matrix = np.ma.masked_invalid(matrix) # mask nan values

i = plot_start # Do the plotting
while i <= plot_stop:
	plot(i,i+plot_window_lenght -1, floor, ceiling)
	i+=plot_window_lenght