import cPickle as pic
from collections import defaultdict
#import pylab
import operator
import sys
import numpy as np
import string
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from numpy import ma
from  matplotlib import cbook

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

transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5, 'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15, 'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}
transform_inv = {0: 'STOP', 1: 'W', 2: 'F', 3: 'Y', 4: 'L', 5: 'I', 6: 'M', 7: 'V', 8: 'C', 9: 'A', 10: 'G', 11: 'P', 12: 'S', 13: 'T', 14: 'N', 15: 'Q', 16: 'H', 17: 'R', 18: 'K', 19: 'D', 20: 'E'}

WT_seq = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL*'


matrix = pic.load(open(sys.argv[1]))
matrix = np.ma.masked_invalid(matrix)

#pic.dump(matrix,open('Parul_data_matrix.pkl','wb'))

labels = np.arange(20)
def plot(start,stop):
    fig=plt.figure(1)
    ax = plt.subplot(111)
    ylabels = range(21)
    for index in ylabels:
        ylabels[index]= transform_inv[index]
    ax.set_xticks(np.arange(0,55, 5))
    ax.tick_params(axis = 'x', which= 'major', length=2, top=0)   
    ax.set_xticklabels(np.arange(start,stop, 5), size = 11)
    ax.set_yticks(np.array(range(len(ylabels))))
    ax.set_yticklabels(ylabels, fontsize = 5, fontweight='bold')
    cmap = plt.get_cmap('RdBu')
    cmap.set_under(cmap(0))
    cmap.set_over(cmap(1.0))
    cmap.set_bad(color = 'grey', alpha = 0.5)
    for i in np.arange(0.5,55.5,1):    
        ax.axvline(i, color='k', lw= 0.8)
    for i in np.arange(0.5,21.5,1):
        ax.axhline(i, color = 'k', lw = 0.8)
    norm = MidPointNorm(midpoint=0, vmin =np.nanmin(matrix), vmax = np.nanmax(matrix))
    #norm = MidPointNorm(midpoint=0, vmin = -3.5, vmax = 1)
    colors = cmap(norm(matrix))
    for index, letter in enumerate(WT_seq[start-1:stop]):
        pos = index
        letter_val = transform[letter]
        circ = plt.Circle((pos,letter_val), radius = 0.15, color = 'g')
        ax.add_patch(circ)
    # for index, letter in enumerate(WT_seq[start-1:stop]):
    #     pos = index
    #     letter_val = transform[letter]

    colors = colors[:,start-1:stop]
    img = plt.imshow(colors, interpolation = 'nearest', cmap = cmap, norm = norm)
    #ax = fig.add_axes([0.05, 0.50, 0.9, 0.2])
    #cb = matplotlib.colorbar.ColorbarBase(ax, extend = 'both',norm =norm, cmap = cmap, ticks = [-3.5,-2.5,-1.5,-0.5,0, 1], orientation = 'horizontal' )       
    #plt.savefig('%s_%s_Parul_WTnorm.png'%(start, stop), dpi=300, bbox='tight')
    plt.show()


step = 55
i = 1
while i < len(WT_seq):
	plot(i,i+step)
	i+=step
#plot(1,55)
# plot(56,110)
# plot(111,165)
# plot(166,220)
fig = plt.figure()
ax = fig.add_axes([0.05, 0.50, 0.9, 0.2])
cmap = plt.get_cmap('RdBu')
cmap.set_under(cmap(0))
cmap.set_over(cmap(1.0))
norm = MidPointNorm(midpoint=0, vmin =np.nanmin(matrix), vmax = np.nanmax(matrix))
cb = matplotlib.colorbar.ColorbarBase(ax, extend = 'both',norm =norm, cmap = cmap, orientation = 'horizontal' )       
plt.show()
#plt.savefig('GSP1_colorbar_no_floor.png', dpi = 300)
#print count_dict