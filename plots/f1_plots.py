
#code to generate Figure S18 A-C
import os
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
#0 is the larger clsuter
# 2 is the unusued clsuter

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


parent = os.path.join(os.path.abspath(__file__), os.pardir)

filename = os.path.abspath(os.path.join(parent,'..', 
                                        'Final_Data', 
                                        'test_DescSelect_0.xlsx'))

data = pd.read_excel(filename, sheet_name = 'F1')



headers = list(data.columns)
headers.remove('Unnamed: 0')
data = data[headers]
aveF1 = data.mean(axis = 1)
stdF1 = data.std(axis = 1)


fig, (ax3,ax4,ax5)  = plt.subplots(1,3, figsize = (12,4))

x3 = np.linspace(10, (len(aveF1)+9), len(aveF1))

#F1 plot with std shading
tprs_upper = np.minimum(aveF1 + stdF1, 1)
tprs_lower = np.maximum(aveF1 - stdF1, 0)
ax3.plot(x3,aveF1, color='red',lw=2, alpha=.8, label = r'Average F1 Score')

ax3.fill_between(x3, tprs_lower, tprs_upper,color='blue', alpha = 0.2, label=r'$\pm$ 1 std. dev.')
ax3.set_xlabel('Number of Descriptors', fontsize = 12)
ax3.set_ylabel('F1 Score', fontsize = 12)
ax3.legend(loc="lower right", fontsize = 10)

ax3.text(0.2,.95,"(A)", fontsize = 12, color = 'black', horizontalalignment='right',
        verticalalignment='top', transform=ax3.transAxes)

ax3.set_ylim([0.6,1])

#F1 smoothed with moving average
x4 = np.linspace(14, (len(aveF1)+4), len(aveF1)-6)
moveAve = moving_average(aveF1, 7)

ax4.plot(x4,moveAve, color='red',lw=2, alpha=.8)
ax4.set_ylabel('F1 Score Moving Average', fontsize = 12)
ax4.set_xlabel('Number of Descriptors', fontsize = 12)
ax4.text(0.2,.95,"(B)", fontsize = 14, color = 'black', horizontalalignment='right',
        verticalalignment='top', transform=ax4.transAxes)
ax4.set_ylim([0.6,1])


#F1 first derivative
derive_move_ave = np.diff(moveAve)
x5 = np.linspace(14, (len(aveF1)+4), len(aveF1)-7)

ax5.plot(x5,derive_move_ave, color='red',lw=2, alpha=.8)
ax5.set_xlabel('Number of Descriptors', fontsize = 12)
ax5.set_ylabel('1st Derivative of Moving Average', fontsize = 12)
ax5.text(0.2,.95,"(C)", fontsize = 14, color = 'black', horizontalalignment='right',
        verticalalignment='top', transform=ax5.transAxes)

ax5.set_ylim([-0.0045,0.02])
ytickVals = np.arange(-0.0025, 0.02, 0.0025)
ax5.set_yticks(ytickVals)


fig.subplots_adjust(wspace =0.5, bottom = 0.2)
plt.show()

fig.savefig(r'F1_plots.png', format='png', dpi=1200)

