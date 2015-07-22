from __future__ import division
import scipy as sp
import matplotlib.pyplot as plt
import os
from scipy import stats
import histogram_functions as hif

# Wrapper functions
def fn(x,a,b):
    return stats.norm.pdf(x,loc=a,scale=b)
    
def fln(x,a,b):
    mu = a
    sigma = b
    return 1/(x*sigma*sp.sqrt(2*sp.pi))*sp.exp(-(sp.log(x)-mu)**2/(2*sigma**2))
#    return stats.lognorm.pdf(x, a, loc=0, scale=b)

#    pdf(x, s, loc=0, scale=1)    
#    lognorm.pdf(x, s) = 1 / (s*x*sqrt(2*pi)) * exp(-1/2*(log(x)/s)**2)
#    for x > 0, s > 0.
#    
#    If log(x) is normally distributed with mean mu and variance sigma**2, 
#    then x is log-normally distributed with shape parameter sigma and scale parameter exp(mu).
    
    

def plot_histogram(bins,bars,color,sub,xlabel,ylabel,text,xtics,yticks):

    alpha=0.1
    bin_width = bins[1]-bins[0]    
    fig = plt.gcf()
    ax = fig.add_subplot(sub)
    ax.bar(bins,bars,align='center',
           width=bin_width,alpha=alpha,color=color)

    if text == 'EXP001':
        ax.plot([1,1],[0,30],'k:')
    else:
        ax.plot([1,1],[0,50],'k:')
    
    if text in ['EXP001','EXP008e3','SUGRA003']:
        plt.title(text)
    else:
        ax.text(0.88,31,text,size=9)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axis([0.85,1.15,0,50])
    plt.xticks([0.9,1,1.1])
    plt.yticks(sp.arange(0,51,10))
    if not xticks:
        plt.setp( ax.get_xticklabels(), visible=False)

    if not yticks:
        plt.setp( ax.get_yticklabels(), visible=False)
    
    return ax
    
os.chdir('/home/io/Dropbox/Projekter/Hubble/VelocityField/code/')

#plotlabel = "otherdists"
#models = ["EXP001", "EXP008e3", "SUGRA003"]
#min_dists = ["5.0","5.0","5.0"]
#texts = ["EXP001", "EXP008e3", "SUGRA003"]
#subs =   [131, 132, 133]
#xtickss = [1,1,1]
#ytickss = [1,0,0]
#xlabel = "$H_{loc}/H_0$"
#xlabels = [xlabel, xlabel, xlabel]
#ylabel = "$P(H_{loc}/H_0)$"
#ylabels = [ylabel, "", ""]
#plt.figure(figsize=(6,2))

plotlabel = "LCDM"
models = ["LCDM", "LCDM", "LCDM","LCDM"]

#plotlabel = "EXP003"
#models = ["EXP003", "EXP003", "EXP003","EXP003"]
#
min_dists = ["0.0", "1.0", "5.0", "10.0"]
texts = ["R="+min_dists[0]+"Mpc/h","R="+min_dists[1]+"Mpc/h",
         "R="+min_dists[2]+"Mpc/h","R="+min_dists[3]+"Mpc/h"]
subs =   [221, 222, 223, 224]
xtickss = [0,0,1,1]
ytickss = [1,0,1,0]
xlabel = "$H_{loc}/H_0$"
xlabels = ["","",xlabel,xlabel]
ylabel = "$P(H_{loc}/H_0)$"
ylabels = [ylabel,"",ylabel,""]
plt.figure(figsize=(6,5))


func = fln
nbins = 50

plt.subplots_adjust(hspace=.001,wspace=.001)

for i in range(len(models)):

    print "-----------------"
    model = models[i]
    print "model =" , model
    min_dist = min_dists[i]
    text = texts[i]
    sub = subs[i]
    xlabel = xlabels[i]
    ylabel = ylabels[i]
    xticks = xtickss[i]
    yticks = ytickss[i]
    
    case = "CoDECS_"+model
    path = '../cases/'+case+'/'
    fname = path+min_dist+'Hubbleconstants.txt'
    H_uncorrected = sp.loadtxt(fname)[:,1]/100
    H_corrected = sp.loadtxt(fname)[:,2]/100
    
    bins, bars = hif.make_normalised_histogram(nbins,H_uncorrected)
    ax = plot_histogram(bins,bars,'b',sub,xlabel,ylabel,text,xticks,yticks)
    bins_fine,fitted_bars = hif.fit_distribution_to_normalised_histogram(bins,bars,func)
    ax.plot(bins_fine,fitted_bars,'b',label='Uncorrected')
    
    things_be_crazy = True
    if things_be_crazy:
        indices = (sp.absolute(H_corrected) < 100) & (H_corrected > 0) 
        H_corrected = H_corrected[indices]
        H_uncorrected = H_uncorrected[indices]
        print "corr:", max(H_corrected), min(H_corrected)
        print "uncorr:", max(H_uncorrected), min(H_uncorrected)
        print "goodness ratio goods/bads = %s/%s = %s" % (len(H_corrected),len(indices),len(H_corrected)/len(indices))
        
    print min_dist
        
    if not min_dist == '0.0':
        bins, bars = hif.make_normalised_histogram(nbins,H_corrected)
        ax = plot_histogram(bins,bars,'g',sub,xlabel,ylabel,text,xticks,yticks)
        bins_fine, fitted_bars = hif.fit_distribution_to_normalised_histogram(bins,bars,func)
        ax.plot(bins_fine,fitted_bars,'g',label='Corrected')

    if not model == 'LCDM':
        case_LCDM = "CoDECS_LCDM"
        path_LCDM = '../cases/'+case_LCDM+'/'
        fname_LCDM = path_LCDM+min_dist+'Hubbleconstants.txt'
        H_uncorrected_LCDM = sp.loadtxt(fname_LCDM)[:,1]/100

        
        bins_LCDM, bars_LCDM = hif.make_normalised_histogram(nbins,H_uncorrected_LCDM)
        bins_fine_LCDM, fitted_bars_LCDM = hif.fit_distribution_to_normalised_histogram(bins_LCDM, bars_LCDM, func)
        ax.plot(bins_fine_LCDM,fitted_bars_LCDM,'m',linewidth=1,alpha=0.5,label='$\Lambda$CDM')

        if not min_dist == '0.0':
            H_corrected_LCDM = sp.loadtxt(fname_LCDM)[:,2]/100        
            bins_LCDM, bars_LCDM = hif.make_normalised_histogram(nbins,H_corrected_LCDM)
            bins_fine_LCDM, fitted_bars_LCDM = hif.fit_distribution_to_normalised_histogram(bins_LCDM,bars_LCDM,func)
            ax.plot(bins_fine_LCDM, fitted_bars_LCDM, 'm',linewidth=1,alpha=0.5)

    if text in ['EXP001',"R="+min_dists[1]+"Mpc/h"]:    
        plt.legend(bbox_to_anchor=(1, 1),prop={'size':8},frameon=False)

    
plt.savefig("/home/io/Dropbox/SharedStuff/Cepheids/Hdists_"+plotlabel+".pdf",bbox_inches='tight')
