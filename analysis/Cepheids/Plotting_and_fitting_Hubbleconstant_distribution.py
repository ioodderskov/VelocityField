from __future__ import division
from scipy import stats # Import the scipy.stats module
from scipy.optimize import curve_fit # Import the curve fitting module
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy as sp
import scipy.integrate as integrate
#import pdb

os.chdir('/home/io/Dropbox/Projekter/Hubble/VelocityField/code/')
#model = 'EXP008e3'
#case = 'CoDECS_'+model
case = 'CoDECS_LCDM'
path = '../cases/'+case+'/'
data = sp.loadtxt(path+'Hubbleconstants.txt')
H_notcorrected = data[:,1]/100
H_corrected = data[:,2]/100



#cut_tails = sp.absolute(H_fullsky-1) < 0.1 
#H_cut_tails = H_fullsky[cut_tails]
#H_fullsky = H_cut_tails
#plot_H_distribution(H_fullsky)



## SET UP INPUT DATA
# Define a lognormal probability distribution from median and geometric
# standard deviation (see link below for tips on using lognormal 
# distributions in Python)
#s = 2.0 # Geometric standard deviation
#M = 40.0 # Geometric mean == median
#shape_in = np.log(s) # Scipy's shape parameter corresponds to log(geo. std. dev.)
#scale_in = M # Scipy's scale parameter corresponds to median
# Calculate histogram for random variates (this is equivalent to my input data)
#xx = stats.lognorm.rvs(shape_in, loc=0, scale=scale_in, size=1000) # Generate data
def make_histogram(ax1,H,color,alpha):
    dist = getattr(stats,'lognorm')
    number_of_bins = 20
    xx = H
    hist, bin_edges = np.histogram(xx, bins=number_of_bins, density=True) # Calculate histogram
    x_hist = bin_edges[1:] # Take the upper edges of the bins
    y_hist = hist.cumsum()/hist.cumsum().max()  # Normalise the cumulative sum
    y_hist = y_hist# add some noise
    
    
    ## FIT THE DISTRIBUTION
    (shape_out, scale_out), pcov = curve_fit(
                lambda xdata, shape, scale: dist.cdf(xdata, shape, loc=0, scale=scale),
                x_hist, y_hist)
    
    ## HOW IT WORKS
    # curve_fit takes a *function*, xdata and ydata as inputs.  It then varies the inputs to
    # the *function* to get the best fit.  In this case, we want to fit the cumulative density
    # function, keeping the location parameter fixed at 0.  We use the lambda function as our
    # input function.  It is a temporary function that is just the cumulative density function,
    # but with the location input already set.
    
    ## PLOT THE GRAPH

    #ax1.hist(xx, bins=25, normed=True,color='green', alpha=0.5, label='Binned data')
    
    
    histogram = sp.histogram(H,bins=number_of_bins,normed=True)
    numbers_in_bins = histogram[0]
#    fractions_in_bins = numbers_in_bins/sp.sum(numbers_in_bins)
    bin_locations = histogram[1]
    bin_width = bin_locations[1]-bin_locations[0]
    bin_centers = (bin_locations[:-1]+bin_locations[1:])/2
    factor = sp.sum(numbers_in_bins)
    plt.bar(bin_centers, numbers_in_bins/factor, align = 'center', 
            width=bin_width, alpha=alpha, color=color)
    #plt.title(model)
    #    plt.hist(H,bins=bin_centers,normed=True)
    
    #pmax = max(numbers_in_bins/factor)+0.05
#    pmax = 0.1
    
    #x = np.linspace(0.7, 1.3, num=300) # values for x-axis
    x = np.linspace(0.9, 1.1, num=300) # values for x-axis
    fitted_distribution = dist.pdf(x, shape_out, loc=0, scale=scale_out)
    #plt.xlim(0,150)
    print "scale = ", scale_out
    print "shape = ", shape_out
    
    def fitted_lognorm(x,mu,sigma):
        return 1./(x*sigma*sp.sqrt(2*sp.pi))*sp.exp(-(sp.log(x)-mu)**2/(2*sigma**2))
        
    mu = sp.log(scale_out)
    sigma = shape_out

    label = '$\exp(\mu)=%0.4f, \sigma=%0.4f$' % (sp.exp(mu),sigma)    
    ax1.plot(x,fitted_distribution/factor, color=color, lw=2, 
             label=label)
    plt.legend()

    
    print "mu, sigma = ", mu,sigma
    print "exp(mu) = ", sp.exp(mu)
    
    print "The probability of measuring a Hubble constant 5% or more too low (high) is",
    print integrate.quad(lambda x: fitted_lognorm(x,mu,sigma), 0, 0.95)
    print "(",integrate.quad(lambda x: fitted_lognorm(x,mu,sigma), 1.05,20),")"
    
    print "Checking that the fitted distribution is normalised"
#    print integrate.quad(lambda x: fitted_lognorm(x,mu,sigma),0,2)
    print integrate.quad(lambda x: dist.pdf(x, shape_out, loc=0, scale=scale_out),0,2) 
    
    print "Checking that the histogram is normalised:", sp.sum(numbers_in_bins/factor)
    print "factor = ", factor
    #leg=ax1.legend()
    #ax2=plt.subplot(122)
    #ax2.plot(x_hist, y_hist, 'ro', label='Binned data points')
    #ax2.plot(x,dist.cdf(x, shape_out, loc=0, scale=scale_out), label='Fitted distribution')
    #plt.xlim(0,150)
    #ax2.set_ylim(0,1.0)
    #plt.xlabel('x')
    #plt.ylabel('y_cdf')
    #leg=ax2.legend(loc='lower right', numpoints=1)
    #results_txt="""M_in=%.2f
    #M_out=%.2f
    #s_in=%.2f
    #s_out=%.2f""" % (M, scale_out, s, np.exp(shape_out))
    #txt=plt.text(0.97, 0.3, results_txt, transform=ax2.transAxes, horizontalalignment='right', fontsize='large')
    #plt.savefig(path+'H_distribution_coma_'+case+'.pdf')
    #plt.savefig('/home/io/Desktop/H_distribution_nocorrections.pdf')
    
    
    

plt.figure(figsize=(6,5))
ax1=plt.subplot(111)

plt.xlabel('$H_{loc}/H_0$',fontsize=16)
plt.ylabel('P($H_{loc}/H_0$)',fontsize=16)
plt.title(case)
#plt.axis([0.85,1.15,0,pmax])
plt.axis([0.9,1.1,0,0.4])
plt.plot([1,1],[0,0.4],'k--',linewidth=1.5)

    
make_histogram(ax1,H_notcorrected,'g',0.4)
make_histogram(ax1,H_corrected,'b',0.4)

plt.savefig("/home/io/Desktop/distribution_of_hubbleconstants"+case+".pdf")