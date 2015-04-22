from __future__ import division
import sys
import scipy as sp
from scipy import stats # Import the scipy.stats module
from scipy.optimize import curve_fit # Import the curve fitting module
import matplotlib.pyplot as plt
import os
import scipy.integrate as integrate


if len(sys.argv) != 2:
    print "Wrong number of arguments"
	
case = sys.argv[1]




#data = []
#for observer in observers:
#    for velocity in observer.observed_radial_peculiar_velocities:
#        data.append(velocity)
path = '../cases/'+case+'/'
pairwise_velocities_file = path+'pairwise_velocities.npy' 
data = sp.load(pairwise_velocities_file)

#data = sp.array(data)
#cut_tails = sp.absolute(data) < 1000
#data_cut = data[cut_tails]

number_of_bins = 200
plt.figure(figsize=(6,5))
ax1=plt.subplot(111)

dist = getattr(stats,'norm')

xx = data
hist, bin_edges = sp.histogram(xx, bins=number_of_bins, density=True) # Calculate histogram
x_hist = bin_edges[1:] # Take the upper edges of the bins
y_hist = hist.cumsum()/hist.cumsum().max()  # Normalise the cumulative sum    


histogram = sp.histogram(data,bins=number_of_bins,normed=True)
numbers_in_bins = histogram[0]
fractions_in_bins = numbers_in_bins/sp.sum(numbers_in_bins)
bin_locations = histogram[1]
bin_width = bin_locations[1]-bin_locations[0]
bin_centers = (bin_locations[:-1]+bin_locations[1:])/2
factor = sp.sum(numbers_in_bins)
plt.bar(bin_centers, numbers_in_bins/factor, align = 'center', width=bin_width, alpha=0.5, color='g')
ymax = 0.08
plt.axis([-2500,2500,0,ymax])

#def cdf(xdata,loc,scale):
#    return stats.norm.cdf(xdata,loc=loc,scale=scale)
#def pdf(xdata,loc,scale):
#    return stats.norm.pdf(xdata,loc=loc,scale=scale)

#(mu_cdf, sigma_cdf), pcov = curve_fit(cdf, x_hist, y_hist,p0=[0,500])
#(mu_pdf, sigma_pdf), pcov = curve_fit(pdf, bin_centers, numbers_in_bins,p0=[0,500])    
mu, sigma = dist.fit(data)
#mu_cut, sigma_cut = dist.fit(data_cut)

x = sp.linspace(-2500, 2500, num=500) # values for x-axis
fitted_distribution = dist.pdf(x, loc=mu,scale=sigma)
#fitted_distribution_pdf = dist.pdf(x, loc=mu_pdf,scale=sigma_pdf)
#fitted_distribution_cdf = dist.pdf(x, loc=mu_cdf,scale=sigma_cdf)
#fitted_distribution_cut = dist.pdf(x, loc=mu_cut,scale=sigma_cut)
#
ax1.plot(x,fitted_distribution/factor, 'b', lw=2,
label='Fit: mu=%0.1f, sigma=%0.1f' %(mu,sigma))
#ax1.plot(x,fitted_distribution_pdf/factor, 'r', lw=2,
#label='Fit to pdf: mu=%0.1f, sigma=%0.1f' %(mu_pdf,sigma_pdf))
#ax1.plot(x,fitted_distribution_cdf/factor, 'black', lw=2,
#label='Fit to cdf: mu=%0.1f, sigma=%0.1f' %(mu_cdf,sigma_cdf))
#ax1.plot(x,fitted_distribution_cut/factor, 'y', lw=2,
#label='Cut data: mu=%0.1f, sigma=%0.1f' %(mu_cut,sigma_cut))
plt.legend()



plt.xlabel('$v_{||}$',fontsize=16)
plt.ylabel('$P(v_{||})$',fontsize=16)
plt.title("%s" % case)

plt.plot([1,1],[0,ymax],'k--',linewidth=1.5)

def fitted_norm(x,mu,sigma):
    return 1./(sp.sqrt(2*sp.pi)*sigma)*sp.exp(-(x-mu)**2/(2*sigma**2))


print "mu, sigma = ", mu,sigma

print "Checking that the histogram is normed"
print "sum(fractions_in_bins) = ", sp.sum(fractions_in_bins)

print "Checking that the fitted distributions are normalised"
mus = [mu]
sigmas = [sigma]
for m,s in zip(mus,sigmas):
    print integrate.quad(lambda x: fitted_norm(x,m,s),-3000,3000)

plt.savefig(path+'pairwise_velocity_distribution.pdf')

