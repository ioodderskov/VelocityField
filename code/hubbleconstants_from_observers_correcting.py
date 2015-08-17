from __future__ import division
import scipy as sp
from scipy.optimize import leastsq
import sys
from scipy import linalg
import pdb

if len(sys.argv) != 3:
    print "Wrong number of arguments"

case = sys.argv[1]
min_dist = sys.argv[2]
path = '../cases/'+case+'/'
#min_dist  = '5.0'


observers_file = path+min_dist+'observers_test_before_HS8.npy'
observers = sp.load(observers_file)

N = len(observers)
vrs_plus_hubbleflow = sp.zeros(N)
vrs_correction_plus_hubbleflow = sp.zeros(N)
vrs_correction = sp.zeros(N)
rs = sp.zeros(N)

observed_velocities = []
velocity_corrections = []


print "len(observers) = ", len(observers)

for index,observer in enumerate(observers):
    if len(observer.chosen_halos) == 0:
        continue
    position_op = observer.chosen_halos[0].position_op
    r = observer.chosen_halos[0].r
    v_peculiar_observed = observer.chosen_halos[0].observed_velocity
    velocity_correction = observer.chosen_halos[0].velocity_correction

    observed_velocities.append(v_peculiar_observed)
    velocity_corrections.append(velocity_correction)
    
    vr_correction = sp.dot(position_op-observer.position,velocity_correction)/r
    vr_peculiar_observed = sp.dot(position_op-observer.position,v_peculiar_observed)/r 
    
    vr_plus_hubbleflow = vr_peculiar_observed+r*100
    vr_correction_plus_hubbleflow = vr_correction + r*100
    
    vrs_plus_hubbleflow[index] = vr_plus_hubbleflow
    vrs_correction_plus_hubbleflow[index] = vr_correction_plus_hubbleflow
    vrs_correction[index] = vr_correction
    rs[index] = r
    
pdb.set_trace()    

    

	
#def error(factor,vr_plus_hubbleflow,vr_correction_plus_hubbleflow):
#    err = vr_plus_hubbleflow-factor*vr_correction_plus_hubbleflow
#    return err
#
#(factor,guess) = leastsq(error,2,args=(vrs_plus_hubbleflow,vrs_correction_plus_hubbleflow))

observed_velocity_norms = sp.array([linalg.norm(observed_velocity) for observed_velocity in observed_velocities])
velocity_correction_norms = sp.array([linalg.norm(velocity_correction) for velocity_correction in velocity_corrections])
#pdb.set_trace()
print "The ratios ranges from", sp.amin(observed_velocity_norms/velocity_correction_norms), "to", sp.amax(observed_velocity_norms/velocity_correction_norms)


#factor = sp.mean(observed_velocity_norms/velocity_correction_norms)
#print "The mean of the ratios is", factor
print "abs(observed_velocities)/abs(velocity_corrections) = ",sp.mean(sp.absolute(observed_velocities),axis=0)/sp.mean(sp.absolute(velocity_corrections),axis=0)

# Correction factor:
omega_m_false = 0.272
omega_m_true = 0.2711

h_false = 0.702
h_true = 0.703

omega_cdm = 0.226
gamma = 0.56
epsilon_c = 2.4

beta_c_dict = {}
beta_c_dict['CoDECS_LCDM'] = 0
beta_c_dict['CoDECS_EXP001'] = 0.05
beta_c_dict['CoDECS_EXP003'] = 0.15
beta_c_dict['CoDECS_EXP008e3'] = 0.4
beta_c_dict['CoDECS_SUGRA003'] = -0.15

beta_c = beta_c_dict[case]
print "beta_c = ", beta_c

growth_false = omega_m_false**0.55
growth_true = omega_m_true**gamma*(1+gamma*omega_cdm/omega_m_true*epsilon_c*beta_c**2)

corr_cosm_param = omega_m_false/omega_m_true*(h_false/h_true)**5
corr_growth = growth_true/growth_false

factor = corr_cosm_param*corr_growth
factor = 1

print "I have set factor =" , factor

Hubbleconstants_notcorrected = vrs_plus_hubbleflow/rs
Hubbleconstants_corrected = (vrs_plus_hubbleflow-factor*vrs_correction)/rs

print "The mean and variance of the non-corrected Hubbleconstants are", sp.mean(Hubbleconstants_notcorrected), sp.std(Hubbleconstants_notcorrected)
print "The mean and variance of the corrected Hubbleconstants are", sp.mean(Hubbleconstants_corrected), sp.std(Hubbleconstants_corrected)

f = open(path+min_dist+'Hubbleconstants_test_before_HS8.txt','w')
f.write("#observer\t not corrected\t corrected\n")
for obs,H_notcorr,H_corr in zip(range(len(observers)),Hubbleconstants_notcorrected,Hubbleconstants_corrected):
    f.write("%s\t%0.3f\t%0.3f\n" % (obs,H_notcorr,H_corr))


f.close()

