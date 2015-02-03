import scipy as sp
import matplotlib.pyplot as plt

data0 = sp.loadtxt('../../cases/Planck512/snapshot/powerspectra.txt')

#data1 = sp.loadtxt('../../cases/Planck512/powerspectra_nside16_lmax16_number_of_SNe7125.txt')
#data2 = sp.loadtxt('../../cases/Planck512/powerspectra_nside16_lmax20_number_of_SNe7125.txt')
#data3 = sp.loadtxt('../../cases/Planck512/powerspectra_nside16_lmax32_number_of_SNe7125.txt')
#
#
#data4 = sp.loadtxt('../../cases/Planck512/powerspectra_nside16_lmax32_number_of_SNe200.txt')
#data5 = sp.loadtxt('../../cases/Planck512/powerspectra_nside16_lmax32_number_of_SNe2000.txt')
#data6 = sp.loadtxt('../../cases/Planck512/powerspectra_nside16_lmax32_number_of_SNeall.txt')
#
#data7 = sp.loadtxt('../../cases/Planck512/powerspectra_nside8_lmax16_number_of_SNe7125.txt')
#data8 = sp.loadtxt('../../cases/Planck512/powerspectra_nside16_lmax64_number_of_SNe7125.txt')
#data9 = sp.loadtxt('../../cases/Planck512/powerspectra_nside32_lmax64_number_of_SNe7125.txt')
#data10 = sp.loadtxt('../../cases/Planck512/powerspectra_nside64_lmax128_number_of_SNe7125.txt')


#for i, data in enumerate([data1, data2, data3, data4, data5, data6, data7, data8, data9, data10]):
for i, data in enumerate([data0]):
    ls = data[0,2:]
    Cls = data[1:,2:]
    
    mean_Cls = sp.mean(Cls,axis=0)
    scaled_Cls = sp.array(sp.sqrt(mean_Cls*ls*(ls+1)))
    plt.figure()
    plt.plot(ls,scaled_Cls)
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlim([0,20])
    plt.ylim([0,800])    
    titel = 'datasaet nummer' + str(i+1)
    plt.title(titel)
