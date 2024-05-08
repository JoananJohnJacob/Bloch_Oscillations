import csv
import numpy as np
import matplotlib.pyplot as plt
import os

def readComplexNumbersFromCSV(filename):
    complex_numbers = []
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            real_part = float(row[0])
            imag_part = float(row[1])
            complex_numbers.append(complex(real_part, imag_part))
    return np.array(complex_numbers)

def readVectorFromCSV(filename):
    numbers=[]
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            number = float(row[0])
            numbers.append(number)
    return np.array(numbers)

P_total = readComplexNumbersFromCSV("Total Polarization.csv")
t = readVectorFromCSV("Time.csv")
E = readVectorFromCSV("Electric Field.csv")

############################################# Fourier Transforms ##############################################

number_column = len(t)
dt = t[1]-t[0]
h = 4.135667696e-15 #Planck's constant in eV s
freq = np.fft.fftfreq(number_column,dt)
Energy_vals = h*freq

E_w = np.fft.ifft(E) #Fourier Transform of Electric field
P_w = np.fft.ifft(P_total)  #Fourier Transform of Total Polarization
alpha_w = P_w/E_w
alpha_w = alpha_w.imag # Absorption
min_index = np.argmin(E_w)
min_value = E_w[min_index]


###################################################################################################################


plt.figure(1)


plt.subplot(3,2,1)
plt.plot(t,E)
plt.xlim([5e-15,25e-15])
plt.grid("True")
plt.xlabel("Time")
plt.ylabel("Electric field")

plt.subplot(3,2,2)
plt.plot(t,P_total)
plt.xlabel("Time")
plt.ylabel("Total Polarization")


plt.subplot(3,2,3)
plt.plot(Energy_vals,abs(E_w))
plt.xlabel('Energy')
plt.xlim([0,5])
plt.ylabel('Electric Field')
plt.grid('True')
# Plot of Electric Field

plt.subplot(3,2,4)
plt.plot(t,P_total,label='Total Polarization')

plt.legend()
plt.xlabel('Time')
plt.ylabel('Total Polarization')
# Plot of Polarization

plt.subplot(3,2,5)
plt.plot(Energy_vals,abs(P_w),label="|P(w)|")
plt.plot(Energy_vals,P_w.real,label="Re{P(w)}")
plt.plot(Energy_vals,P_w.imag,label="Img{P(w)}")

plt.xlabel('Energy in eV')
plt.ylabel('P($\Omega$)')
plt.title('Energy vs P($\omega$)')
plt.xlim([0,5])
# plt.ylim([-2e-39,6e-39])
plt.grid('True')
plt.legend()
#plot of P(w)
plt.subplot(3,2,6)
plt.plot(Energy_vals,alpha_w,'m-',label='Absorption')
plt.xlabel('Energy')
plt.ylabel('Absorption')
plt.xlim([1,3])
plt.ylim([0,2e-34]) 
plt.grid('True')

plt.tight_layout()
""" 
current_dir = os.getcwd()

plt.savefig(os.path.join(current_dir,'Results.pdf'),bbox_inches='tight')
 """
plt.show()

