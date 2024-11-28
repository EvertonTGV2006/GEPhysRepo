
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv("2_1.csv");
print(data)

vIn = data[data.columns[0]];
r1In = data[data.columns[2]];
r2In = data[data.columns[2]];



# iIn = iIn * np.power(10.0, 0);

# vkVrt = vIn;

# absErrVin = 0.05
# absErriIn = 0.0005

# percErrVin = absErrVin / vIn
# percErrIin = absErriIn / iIn

# print(vkVrt)

# mu0 = 4 * np.pi * np.power(10.0 , -7)
# a = 0.147
# N = 124

# B = mu0 * N * iIn / (a * np.power(1.25, 1.5))

# R = 0.02

# percErrR = 0.0005 / R

# percErrBrSquared = 2 * (percErrIin + percErrR)

# brSquared = np.power(B * R, 2) / 2

print(vIn)

absErrV = 0.05
perErrV = absErrV / vIn * 0.5
vIn = vIn * 1000

print(vIn)

vkVrt = np.power(vIn, -0.5)

r1Err = 0.1
r2Err = 0.1
r1PerErr = r1Err / r1In
r2PerErr = r2Err / r2In

def lin(x,a,b):
    return a * x + b

parameters, covmat = curve_fit(lin, vkVrt, r1In, [1, 1])

print(parameters)

yFit = lin(vkVrt, parameters[0], parameters[1])

plt.errorbar(vkVrt, r1In, r1In * r1PerErr, perErrV * vkVrt, '.')
plt.plot(vkVrt, yFit, 'r')
plt.grid()
plt.xlabel("1 / V^0.5 / V^-0.5")
plt.ylabel("r1In")
parerr=np.sqrt(covmat.diagonal())
plt.annotate("y = " + str(parameters[0]) + " +- "  + str(parerr[0]) + " x + \n" + str(parameters[1]) + " +- " + str(parerr[1]), [0.014, 2.6])

print(covmat)
plt.show()


