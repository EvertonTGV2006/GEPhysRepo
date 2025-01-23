
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv("1_1.csv");
print(data)

vIn = data[data.columns[0]];
r1In = data[data.columns[1]];
#r2In = data[data.columns[2]];



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
vIn = vIn


print(vIn)

vkVrt = np.power(vIn, -0.5)

r1Err = 0.002
r2Err = 0.05
r1PerErr = r1Err / r1In
#r2PerErr = r2Err / r2In

vkVrt = vIn

r1In=np.pow(r1In,2)

def lin(x,a,b):
    return a * x + b

parameters, covmat = curve_fit(lin, vkVrt, r1In, [1, 1])

print(parameters)

yFit = lin(vkVrt, parameters[0], parameters[1])

plt.errorbar(vkVrt, r1In, r1In * r1PerErr, perErrV * vkVrt, '.',color='tab:blue')
plt.plot(vkVrt, yFit, color='tab:cyan')
#parameters2, covmat2 = curve_fit(lin, vkVrt, r2In, [1, 1])

#print(parameters2)

#yFit = lin(vkVrt, parameters2[0], parameters2[1])

#plt.errorbar(vkVrt, r2In, r2In * r1PerErr, perErrV * vkVrt, '.',color='tab:red')
#plt.plot(vkVrt, yFit, color='tab:orange')
plt.grid()
plt.xlabel("Accelerating Voltage / V")
plt.ylabel("I^2 / A^2")
plt.title("Figure 2")
#plt.xlim(0.012, 0.021)
#plt.ylim(0.75, 2.75)
parerr=np.sqrt(covmat.diagonal())
#parerr2=np.sqrt(covmat2.diagonal())
plt.annotate("y = (" + "{:.2g}".format(parameters[0]) + "+-"  + "{:.1g}".format(parerr[0]) + ") * x", [200, 11])
#plt.annotate("y = (" + "{:.0f}".format(parameters2[0]) + "+-"  + "{:.0f}".format(parerr2[0]) + ") * x + (" + "{:.2f}".format(parameters2[1]) + "+-" + "{:.0g}".format(parerr2[1]) + ")", [0.014, 2.6])

print(covmat)
plt.show()


