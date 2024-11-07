import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def lin(x, a, b):
    return x * a + b;

x = np.array([30, 60, 90, 120, 150, 180, 210])
y = np.array([4.72, 18.44, 46.50, 88.65, 145.0, 212.1, 294.2])

lnx = np.log(x)
lny = np.log(y)

a, b = curve_fit(lin, lnx, lny)


print("A: " + str(a[0]) + " B: " + str(a[1]))


plt.xlabel("Log Number of Turns")
plt.ylabel("Log Inductance")
plt.xlim(3, 6)
plt.ylim(0, 10)
plt.plot(lnx, lny, '.')
plt.plot(lnx, lin(lnx, a[0], a[1]))
plt.show()

print(np.pi)


