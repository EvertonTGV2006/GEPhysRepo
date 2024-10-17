
from email.utils import collapse_rfc2231_value
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def exp(x, a, b, c):
    return a * np.exp(-b * x) + c

def lin(x, a, b):
    return a * x + b

def poly(x, a, b, c):
    return a * (x **2) + b * x + c

def ex_4_1(df):
    plt.figure(4)
    x = df[df.columns[0]].values
    x_err = df[df.columns[1]].values
    y = df[df.columns[2]].values
    y_err = df[df.columns[3]].values
    plt.errorbar(x, y, y_err, fmt = '.', label = 'Raw Data', capsize=5)
    par_lin, cov_mat = curve_fit(lin, x, y, sigma=y_err, absolute_sigma=True)
    par_err = np.sqrt(cov_mat.diagonal())
    a_fit = par_lin[0]
    a_err = par_err[0]
    b_fit = par_lin[1]
    b_err = par_lin[1]
    plt.plot(x, lin(x, a_fit, b_fit), label = 'Fitted Line')

    plt.ylabel('Velocity / km/s')
    plt.xlabel('Distance / M-Parsecs')
    plt.title("Calculation of Hubble's Constant")

    plt.annotate('Fitted Line: A * x + B:\nA: ' + str(a_fit) + ' +- ' + str(a_err) + '\n' + 'B: ' + str(b_fit) + ' +- ' + str(b_err) + '\nSo Hubble\'s Constant = ' + str(a_fit) + '\n+- ' + str(a_err) + ' km/s / M-parsec', [10, 25000])
    plt.savefig('ex5_2.png', dpi=200)

    plt.show()

def ex_5_1(df):
    #Exponential Fitting
    # 5.1 exp fitting
    x = df['time/s'].values
    y = df['counts'].values

    
    
    A = 1200 #guesses of A
    B = 2
    C = 0.5

    parameters, covariant_matrix = curve_fit(exp, x, y, p0=(A, B, C))

    parameter_uncert = np.sqrt(covariant_matrix.diagonal())


    a_fit = parameters[0]
    a_err = parameter_uncert[0]
    b_fit = parameters[1]
    b_err = parameter_uncert[1]
    c_fit = parameters[2]
    c_err = parameter_uncert[2]
    plt.figure(1);
    plt.plot(x, y, '.', label = 'Raw Data')
    plt.plot(x, exp(x, a_fit, b_fit, c_fit), 'r', label = 'Fitted Line')
    plt.xlabel=df.columns[0]
    plt.ylabel = df.columns[1]


    plt.xlim(left = 0)
    plt.ylim(bottom = 400, top = 1300)
    plt.legend()
       
    print('y = ' + str(a_fit) + ' e ^(-' + str(b_fit) + 'x) + ' + str(c_fit))
    print('A: ' + str(a_fit) + ' +- ' + str(a_err))
    print('B: ' + str(b_fit) + ' +- ' + str(b_err))
    print('C: ' + str(c_fit) + ' +- ' + str(c_err))

    plt.annotate('Fitted Line: A * e^(-B * x) + C \nA: ' + str(a_fit) + ' +- ' + str(a_err) + '\n' + 'B: ' + str(b_fit) + ' +- ' + str(b_err) + '\n' + 'C: ' + str(c_fit) + ' +- ' + str(c_err), [0.5, 1100])

    plt.savefig('ex5_1.png', dpi = 200)

    return a_fit, b_fit, c_fit

def ex_5_2(df, a_fit, b_fit, c_fit):
    x = df['time/s'].values
    y = df['counts'].values
    df['counts corrected'] = df['counts'].values - c_fit
    df['counts log'] = np.log(df['counts corrected'])

    plt.figure(2)
    plt.plot(x, df['counts log'], '.', label = 'Raw Data')

    A = -0.5
    B = 7

    par_lin, cov_lin = curve_fit(lin, x, df['counts log'], p0 = (A, B))
    a2_fit = par_lin[0]
    b2_fit = par_lin[1]
    par_err = np.sqrt(cov_lin.diagonal())
    a2_err = par_err[0]
    b2_err = par_err[1]
    plt.plot(x, lin(x, a2_fit, b2_fit), label = 'Fitted Line')
    plt.xlabel=df.columns[0]
    plt.ylabel = 'Logarithmic Corrected Counts'


    plt.xlim(left = 0)
    plt.ylim(bottom = 5, top = 7)
    plt.legend()

    plt.annotate('Fitted Line: A * x + B:\nA: ' + str(a2_fit) + ' +- ' + str(a2_err) + '\n' + 'B: ' + str(b2_fit) + ' +- ' + str(b2_err), [0.5, 6.7])
    plt.savefig('ex5_2.png', dpi=200)


def ex_6_1():
    x = [3, 6, 9, 12, 15, 18]
    y = [4.2, 22.1, 51.6, 91.8, 141.6, 200.6]
    df = pd.DataFrame({"time of flight (ms)":x, "displacement (pixels)":y})
    df["time of flight (s)"] = df["time of flight (ms)"]*1e-3 #convert to ms to seconds

    par_poly, cov_poly = curve_fit(poly, x, y)
    par_poly_err = np.sqrt(cov_poly.diagonal())
    a3_fit = par_poly[0]
    b3_fit = par_poly[1]
    c3_fit = par_poly[2]
    a3_err = par_poly_err[0]
    b3_err = par_poly_err[1]
    c3_err = par_poly_err[2]

    plt.figure(3)
    plt.plot(x, y, '.', label = 'Raw Data')
    y_fit = [];
    for i in range(0, len(x)):
        y_fit.append(poly(x[i], a3_fit, b3_fit, c3_fit))
    plt.plot(x, y_fit, 'r', label = 'Fitted Polynomial')


    plt.xlabel='time / s'
    plt.ylabel = 'Logarithmic Corrected Counts'


    plt.xlim(left = 0)
    plt.ylim(bottom = 0, top = 225)

    plt.legend()



    plt.annotate('Fitted Line: A x^2 + B x + C \nA: ' + str(a3_fit) + ' +- ' + str(a3_err) + '\n' + 'B: ' + str(b3_fit) + ' +- ' + str(b3_err) + '\n' + 'C: ' + str(c3_fit) + ' +- ' + str(c3_err), [0.2, 185])
    plt.savefig('Ex_6_1.png', dpi = 200);


if __name__ == '__main__':
    df=pd.read_csv('hubble.csv')
    ex_4_1(df)

    df = pd.read_csv('protein_fluorescence.csv')

    a_fit, b_fit, c_fit = ex_5_1(df)
    ex_5_2(df, a_fit, b_fit, c_fit)
    ex_6_1()

    #5.2
    #make new col for corrected counts


    #6 now , input data by hand as a small dataset





    plt.show();




