import matplotlib.pyplot as plot
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from statsmodels.tsa.stattools import acf

def f(x,a,b):  # defining the fitting function
    return 1 * np.exp(-b*x)

#x = np.linspace(1,15,16,endpoint=True)

data=np.loadtxt("output_press.dat", usecols=(0), delimiter=' ', unpack='true')
#y=np.array(acf(data,unbiased=True,nlags=15))
pd.plotting.autocorrelation_plot(data)

#p_opt, p_cov = curve_fit(f, x, y)
#y_fit = f(x,p_opt[0],p_opt[1])
#print(p_opt[0],p_opt[1])
#plot.plot(x,y_fit) # plotting fitted function
#plot.plot(x,y)
plot.xlim(0,140)
plot.show()
