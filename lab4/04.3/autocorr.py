import matplotlib.pyplot as plot
import numpy as np
import pandas as pd

data=np.loadtxt("output_press.dat", usecols=(0), delimiter=' ', unpack='true')
pd.plotting.autocorrelation_plot(data);
plot.show()
