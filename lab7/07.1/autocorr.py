import matplotlib.pyplot as plot
import numpy as np
import pandas as pd

press=np.loadtxt("gas/output.pres.0", usecols=(1), unpack='true')
ener=np.loadtxt("gas/output.epot.0", usecols=(1), unpack='true')
pd.plotting.autocorrelation_plot(press)
pd.plotting.autocorrelation_plot(ener)

plot.xlim(0,2000)
plot.show()
