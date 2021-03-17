import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import KernelDensity


df = pd.read_csv('split_y.csv')

X = df.split.values

X_plot = np.linspace(min(X),max(X),1000)

X_plot = X_plot.reshape(-1,1)

X = X.reshape(-1,1)

kde = KernelDensity(kernel='gaussian', bandwidth=5).fit(X)

log_dens = kde.score_samples(X_plot)

#plt.plot(X_plot,np.exp(log_dens))

data = np.exp(log_dens)

xout = X_plot.flatten()

df = pd.DataFrame({'sse':xout,'kde':data})
df.to_csv('kde_data.csv')
