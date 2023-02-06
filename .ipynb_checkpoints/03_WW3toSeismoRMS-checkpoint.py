from IPython.core.display import display, HTML
display(HTML("<style>.container { width:80% !important; }</style>")) 
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "last_expr"
import os
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset as netcdf_dataset
import netCDF4
import numpy as np
from scipy import interpolate
import scipy.integrate
from obspy.geodetics.base import gps2dist_azimuth
from obspy import UTCDateTime
from obspy.geodetics.base import kilometer2degrees, locations2degrees
import numpy.matlib
import warnings
warnings.filterwarnings('ignore')
from cartopy import config
import cartopy.crs as ccrs
os.environ["CARTOPY_USER_BACKGROUNDS"] = "BG/"

import ww32seismo
from ww32seismo import *


# Import configurations and Settings

import dynamic_yaml
with open("config.yml", 'r') as f:
    configs = dynamic_yaml.load(f)

target = configs.params.station
target_lat = configs.params.station_lat
target_lon = configs.params.station_lon
rhos = configs.params.rhos
beta = configs.params.beta
Rg = configs.params.Rg
Q = configs.params.Q
Re = 4.0e7/(2*np.pi)
depth_file = configs.files.depth

#  needed files: Rayleigh_source.txt, depth file, and 2 models with and without reflection /!\

dataset = netcdf_dataset(r"{}".format(depth_file))
dpt = pd.DataFrame(np.asarray(dataset.variables["dpt"])[50,:,:], columns=dataset.variables["longitude"], index=dataset.variables["latitude"])
dpt[dpt==-32767] *= 0.0
dpt[dpt<=0.0] = 0.0

if not os.path.isdir("DATA"):
    os.mkdir("DATA")
if not os.path.isdir("FIGURES"):
    os.mkdir("FIGURES")    


# # Plot depth and distance

lats, lons, distance_df = ww32seismo.get_distance(configs, dataset, dpt, plot=False)

df = pd.read_csv(r"DATA/Rayleigh_source.txt", header=None, delim_whitespace=True, index_col=0)
df.index *= np.pi
df = df.fillna(0.0)
C_base = (df[:8]**2).sum(axis=1)
C_base.at[C_base.index[-1]+0.01] = 0.0
C_base.at[-1.0] = 0.0
C_base.at[20.0] = 0.0
C_base = C_base.sort_index()

Cf = interpolate.interp1d(C_base.index, C_base)


# In[7]:


from ww32seismo import *
dfF_fs = get_ww3(configs, 10, lats, lons, Re, dpt, Cf, distance_df, plot=False)

fig = plt.figure(figsize=(16,4), facecolor="w")
cmap = plt.get_cmap('viridis')
psd = 10* np.log10(dfF_fs)
plt.pcolormesh(dfF_fs.columns, 1./dfF_fs.index, psd, cmap=cmap, vmin = -150, vmax =-110)
plt.colorbar().set_label("Amplitude (dB)")
plt.ylabel("Period (s)")
plt.yscale('log')
fig.autofmt_xdate()
#plt.ylim(10e-2,10e2)
plt.show()

fig = plt.figure(figsize=(16,5), facecolor="w")
integ = np.sqrt(scipy.integrate.trapz(dfF_fs.fillna(0), dfF_fs.index, axis=0))
plt.plot(dfF_fs.columns, integ)
plt.ylabel("Amplitude")
fig.autofmt_xdate()
plt.xlim(dfF_fs.columns[0],dfF_fs.columns[-1])
plt.show()

