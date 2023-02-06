import ftplib
import os
import tqdm
import dynamic_yaml
import numpy as np


    
def download_p2l():
    """
    Downloads the P2L files directly from the IFREMER ftp for the year in the config.yml file
    """
    with open("config.yml", 'r') as f:
        configs = dynamic_yaml.load(f)
    # Fill Required Information
    HOSTNAME = "ftp.ifremer.fr"
    # Connect FTP Server
    ftp_server = ftplib.FTP(HOSTNAME)
    ftp_server.login() # Anonymous login
    # force UTF-8 encoding
    ftp_server.encoding = "utf-8"
    #from config
    sismo_path = configs.download.sismo_path
    year = configs.download.year
    p2l_save = configs.download.p2l_save
    if not os.path.isdir(p2l_save):
        os.mkdir(p2l_save)

    ref_paths = []

    try:
        ref_paths = ftp_server.nlst(sismo_path) # Get initial work directories with and without reflection
    except ftplib.error_perm as resp:
        if str(resp) == "550 No files found":
            print("No files in this directory")
        else:
            raise

    for p in ref_paths:
        print(p)
        f_list = ftp_server.nlst(os.path.join(p,"{}".format(year), "FIELD_NC"))
        p2l_list = [i for i in f_list if "_p2l.nc" in i] # Only list the p2l files
        pbar = tqdm.tqdm(p2l_list, desc="Downloading {}".format(p))
        for p2l in pbar:
            fn = p2l.split("/SISMO/")[1].split("FIELD_NC/")[1]
            save_dir = os.path.join(p2l_save, p2l.split("/SISMO/")[1].split("/LOPS")[0])
            if not os.path.isdir(save_dir):
                os.makedirs(save_dir)
            fpath = os.path.join(save_dir, fn)
            if not os.path.isfile(fpath):
                ftp_server.retrbinary('RETR {}'.format(p2l),open(fpath, 'wb').write)
    ftp_server.quit()
    
def plot_spec(dfF_fs, station):
    """
    Plots the spectroram for the station if the corresponding dataframe if already available
    """
    plt.rcParams['figure.figsize'] = (16,6)
    plt.rcParams['axes.facecolor'] = "w"
    fig, ax = plt.subplots()


    cmap = plt.get_cmap('viridis')
    cmin=-160
    cmax=-110
    ymin=1
    ymax=12

    psd = 10* np.log10(dfF_fs)
    plt.pcolormesh(dfF_fs.columns, 1./dfF_fs.index, psd, cmap=cmap, vmin = cmin, vmax =cmax)
#    plt.axvline(pd.to_datetime("2021-10-25"), color="w", ls="--")
#    plt.axvline(pd.to_datetime("2021-11-05"), color="w", ls="--")
    cb = plt.colorbar(ax=ax).set_label("Amplitude [$m^2/s^4/Hz$] [dB]")
    plt.ylabel("Period (s)")
    plt.yscale('log')
    fig.autofmt_xdate()
    plt.ylim(0.1,ymax)
    plt.title(station)
    plt.show()

def plot_rms(dfF_fs, station):
    """
    Plots the seismic RMS for the station if the corresponding dataframe if already available
    """    
    fig = plt.figure(figsize=(16,6), facecolor="w")
    integ = np.sqrt(scipy.integrate.trapz(dfF_fs.fillna(0), dfF_fs.index, axis=0))
    plt.plot(dfF_fs.columns, integ)
    plt.ylabel("Amplitude")
    fig.autofmt_xdate()
    plt.title(station)
    plt.xlim(dfF_fs.columns[0],dfF_fs.columns[-1])
    plt.show()    
    
    
    
def dispNewtonTH(f,dep):
    """inverts the linear dispersion relation (2*pi*f)^2=g*k*tanh(k*dep) to get
    k from f and dep. 2 Arguments: f and dep.
    """
    eps = 0.000001
    g = 9.81
    sig = 2.*np.pi*f
    Y = dep*sig**2./g
    X = np.sqrt(Y)
    I=1
    F = np.ones(dep.shape)
    
    while np.abs(np.max(F)) > eps:
        H = np.tanh(X)
        F = Y-X*H
        FD = -H-X/np.cosh(X)**2
        X -= F/FD
        F = F.values
    return X/dep

def compute_depth_correction(f, dep):
    wnum = dispNewtonTH(f, dep)
    kd = wnum*dep
    depth_correction = (np.tanh(kd))**2.0*(1.0 + 2.0*(kd/np.sinh(2*kd)))
    return depth_correction

def get_omegaoverbeta(fs, dpt):
    omega = 2 * np.pi * fs * 2
    omehoverbeta=omega * dpt / beta     # omega*h/beta
    omehoverbeta[omehoverbeta>=10] = 10.0
    omehoverbeta[omehoverbeta<=0] = 0.0
    return omehoverbeta

def get_C(fs, dpt):
    """Returns the value of ^C for fs
    """
    ob = omegaoverbeta(fs, dpt)
    return Cf(ob)