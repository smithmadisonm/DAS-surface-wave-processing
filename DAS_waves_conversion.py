def DAS_wave_conversion(das_data,f_cutoff,strain_factor,fs):
    #function to use simple pwelch to estimate wave spectra and bulk wave parameters for DAS data
    #assume some arbitrary strain factor that results in decent looking results, for now...
    
    import pandas as pd
    import numpy as np
    from scipy import signal
    
    #pwelch - defualt is 50% overlap
    window = 128
    nfft = 256
    f_psd, ds_psd = signal.welch((das_data)/strain_factor,fs=fs,nfft=nfft,nperseg=window)#,window=[128]
    
    #depth attenuation correction (using depth from netcdf)
    depth = 1
    k = (2*np.pi*f_psd)**2 / 9.8
    attenuation = np.exp(k*depth)
    attenuation = attenuation**2; # square for energy 
    ds_psd_corr = ds_psd*attenuation
    
    #calculate bulk wave characteristics
    max_i = np.argmax(ds_psd)
    Tp = 1/(f_psd[max_i])
    
    psd_fwaves = ((f_psd > 0.04) & (f_psd < f_cutoff))
    fe = ((ds_psd_corr[psd_fwaves] * f_psd[psd_fwaves]) /ds_psd_corr[psd_fwaves].sum() ).sum() #(f*E)/E
    Te = 1/fe
    
    bandwidth = (f_psd[1::] - f_psd[0:-1]).mean()
    Hs = 4*np.sqrt( ds_psd_corr[psd_fwaves].sum() * bandwidth ) 
    
    return f_psd, ds_psd, ds_psd_corr, Tp, Te, Hs



