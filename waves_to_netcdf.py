
def CODAS_waves_to_netcdf(das_directory, output_directory,strain_factor,fs):
    import glob
    import xarray as xr
    import pandas as pd
    import numpy as np

    from datetime import date, datetime

    das_chanfolders = glob.glob(das_directory + '*')
    das_channels = [round(float(x.split('/')[-1])) for x in das_chanfolders]
    das_channels.sort()

    #get all times 
    das_timestr = sorted([x[60:75] for x in glob.glob(das_chanfolders[0]+'/*.ncdf')])
    das_times = sorted([datetime.strptime(x[60:75],'%Y%m%d_%H%M%S') for x in glob.glob(das_chanfolders[0]+'/*.ncdf')])
    
    
    f_cutoff = .3
    
    #create stacked arrays of outputs for all times, channels
    das_Hs = []
    das_Tp = []
    das_Te = []
    das_psd = []
    das_psd_corr = []
    for di,das_channel in enumerate(das_channels):
        Hs_inner = []
        Tp_inner = []
        Te_inner = []

        psd_inner = []
        psd_corr_inner = []
        for ti,das_time in enumerate(das_times):
            das_file = glob.glob(das_directory + str(das_channels[di]) + '/' + 'CODAS.D*__' + das_timestr[ti]+'.*__chn*'+str(das_channels[di])+'.ncdf')[0]
            ds_disk = xr.open_dataset(das_file)
            f_psd, ds_psd, ds_psd_corr, Tp_psd, Te_psd, Hs_psd = DAS_wave_conversion(ds_disk,f_cutoff,strain_factor,fs)
            Hs_inner.append(Hs_psd)
            Tp_inner.append(Tp_psd)
            Te_inner.append(Te_psd)

            psd_inner.append(ds_psd)
            psd_corr_inner.append(ds_psd_corr)
        das_Hs.append(Hs_inner)
        das_Tp.append(Tp_inner)
        das_Te.append(Te_inner)

        das_psd.append(psd_inner)
        das_psd_corr.append(psd_corr_inner)
        
    # define data with variable attributes
    data_Hs = {'Hs':(['channels','time'], das_Hs, 
                             {'units': 'm', 
                              'long_name':'significant wave height'})}
    data_Tp = {'Tp':(['channels','time'], das_Tp, 
                             {'units': 's', 
                              'long_name':'peak wave period'})}
    data_Te = {'Te':(['channels','time'], das_Te, 
                             {'units': 's', 
                              'long_name':'energy-weighted wave period'})}

    data_E = {'E':(['channels','time','frequency'], das_psd, 
                             {'units': 'm', 
                              'long_name':'energy spectrum'})}
    data_E_corr = {'E_corr':(['channels','time','frequency'], das_psd_corr, 
                             {'units': 'm', 
                              'long_name':'corrected energy spectrum'})}

    # define coordinates
    coords = {'time': (['time'], das_times),
              'channels': (['channels'], das_channels),
             'frequency': (['frequency'], f_psd)}

    # define global attributes
    attrs = {'creation_date':datetime.now(), 
             'author':'M Smith', 
             'email':'madisonmsmith@whoi.edu'}

    # create dataset
    ds_Hs = xr.Dataset(data_vars=data_Hs, 
                    coords=coords, 
                    attrs=attrs)
    ds_Tp = xr.Dataset(data_vars=data_Tp, 
                    coords=coords, 
                    attrs=attrs)
    ds_Tp = xr.Dataset(data_vars=data_Tp, 
                    coords=coords, 
                    attrs=attrs)

    ds_E = xr.Dataset(data_vars=data_E, 
                    coords=coords, 
                    attrs=attrs)
    ds_Ecorr = xr.Dataset(data_vars=data_E_corr, 
                    coords=coords, 
                    attrs=attrs)
    
    #merge datasets of each variable
    ds_DAS = xr.merge((ds_Hs,ds_Tp,ds_Tp,ds_E,ds_Ecorr))

    # SAVE AS NETCDF 
    ds_DAS.to_netcdf(output_directory +'/uw_'+datetime.strftime(das_times[0],'%Y-%m')+'_waveoutputs_strain'+str(strain_factor)+'.nc')