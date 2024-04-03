function RBRresults2NC(RBR_in,filename)
%
% creates a netCDF file using existing RBR structure and writes it into 'filename'
% (must include .nc)
%
%   >> RBR2NC(RBR, filename)
%
% Use SWIFT ID to determine v3 or v4 (sensors at different depths) 
% and skip substructures that are not supported yet
%
% original by L. Hosekova in 2020 for SWIFT data structures
%   Edited on May 1, 2020 by Suneil Iyer for ATOMIC with compliant names
%   Feb 2023 by J. Thomson for SASSIE and general use
%   Apr 2024, J. Thomson adapted to bottom moorings

RBR = RBR_in;

if isfield(RBR,'wavehistogram')
    RBR=rmfield(RBR,'wavehistogram');
end

RBR=rmfield(RBR,'IGwaveheight');
RBR=rmfield(RBR,'IGperiod');

if isfield(RBR,'salinity') && length(RBR(1).salinity)>1
    for si=1:length(RBR)
        RBR(si).salinity = nanmean(RBR(si).salinity);
        RBR(si).watertemp = nanmean(RBR(si).watertemp);
    end
end

%% loading variables
% extract dimension sizes: time, freq, z, zHR (if available)

ncid=netcdf.create(filename,'CLOBBER');
t_dim=netcdf.defDim(ncid,'time', length(RBR));
full_names=fieldnames(RBR);

if isfield(RBR,'wavespectra') && min(RBR(1).wavespectra.freq)>0
    f_dim = netcdf.defDim(ncid,'freq', length(RBR(1).wavespectra.freq));
    spec_names=fieldnames(RBR(1).wavespectra);
end
if isfield(RBR,'uplooking')
    z_dim = netcdf.defDim(ncid,'z', length(RBR(1).uplooking.z));
    z_names = fieldnames(RBR(1).uplooking);
end
if isfield(RBR,'downlooking')
    z_dim = netcdf.defDim(ncid,'z', length(RBR(1).downlooking.z));
    z_names = fieldnames(RBR(1).downlooking);
end
if isfield(RBR,'signature')
    sig_names = fieldnames(RBR(1).signature)
    if isfield(RBR(1).signature,'HRprofile')
        zHR_dim = netcdf.defDim(ncid,'zHR', length(RBR(1).signature.HRprofile.z));
        zHR_names = fieldnames(RBR(1).signature.HRprofile);
    end
    if isfield(RBR(1).signature,'profile')
        z_dim = netcdf.defDim(ncid,'z', length(RBR(1).signature.profile.z));
        z_names = fieldnames(RBR(1).signature.profile);
    end
end



j=1;
for i=1:length(full_names), 
    if ~strcmp(full_names{i},'ID') && ~strcmp(full_names{i},'date')
        if strcmp(full_names{i},'signature')
            for t=1:length(RBR)
                for iz=1:length(z_names)
                    eval(strcat('S.signature.profile.',z_names{iz},'(t,:)=RBR(t).signature.profile.',z_names{iz},'(:)'))
                end
                for iz=1:length(zHR_names)
                    eval(strcat('S.signature.HRprofile.',zHR_names{iz},'HR(t,:)=RBR(t).signature.HRprofile.',zHR_names{iz},'(:)'))
                end
            end
        elseif strcmp(full_names{i},'time')
            S.time= [RBR.time]-datenum(1970,1,1,0,0,0);
        else
            eval(strcat('S.',full_names{i},'=[RBR.',full_names{i},']')); % errors here if check factor in some but not all
        end
        names{j} = full_names{i};
        j = j+1;
    end
end


%% creating netcdf variables

for i=1:length(names)
    if strcmp(names{i},'wavespectra')
        for j=1:length(spec_names)
            if strcmp(spec_names{j},'freq')
                eval(strcat(spec_names{j},'_id = netcdf.defVar(ncid,''',spec_names{j},''',''NC_DOUBLE'',[f_dim])'));
            else
                eval(strcat(spec_names{j},'_id = netcdf.defVar(ncid,''',spec_names{j},''',''NC_DOUBLE'',[t_dim f_dim])'));
            end
        end
    elseif strcmp(names{i},'uplooking') || strcmp(names{i},'downlooking')
        for j=1:length(z_names)
            if strcmp(z_names{j},'z')
                eval(strcat(z_names{j},'_id = netcdf.defVar(ncid,''',z_names{j},''',''NC_DOUBLE'',[z_dim])'));
                %             elseif strcmp(z_names{j},'tkedissipationrate')
                %                 eval(strcat(z_names{j},'_id = netcdf.defVar(ncid,''',z_names{j},''',''NC_DOUBLE'',[t_dim])'));
            else
                eval(strcat(z_names{j},'_id = netcdf.defVar(ncid,''',z_names{j},''',''NC_DOUBLE'',[t_dim z_dim])'));
            end
        end
    elseif strcmp(names{i},'signature')
        for j=1:length(z_names)
            if strcmp(z_names{j},'altimeter')
                eval(strcat(z_names{j},'_id = netcdf.defVar(ncid,''',z_names{j},''',''NC_DOUBLE'',[t_dim])'));
            elseif strcmp(z_names{j},'z')
                eval(strcat(z_names{j},'_id = netcdf.defVar(ncid,''',z_names{j},''',''NC_DOUBLE'',[z_dim])'));
            else
                eval(strcat(z_names{j},'_id = netcdf.defVar(ncid,''',z_names{j},''',''NC_DOUBLE'',[t_dim z_dim])'));
            end
        end
        for j=1:length(zHR_names)
            if strcmp(zHR_names{j},'z')
                zHR_id = netcdf.defVar(ncid,'zHR','NC_DOUBLE',[zHR_dim]);
            else
                eval(strcat(zHR_names{j},'HR_id = netcdf.defVar(ncid,''',zHR_names{j},'HR'',''NC_DOUBLE'',[t_dim zHR_dim])'));
            end
        end
        %edit variable names to CF convention names - do this for all vars
        %with different names than the RBR defaults
    elseif strcmp(names{i},'lon')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''lon'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'lat')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''lat'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'watertemp')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_water_temperature'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'watertemp_d2')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_water_temperature_at_depth'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'airtemp')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''air_temperature'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'salinity')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_water_salinity'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'salinity_d2')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_water_salinity_at_depth'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'qa')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''specific_humidity'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'peakwavedirT')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_surface_wave_from_direction_at_variance_spectral_density_maximum'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'peakwaveperiod')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_surface_wave_period_at_variance_spectral_density_maximum'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'centroidwaveperiod')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_surface_wave_mean_period'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'sigwaveheight')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_surface_wave_significant_height'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'mss')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_surface_wave_mean_square_slope_normalized_by_frequency_width'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'ustar')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''friction_velocity_in_air_from_wave_spectra'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'airtempstddev')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''air_temperature_stddev'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'windspd')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''wind_speed'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'windspdstddev')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''wind_speed_stddev'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'winddirT')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''wind_direction'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'winddirTstddev')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''wind_direction_stddev'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'driftspd')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''drift_speed'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'driftdirT')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''drift_direction'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'relhumidity')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''relative_humidity'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'qsea')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''sea_surface_saturation_specific_humidity'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'qair')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''specific_humidity'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'relhumiditystddev')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''relative_humidity_stddev'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'airpres')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''air_pressure'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'airpresstddev')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''air_pressure_stddev'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'flag_values_watertemp')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''flag_values_watertemp'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'flag_values_airtemp')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''flag_values_airtemp'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'flag_values_humidity')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''flag_values_humidity'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'flag_values_windpsd')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''flag_values_windspd'',''NC_DOUBLE'',[t_dim])'));
    elseif strcmp(names{i},'flag_values_salinity')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''flag_values_salinity'',''NC_DOUBLE'',[t_dim])'));

    elseif ~strcmp(names{i},'ID')
        eval(strcat(names{i},'_id = netcdf.defVar(ncid,''',names{i},''',''NC_DOUBLE'',[t_dim])'));

    end
end


% define some global attributes
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid,varid,'creator','Jim Thomson (APL-UW) and Madison Smith (WHOI)');
netcdf.putAtt(ncid,varid,'please_acknowledge:','investigators above');
netcdf.putAtt(ncid,varid,'institution','Applied Physics Laboratory at the University of Washington (APL-UW) and Woods Hole Oceanographic Inst. (WHOI)');
netcdf.putAtt(ncid,varid,'contact_email_1','jthomson@apl.washington.edu');
netcdf.putAtt(ncid,varid,'description',['Seafloor mooring measuring pressure and temperature.  Pressure data are processed for surface wave statistics.']);
netcdf.putAtt(ncid,varid,'comment','Wave spectra were corrected for depth attenuation following linear wave theory.  The reported spectra are surface wave energy spectra.  The spectra are truncated, such that higher frequencies with noise contamination are removed.');
netcdf.putAtt(ncid,varid,'reference1','https://github.com/smithmadisonm/DAS-surface-wave-processing');
netcdf.putAtt(ncid,varid,'level','Version 1');
netcdf.putAtt(ncid,varid,'history','Version 1');
netcdf.putAtt(ncid,varid,'missing_data_flag','-999');


netcdf.endDef(ncid);
%% filling them with values

for i=1:length(names)
    if strcmp(names{i},'wavespectra')
        for j=1:length(spec_names)
            if strcmp(spec_names{j},'freq')
                netcdf.putVar(ncid,freq_id, S.wavespectra(1).freq);
            else
                eval(strcat('netcdf.putVar(ncid,',spec_names{j},'_id, [S.wavespectra.',spec_names{j},']'')'));
            end
        end
    elseif strcmp(names{i},'uplooking')
        for j=1:length(z_names)
            if strcmp(z_names{j},'z')
                netcdf.putVar(ncid,z_id, S.uplooking(1).z);
            else
                eval(strcat('netcdf.putVar(ncid,',z_names{j},'_id, [S.uplooking.',z_names{j},'])'));
            end
        end
    elseif strcmp(names{i},'downlooking')
        for j=1:length(z_names)
            if strcmp(z_names{j},'z')
                netcdf.putVar(ncid,z_id, S.downlooking(1).z);
            else
                eval(strcat('netcdf.putVar(ncid,',z_names{j},'_id, [S.downlooking.',z_names{j},']'')'));
            end
        end
    elseif strcmp(names{i},'signature')
        for j=1:length(z_names)
            if strcmp(z_names{j},'z')
                netcdf.putVar(ncid,z_id, S.signature.profile.z(1,:));
            else
                eval(strcat('netcdf.putVar(ncid,',z_names{j},'_id, [S.signature.profile.',z_names{j},']'')'));
            end
        end
        for j=1:length(zHR_names)
            if strcmp(zHR_names{j},'z')
                netcdf.putVar(ncid,zHR_id, S.signature.HRprofile.zHR(1,:));
            else
                eval(strcat('netcdf.putVar(ncid,',zHR_names{j},'HR_id, [S.signature.HRprofile.',zHR_names{j},'HR]'')'));
            end
        end
    else

        eval(strcat('netcdf.putVar(ncid,',names{i},'_id, S.',names{i},')'));

    end
end


netcdf.close(ncid)


%% units and descriptions
for i=1:length(names)
    if strcmp(names{i},'wavespectra')
        for j=1:length(spec_names)
            if strcmp(spec_names(j),'energy')
                ncwriteatt(filename,'energy','units','m^2/Hz')
                ncwriteatt(filename,'energy','long_name','wave energy spectral density as a function of frequency')
                ncwriteatt(filename,'energy','standard_name','sea_surface_wave_variance_spectral_density')
                ncwriteatt(filename,'energy','instrument','RBR Duo pressure and temperature')
                ncwriteatt(filename,'energy','method','https://github.com/smithmadisonm/DAS-surface-wave-processing/blob/main/processRBRduo.m')
            end
            if strcmp(spec_names(j),'freq')
                ncwriteatt(filename,'freq','units','Hz')
                ncwriteatt(filename,'freq','long_name','spectral frequencies')
                ncwriteatt(filename,'freq','standard_name','wave_frequency')
                ncwriteatt(filename,'energy','instrument','RBR Duo pressure and temperature')
                ncwriteatt(filename,'energy','method','https://github.com/smithmadisonm/DAS-surface-wave-processing/blob/main/processRBRduo.m')
            end
        end

   elseif ~strcmp(names{i},'ID')
        if strcmp(names(i),'time')
            ncwriteatt(filename,'time','units','days since 1970-01-01 00:00:00');
            ncwriteatt(filename,'time','long_name','Days since 1 January 1970');
            ncwriteatt(filename,'time','standard_name','time')
        end
        if strcmp(names(i),'lat')
            ncwriteatt(filename,'lat','units','degree_north')
            ncwriteatt(filename,'lat','long_name','latitude')
            ncwriteatt(filename,'lat','standard_name','latitude')
            ncwriteatt(filename,'lat','instrument','GPS')
        end
        if strcmp(names(i),'lon')
            ncwriteatt(filename,'lon','units','degree_east')
            ncwriteatt(filename,'lon','long_name','longitude')
            ncwriteatt(filename,'lon','standard_name','longitude')
            ncwriteatt(filename,'lon','instrument','GPS')
        end
        if strcmp(names(i),'watertemp')
            ncwriteatt(filename,'sea_water_temperature','units','degree_C')
            ncwriteatt(filename,'sea_water_temperature','long_name','sea water temperature at depth')
            ncwriteatt(filename,'sea_water_temperature','standard_name','sea_water_temperature')
            ncwriteatt(filename,'sea_water_temperature','instrument','RBR duo')
            ncwriteatt(filename,'sea_water_temperature','_FillValue',-999)
        end

        if strcmp(names(i),'sigwaveheight')
            ncwriteatt(filename,'sea_surface_wave_significant_height','units','m')
            ncwriteatt(filename,'sea_surface_wave_significant_height','long_name','significant wave height')
            ncwriteatt(filename,'sea_surface_wave_significant_height','standard_name','sea_surface_wave_significant_height')
            ncwriteatt(filename,'sea_surface_wave_significant_height','instrument','RBRduo')
        end
        if strcmp(names(i),'peakwaveperiod')
            ncwriteatt(filename,'sea_surface_wave_period_at_variance_spectral_density_maximum','units','s')
            ncwriteatt(filename,'sea_surface_wave_period_at_variance_spectral_density_maximum','long_name','peak of period orbital velocity spectra')
            ncwriteatt(filename,'sea_surface_wave_period_at_variance_spectral_density_maximum','standard_name','sea_surface_wave_period_at_variance_spectral_density_maximum')
            ncwriteatt(filename,'sea_surface_wave_period_at_variance_spectral_density_maximum','instrument','RBRduo')
        end
        if strcmp(names(i),'centroidwaveperiod')
            ncwriteatt(filename,'sea_surface_wave_mean_period','units','s')
            ncwriteatt(filename,'sea_surface_wave_mean_period','long_name','centroid (mean) period orbital velocity spectra')
            ncwriteatt(filename,'sea_surface_wave_mean_period','standard_name','sea_surface_wave_mean_period')
            ncwriteatt(filename,'sea_surface_wave_mean_period','instrument','RBRduo')
        end

        if strcmp(names(i),'peakwaveperiod')
            ncwriteatt(filename,'sea_surface_wave_period_at_variance_spectral_density_maximum','units','s')
            ncwriteatt(filename,'sea_surface_wave_period_at_variance_spectral_density_maximum','long_name','peak of period orbital velocity spectra')
            ncwriteatt(filename,'sea_surface_wave_period_at_variance_spectral_density_maximum','standard_name','sea_surface_wave_period_at_variance_spectral_density_maximum')
            ncwriteatt(filename,'sea_surface_wave_period_at_variance_spectral_density_maximum','instrument','RBR Duo')
        end
     
        if strcmp(names(i),'depth')
            ncwriteatt(filename,'depth','units','m')
            ncwriteatt(filename,'depth','long_name','measured water depth at mooring location')
            ncwriteatt(filename,'depth','instrument','RBR Duo pressure sensor')
        end

        if strcmp(names(i),'')
            ncwriteatt(filename,'','units','')
            ncwriteatt(filename,'','long_name','')
        end




    end
end





end



