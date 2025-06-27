% process RBR Duo pressure records using SWIFT wave codes
% assumes continous recording of pressure and temp at 1 Hz
% J. Thomson, Nov 2019

clear all, close all

%filename = 'CODA_S2P1_BottomPT_Nov2019-Sep2020'; lat = 70.619650; lon = -149.575383;
%filename = 'CODA_S2-P3_BottomPT_Aug-Nov2019'; lat = 70.718450; lon = -149.508933;
%filename = 'CODA_S3P1_BottomPT_Nov2019-Sep2020'; lat = 70.252550; lon = -145.994867;
%filename = 'CODA_S1P1_BottomPT_Nov2019-Sep2020'; lat = 70.34613; lon = -162.05704;
%filename = 'CODA_S1-P2_BottomPT_Nov2019';  lat = 70.39473;  lon = -162.13675;
%filename = 'CODA_S1-P3_BottomPT_Nov2019'; lat = 70.43995; lon = -162.20300;
%filename = 'CODA_EluitkakPass_BottomPT_Aug2019-Aug2020'; lat = 71.36; lon = -156.35;
%filename = '203193_20220903_0202_T1C'; lat = NaN; lon = NaN;
%filename = 'ArcticCable_Site1_Apr-Sep2023'; lat = 70.5594; lon = -150.0071;
%filename = 'ArcticCable_Site2_Apr-Sep2023'; lat = 70.6546; lon = -150.0001;
%filename = 'ArcticCable_Site3_Apr-Sep2023'; lat = 70.7394; lon = -150.003;
%filename = 'MBanchor_RBR_16-25Jun2025'; lat = 46.690149; lon = -124.003540; 
%filename = 'MCanchor_RBR_17-26Jun2025'; lat = 46.717835; lon = -124.084016; 
%filename = 'OSanchor_RBR_17-26Jun2025'; lat = 46.72183; lon = -124.15081; 
%filename = 'SBanchor_RBR_16-27Jun2025'; lat = 46.658626; lon = -123.992052; 
%filename = 'MBanchor_RBR_16-27Jun2025'; lat = 46.690149; lon = -124.003540; 
filename = 'EBanchor_RBR_16-27Jun2025'; lat = 46.693024; lon = -123.959714; 



t1 = datenum(2025,6,16); % start time
t2 = datenum(2025,6,26); % end time

plotspectra = false;

%% fixed params

burstsecs = 1800; % burst length in seconds
atm = 10.13; % atmospheric pressure offset (dB)
mindepth = 1; % detect out-of-water times
minwaveheight = 0.05; % minimum wave height observable (signal to noise issue)
maxwaveheight = 6;
maxwaveperiod = 16;


%% read raw data, using a start and end time to split up big files
fid = RSKopen( [ filename '.rsk' ] );
RSK = RSKreaddata(fid,'t1',t1,'t2',t2);
fs = 1 / ( (RSK.data.tstamp(2)-RSK.data.tstamp(1))*24*3600 );
burstlength = floor( fs * burstsecs ); 

%% initialize movie
if plotspectra,
    % initialize video
    vidObj = VideoWriter([filename '_allspectra.avi'],'MPEG-4');
    open(vidObj);
end

%% parse into bursts
% would be better to adjust start to top of first hour
startindex = 1;
totallength = length(RSK.data.tstamp);

for bi = 1 : floor((totallength-startindex) ./ burstlength);
    
    RBR(bi).time = RSK.data.tstamp(startindex);
    RBR(bi).lat = lat; 
    RBR(bi).lon = lon;
    
    watertemp = RSK.data.values( [startindex:(startindex+burstlength)] ,1);
    RBR(bi).watertemp = mean(watertemp);
    
    pres = RSK.data.values( [startindex:(startindex+burstlength)] ,2) - atm;
    mooringdepth = nanmean(pres);
    waterdepth = mooringdepth; 

    % wave spectra and stats
    [ Hs, Tp, Hig, Tig, E, f ] = Pwaves(pres,fs); 
    
    RBR(bi).sigwaveheight = Hs;
    RBR(bi).peakwaveperiod = Tp;
    RBR(bi).IGwaveheight = Hig;
    RBR(bi).IGperiod = Tig;
    RBR(bi).wavespectra.energy = E';
    RBR(bi).wavespectra.freq = f';
    RBR(bi).depth = mooringdepth;
    
    
    % correct for depth attenuation
    RBRraw(bi) = RBR(bi);
    clear Ecorr*
    for j=1:length(f),
        k(j) = wavenumber( f(j), waterdepth );
    end
    attenuation = cosh( k .* waterdepth ) ./ cosh( k.*(waterdepth-mooringdepth) ) ;
    %attenuation = exp(k*depth); % should be same as cosh for deep water
    attenuation = attenuation.^2; % square for energy
    noise = attenuation > 100 | isnan(attenuation); % limit the size of the attenuation correction
    Ecorr = E.*attenuation;
    % extrapolate the equilibruim range as another product
    if find(noise,1) > 2,
        Ecorrextrap = Ecorr;
        Ecorrextrap( noise ) = ( Ecorr( find(noise,1) - 1 ) .* f( find(noise,1) - 1 ).^4 ) .* f(noise).^-4;
    else
        Ecorrextrap = NaN(size(Ecorr));
    end
    Ecorr( noise ) = NaN; % cut it off when correction too big, don't amplify noise
    Hscorr = 4 * sqrt( nansum(Ecorr) * median(diff(f)) ); 
    Hscorrextrap = 4 * sqrt( nansum(Ecorrextrap) * median(diff(f)) ); 
    [Emax fpindexcorr ] = max(Ecorr);
    [Emax fpindexcorrextrap ] = max(Ecorrextrap);
    Tpcorr = f(fpindexcorr).^-1;
    Tpcorrextrap = f(fpindexcorrextrap).^-1;
    
    RBRcorr(bi) = RBR(bi); % initilaize 
    RBRcorr(bi).sigwaveheight = Hscorr;
    RBRcorr(bi).peakwaveperiod = Tpcorr;
    RBRcorr(bi).wavespectra.energy = Ecorr';
    
    RBRcorrextrap(bi) = RBR(bi); % initilaize 
    RBRcorrextrap(bi).sigwaveheight = Hscorrextrap;
    RBRcorrextrap(bi).peakwaveperiod = Tpcorrextrap; % Tp must be the same, because extrap decreases from peak
    RBRcorrextrap(bi).wavespectra.energy = Ecorrextrap';
    
    startindex = startindex + burstlength;
    
    if plotspectra
        figure(1), clf
        subplot(2,1,1)
        plot(pres)
        ylabel('pressure [db]'), xlabel('index []')
        title( [filename '     ' datestr(RBR(bi).time) ],'interp','none')
        subplot(2,1,2)
        loglog(f,E, f,Ecorr, f, Ecorrextrap, '--' ,'linewidth',2)
        axis([1e-3 1 1e-6 1e1])
        ylabel('Energy density [m^2/Hz]')
        xlabel('Frequency [Hz]')
        legend('Raw (bottom)','Corrected to surface','Corrected and extrapolated','Location','NorthEastOutside')
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
    end
    
end

if plotspectra
    close(vidObj);
end

%% QC waves 

for bi=1:length(RBR),

    if RBR(bi).sigwaveheight > maxwaveheight | RBRcorr(bi).sigwaveheight > maxwaveheight | RBRcorrextrap(bi).sigwaveheight > maxwaveheight | ... 
             RBR(bi).peakwaveperiod > maxwaveperiod | RBRcorr(bi).peakwaveperiod > maxwaveperiod | RBRcorrextrap(bi).peakwaveperiod > maxwaveperiod
        RBR(bi).sigwaveheight = NaN;
        RBRcorr(bi).sigwaveheight = NaN;
        RBRcorrextrap(bi).sigwaveheight = NaN;
        RBR(bi).peakwaveperiod = NaN;
        RBRcorr(bi).peakwaveperiod = NaN;
        RBRcorrextrap(bi).peakwaveperiod = NaN;
        RBR(bi).IGwaveheight = NaN;
        RBR(bi).IGperiod = NaN;
    end

    if RBR(bi).sigwaveheight < minwaveheight | RBRcorr(bi).sigwaveheight < minwaveheight | RBRcorrextrap(bi).sigwaveheight < minwaveheight
        %RBR(bi).sigwaveheight = NaN;
        %RBRcorr(bi).sigwaveheight = NaN;
        %RBRcorrextrap(bi).sigwaveheight = NaN;
        RBR(bi).peakwaveperiod = NaN;
        RBRcorr(bi).peakwaveperiod = NaN;
        RBRcorrextrap(bi).peakwaveperiod = NaN;
    end
    
    RBR(bi).peakwavedirT = NaN; % remove all wave directions (no vel data with which to make dir est)
    
end

%% QC for in water 


RBR(find([RBR.depth]<mindepth)) = [];
RBRcorr(find([RBRcorr.depth]<mindepth)) = [];
RBRcorrextrap(find([RBRcorrextrap.depth]<mindepth)) = [];


%% plot and save

figure(2), 

ax(1) = subplot(3,1,1);
plot([RBR.time],[RBR.sigwaveheight],[RBRcorr.time],[RBRcorr.sigwaveheight],...
    [RBRcorrextrap.time],[RBRcorrextrap.sigwaveheight],'--')
datetick
ylabel('H_s [m]')
legend('uncorrected','corrected','correct and extrapolated','Location','Northwest')
title([filename],'interpreter','none')


ax(2) = subplot(3,1,2);
plot([RBR.time],[RBR.peakwaveperiod],[RBRcorr.time],[RBRcorr.peakwaveperiod],...
    [RBRcorrextrap.time],[RBRcorrextrap.peakwaveperiod],'--')
datetick
ylabel('T_p [m]')
legend('uncorrected','corrected','correct and extrapolated','Location','Northwest')
linkaxes(ax,'x')
datetick

ax(3) = subplot(3,1,3);
plot([RBR.time],[RBR.IGwaveheight],'k.')
datetick
ylabel('H_{IG} [m]')

print('-dpng',[filename  '_' datestr(t1,1) '_' datestr(t2,1) '_wavestats.png'])


save( [filename '_' datestr(t1,1) '_' datestr(t2,1) '_uncorrected' ], 'RBR')

save( [filename '_' datestr(t1,1) '_' datestr(t2,1) '_depthcorrected' ], 'RBRcorr')

save( [filename '_' datestr(t1,1) '_' datestr(t2,1) '_depthcorrectedextrapolated' ], 'RBRcorrextrap')

%% spectrogram

figure(3), clf

for bi=1:length(RBR)
    E(bi,:) = RBR(bi).wavespectra.energy;
end

pcolor([RBR.time],f,log10(E'))
shading flat
datetick
ylabel('f [Hz]')
set(gca,'yscale','log')
title([filename],'interpreter','none')

print('-dpng',[filename  '_' datestr(t1,1) '_' datestr(t2,1) '_spectrogram.png'])





