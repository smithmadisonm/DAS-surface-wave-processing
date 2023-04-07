% matlab script to process DAS data and compare with SWIFT buoy data
%
% J. Thomson, 12/2022
%
clear all, close all

tic

DASdatapath = '/Users/jthomson/Desktop/ArcticCable_local/dasncdf/';
%DASdatapath = '/Volumes/Data/BeaufortChukchi/OliktokDAS/DAS_Aug2022';

channel = 7960%[100:20:18400]; %7960
plotflag = true;
movieflag = false;
minHs = 0.05; % min wave height
fs_DAS = 2; 
fs_SWIFT = 25;

%% depth attenuation
%load('/Volumes/Data/BeaufortChukchi/OliktokDAS/depthattenuation.mat','f','depth','channel','attenuation');
%load('/Volumes/Data/BeaufortChukchi/OliktokDAS/depthattenuation.mat','attenuation');  % freq x channel (42 x 916 or 85 x 916 for extended freq)
load('~/Dropbox/Projects/ArcticCable/2022/DAS/depthattenuation.mat','attenuation');  % freq x channel (42 x 916 or 85 x 916 for extended freq)


%% SWIFT data, with option to reprocess for higher freq (max 1 Hz instead of 0.5 Hz)

%load('/Users/jthomson/Dropbox/Projects/ArcticCable/2022/SWIFT18_postprocessed/SWIFT18_CODAS_Aug2022_reprocessedIMU_RCprefitler_cleaned.mat')
%load('/Volumes/Data/BeaufortChukchi/OliktokDAS/SWIFT18_OliktokPt_Aug2022/SWIFT18_OliktokPt_Aug2022_reprocessedIMU_RCprefilter_displacements.mat')
load('~/Dropbox/Projects/ArcticCable/2022/SWIFT18_postprocessed/SWIFT18_OliktokPt_Aug2022_reprocessedIMU_RCprefilter_displacements.mat')

f = SWIFT(1).wavespectra.freq;
df = median(diff(f));
newf = [f(1):df:1]; 
f = newf;

for si=1:length(SWIFT)
    if length(isfinite(SWIFT(si).z)) > 10000
        %E_SWIFT(si,:) = SWIFT(si).wavespectra.energy;
        [newE f ] = pwelch(SWIFT(si).z(isfinite(SWIFT(si).z) ),[],[],f,fs_SWIFT);
        E_SWIFT(si,:) = newE;  
        Hs_SWIFT(si) = 4 * nansum(newE*df).^0.5; SWIFT(si).sigwaveheight = 4 * nansum(newE*df).^0.5;
        var_SWIFT(si) = var( SWIFT(si).z(isfinite(SWIFT(si).z) ) );
        SWIFT(si).wavespectra.energy = newE'; 
        SWIFT(si).wavespectra.freq = f';
    else
      E_SWIFT(si,:) = NaN(1,length(f)); 
      Hs_SWIFT(si) = NaN;
      var_SWIFT(si) = NaN;
      SWIFT(si).wavespectra.energy = NaN(1,length(f))'; 
      SWIFT(si).wavespectra.freq = f';
    end
      SWIFT(si).wavespectra.a1 = NaN(1,length(f))'; 
      SWIFT(si).wavespectra.b1 = NaN(1,length(f))'; 
      SWIFT(si).wavespectra.a2 = NaN(1,length(f))'; 
      SWIFT(si).wavespectra.b2 = NaN(1,length(f))'; 
      SWIFT(si).wavespectra.check = NaN(1,length(f))';     
end

SWIFT = rmfield(SWIFT,'x');
SWIFT = rmfield(SWIFT,'y');
SWIFT = rmfield(SWIFT,'z');
SWIFT = rmfield(SWIFT,'u');
SWIFT = rmfield(SWIFT,'v');
SWIFT = rmfield(SWIFT,'rawlon');
SWIFT = rmfield(SWIFT,'rawlat');

time_SWIFT = [SWIFT.time];
%Hs_SWIFT = [SWIFT.sigwaveheight];

%% loop thru DAS data, in directories by channel (location) and files by time

for ci=1:length(channel) % position loop

    vidObj = VideoWriter([DASdatapath '/' num2str(channel(ci))  '/DASchannel'  num2str(channel(ci)) '.mpeg' ],'MPEG-4');
    open(vidObj);

    flist = dir([DASdatapath '/' num2str(channel(ci)) '/*.ncdf']);

    for fi = 1:length(flist) % time loop

        data = ncread([DASdatapath  '/' num2str(channel(ci)) '/' flist(fi).name],'data');

        year = str2num( flist(fi).name(11:14) );
        month = str2num( flist(fi).name(15:16) );
        day = str2num( flist(fi).name(17:18) );
        hour = str2num( flist(fi).name(20:21) );
        minute = str2num( flist(fi).name(22:23) );
        second = str2num( flist(fi).name(24:25) );
        time_DAS(fi) = datenum(year, month, day, hour, minute, second);
        clear year month day hour minute second

        % spectra
        cleandata = filloutliers( detrend(data), 'linear');
        [thisE thisf ] = pwelch(data,[],[],[],fs_DAS);
        E_DAS_raw(fi,:) = interp1(thisf, thisE, f); % interpolate onto SWIFT frequencies
        E_DAS_depthcorrected(fi,:) = E_DAS_raw(fi,:) .* attenuation(:,fi)';
        var_DAS(fi) = var(data);

        % compare with SWIFT
        [tdiff matchedindex] = min( abs( time_DAS(fi) - time_SWIFT ) );
        if tdiff < 1/24 & Hs_SWIFT(matchedindex) > minHs
            E_ratio(fi,:) = E_SWIFT(matchedindex,:) ./ E_DAS_depthcorrected(fi,:);
            var_ratio(fi) = var_SWIFT(matchedindex) ./ var_DAS(fi);
            Hs_SWIFT_atDAStime(fi) = Hs_SWIFT(matchedindex);
        else
            E_ratio(fi,:) = NaN * E_SWIFT(matchedindex,:) ./ E_DAS_depthcorrected(fi,:);
            var_ratio(fi) = NaN;
            Hs_SWIFT_atDAStime = NaN;
        end

        if movieflag
            figure(1), clf       
            subplot(2,1,1)
            plot(data,'b'), hold on, plot(cleandata,'r')
            title(flist(fi).name)
            subplot(2,1,2)
            loglog(f,E_SWIFT(matchedindex,:),'k:',thisf,thisE,'b', f,E_DAS_raw(fi,:),'r',f,E_DAS_depthcorrected(fi,:),'g');
            axis square
            currFrame = getframe(gcf);
            writeVideo(vidObj,currFrame);
        end
        
        clear this* data
        
    end

    close(vidObj);
    
    E_coef = nanmedian(E_ratio); % time average of the ratio at each frequency
    var_coef = nanmedian(var_ratio); % time average of the variance ratio
    E_coefstddev = nanstd(E_ratio); % standard deviation of the ratio at each frequency
    var_coefstdev = nanstd(var_ratio);
    E_DAS = E_DAS_depthcorrected.*E_coef; % calibrated sea surface elevation energy spectra from the DAS
    var_DAS = var_DAS .* var_coef; % calibrated variance

    Hs_DAS = 4*nansum(E_DAS.*df,2).^.5;
    Hs_DAS_var = 4*(var_DAS).^.5;

    save([DASdatapath '/' num2str(channel(ci))  '/DASspecta_channel'  num2str(channel(ci)) '.mat'],'E*','f*','time*','Hs*')

    %% plots for this channel

    if plotflag
        
    figure(2), clf

    subplot(3,1,1)
    plot(time_DAS,Hs_DAS,'x',time_SWIFT,Hs_SWIFT)
    datetick
    title(['DAS channel ' num2str(channel(ci))])
    legend('DAS','SWIFT')
    ax = axis;
    ylabel('Hs [m]')

    subplot(3,1,2)
    pcolor(time_DAS,f,log10(E_DAS'))
    shading flat
    axis([ax(1) ax(2) 0 1])
    datetick('x','keeplimits')
    legend('DAS')
    ylabel('f [Hz]')


    subplot(3,1,3)
    pcolor([SWIFT.time],f,log10(E_SWIFT'))
    shading flat
    datetick
    axis([ax(1) ax(2) 0 1])
    legend('SWIFT')
    ylabel('f [Hz]')

    print('-dpng',[DASdatapath  '/' num2str(channel(ci)) '/DAS-SWIFTcompare_channel'  num2str(channel(ci)) '.png'])

    figure(3), clf
    loglog(f,E_ratio,'color',[.5 .5 .5]), hold on
    loglog(f,nanmean(E_ratio),'k','linewidth',3)
    loglog(f,nanmedian(E_ratio),'c','linewidth',3)
    xlabel('f [Hz]'), ylabel('E ratio')
    title(['DAS channel ' num2str(channel(ci))])
    print('-dpng',[DASdatapath  '/' num2str(channel(ci)) '/Eratio_channel'  num2str(channel(ci)) '.png'])
  
    figure(4), clf
    plot(Hs_SWIFT_atDAStime, Hs_DAS, 'kx', Hs_SWIFT_atDAStime, Hs_DAS_var, 'rx')
    xlabel('SWIFT H_s [m]'), ylabel('DAS H_s [m]')
    legend(['spectral Hs, R^2 = ' num2str(corr(Hs_SWIFT_atDAStime(isfinite(Hs_DAS))', Hs_DAS(isfinite(Hs_DAS))).^2)]... 
        ,['variance Hs, R^2 = ' num2str(corr(Hs_SWIFT_atDAStime(isfinite(Hs_DAS_var))', Hs_DAS_var(isfinite(Hs_DAS_var))').^2)])
    print('-dpng',[DASdatapath  '/' num2str(channel(ci)) '/DAS-SWIFTcompare_channel'  num2str(channel(ci)) '_Hs.png'])

    else
    end

end

toc


%% combined results (all chanels, all times)
% 

for ci=1:length(channel) 
    
    load([DASdatapath '/' num2str(channel(ci))  '/DASspecta_channel'  num2str(channel(ci)) '.mat'])
    
    figure(5), 
    cmap = colormap;
    cindex = ceil(ci./length(channel)*length(cmap));
    plot(Hs_SWIFT_atDAStime,Hs_DAS,'.','color',cmap(cindex,:))
    hold on
    
    figure(6),
    cmap = colormap;
    cindex = ceil(ci./length(channel)*length(cmap));
    loglog(f,E_coef,'color',cmap(cindex,:) )
    hold on
    
end

figure(5)
xlabel('SWIFT H_s [m]'), ylabel('DAS H_s [m]')
plot([0 1],[0 1],'k:')
print('-dpng',[ DASdatapath '/Hscatter_allchannels.png' ])

figure(6)
xlabel('f [Hz]'), ylabel('E coef')
print('-dpng',[ DASdatapath '/Ecoef_allchannels.png' ])



