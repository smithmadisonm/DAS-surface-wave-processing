% matlab script to process DAS data and compare with SWIFT buoy data
%
% J. Thomson, 12/2022
%
clear all 

DASdatapath = '/Users/jthomson/Desktop/dasncdf/';

channel = [100:100:18400];

%% SWIFT data

load('/Users/jthomson/Dropbox/Projects/ArcticCable/2022/SWIFT18_postprocessed/SWIFT18_CODAS_Aug2022_reprocessedIMU_RCprefitler_cleaned.mat')
for si=1:length(SWIFT)
    E_SWIFT(si,:) = SWIFT(si).wavespectra.energy;
end
f = SWIFT(1).wavespectra.freq;
df = median(diff(f));
time_SWIFT = [SWIFT.time];
Hs_SWIFT = [SWIFT.sigwaveheight];

%% loop thru DAS data, in directories by channel (location) and files by time

for ci=1:length(channel) % position loop

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
        [thisE thisf ] = pwelch(data,[],[],[],2);
        E_DAS(fi,:) = interp1(thisf, thisE, f); % interpolate onto SWIFT frequencies
        clear this* data

        % compare with SWIFT
        [tdiff matchedindex] = min( abs( time_DAS(fi) - time_SWIFT ) );
        if tdiff < 1/24
            E_ratio(fi,:) = E_SWIFT(matchedindex,:) ./ E_DAS(fi,:);
        else
            E_ratio(fi,:) = NaN * E_SWIFT(matchedindex,:) ./ E_DAS(fi,:);
        end

        %figure(1),
        %loglog(f,E);

    end

    E_coef = nanmean(E_ratio); % time average of the ratio at each frequency
    E_coefstddev = nanstd(E_ratio); % standard deviation of the ratio at each frequency
    E_DAS = E_DAS.*E_coef; % calibrated sea surface elevation energy spectra from the DAS

    Hs_DAS = 4*nansum(E_DAS*df,2).^.5;

    save([DASdatapath '/' num2str(channel(ci))  '/DASspecta_channel'  num2str(channel(ci)) '.mat'],'E*','f*','time*')

    %% plots for this channel

    figure(2)

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
    axis([ax(1) ax(2) 0 0.5])
    datetick('x','keeplimits')
    legend('DAS')
    ylabel('f [Hz]')


    subplot(3,1,3)
    pcolor([SWIFT.time],f,log10(E_SWIFT'))
    shading flat
    datetick
    axis([ax(1) ax(2) 0 0.5])
    legend('SWIFT')
    ylabel('f [Hz]')

    print('-dpng',[DASdatapath  '/' num2str(channel(ci)) '/DAS-SWIFTcompare_channel'  num2str(channel(ci)) '.png'])


end

%% combined results (all chanels, all times)
% 
% save('/Users/jthomson/Dropbox/Projects/ArcticCable/2022/DASresults.mat','time','channel','Hs*')
% 
% %% more plots (all channels)
% 
% figure(3),clf
% pcolor(time,channel,Hs_DAS_all'./1000)
% datetick
% ylabel('channel (distance)')
% cb = colorbar;
% cb.Label.String = 'DAS wave height';
% 
% print('-dpng','/Users/jthomson/Dropbox/Projects/ArcticCable/2022/DAS_Hs_channel.png')
% 



