% matlab script to process DAS data and compare with SWIFT buoy data
% 
% J. Thomson, 12/2022
%

DASdatapath = '/Users/jthomson/Desktop/dasncdf/';

channel = [100:20:1840];

%% loop thru DAS data

for ci=1:length(channel)

flist = dir([DASdatapath '/' num2str(channel(ci)) '/*.ncdf']);

for fi = 1:length(flist)

    data = ncread([DASdatapath  '/' num2str(channel(ci)) '/' flist(fi).name],'data');
    
    year = str2num( flist(fi).name(11:14) );
    month = str2num( flist(fi).name(15:16) );
    day = str2num( flist(fi).name(17:18) );
    hour = str2num( flist(fi).name(20:21) );
    minute = str2num( flist(fi).name(22:23) );
    second = str2num( flist(fi).name(24:25) );
    time(fi) = datenum(year, month, day, hour, minute, second);

    % spectra
    [E f ] = pwelch(data,[],[],[],2);
    
    figure(1), 
    loglog(f,E);

    allE_DAS(fi,:) = E;
    df = median(diff(f));
    Hs_DAS(fi) = 4*sum(E*df).^.5;
    Hs_DAS_all(fi,ci) = 4*sum(E*df).^.5;


end

f_DAS = f;
save([DASdatapath '/' num2str(channel(ci))  '/DASspecta.mat'],'*DAS')

%% SWIFT data

load('/Users/jthomson/Dropbox/Projects/ArcticCable/2022/SWIFT18_postprocessed/SWIFT18_CODAS_Aug2022_reprocessedIMU_RCprefitler_cleaned.mat')
for si=1:length(SWIFT)
    allE_SWIFT(si,:) = SWIFT(si).wavespectra.energy;
end
f_SWIFT = SWIFT(1).wavespectra.freq;

%% channel plots

figure(2)

subplot(3,1,1)
plot(time,Hs_DAS/1000,'x',[SWIFT.time],[SWIFT.sigwaveheight])
datetick
title(['DAS channel ' num2str(channel(ci))])
legend('DAS','SWIFT')
ax = axis;
ylabel('Hs [m]')

subplot(3,1,2)
pcolor(time,f_DAS,log10(allE_DAS'))
shading flat
axis([ax(1) ax(2) 0 0.5])
datetick('x','keeplimits')
legend('DAS')
ylabel('f [Hz]')


subplot(3,1,3)
pcolor([SWIFT.time],f_SWIFT,log10(allE_SWIFT'))
shading flat
datetick
axis([ax(1) ax(2) 0 0.5])
legend('SWIFT')
ylabel('f [Hz]')

print('-dpng',[DASdatapath  '/' num2str(channel(ci)) '/DAS-SWIFTcompare.png'])


end

save('/Users/jthomson/Dropbox/Projects/ArcticCable/2022/DASresults.mat','time','channel','Hs_DAS_all')

%% more plots (all channels)

figure(3),clf
pcolor(time,channel,Hs_DAS_all'./1000)
datetick
ylabel('channel (distance)')
cb = colorbar;
cb.Label.String = 'DAS wave height';

print('-dpng','/Users/jthomson/Dropbox/Projects/ArcticCable/2022/DAS_Hs_channel.png')

