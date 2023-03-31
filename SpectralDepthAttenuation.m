% spectral energy correction for depth attenuation of surface waves along
% Oliktok Point DAS cable run

%% cable info
cableinfo = importdata('/Volumes/Data/BeaufortChukchi/OliktokDAS/CODAS_info.csv');
depth = cableinfo.data(:,6);
channel = cableinfo.data(:,7);

plot(channel, depth)
xlabel('DAS channel'), ylabel('Depth [m]')
print('-dpng','/Volumes/Data/BeaufortChukchi/OliktokDAS/CODAS_depthprofile.png');

%% frequencies to use (match SWIFT)
load('/Users/jthomson/Dropbox/Projects/ArcticCable/2022/SWIFT18_postprocessed/SWIFT18_CODAS_Aug2022_reprocessedIMU_RCprefitler_cleaned.mat')
f = SWIFT(1).wavespectra.freq;

%%

attenuation = NaN(length(f), length(channel));

tic
for ci=1:length(channel)
    
    ci
    
    if depth(ci) > 0
        for fi=1:length(f),
            k(fi) = wavenumber( f(fi), depth(ci) );
        end
        
        attenuation(:, ci) = cosh( k .* depth(ci) ) ;

    end
    
end

toc

attenuation = attenuation.^2; % square for energy

save('/Volumes/Data/BeaufortChukchi/OliktokDAS/depthattenuation.mat','f','depth','channel','attenuation');


%% plot
figure(1), clf
cmap = colormap;

for ci=1:100:length(channel) 
    ci
    cindex = ceil(ci./length(channel)*length(cmap));
    loglog(f,attenuation,'color',cmap(cindex,:) )
    hold on
    
end

print('-dpng','/Volumes/Data/BeaufortChukchi/OliktokDAS/depthattenuation_spectral.png');


