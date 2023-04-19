% calculate analytic shoaling and fetch dependence along Oliktok cable  
% to show possible wave gradients in open water
%
% J. Thomson, Apr 2023
clear all

%% input parameters
U10 = 5; % wind speed
H0 = 1; % offshore (incident) wave height
Tp = 6; % peak period

%% channel info 
load('/Users/jthomson/Dropbox/Projects/ArcticCable/2022/DAS/depthattenuation.mat')
xoffset = -120; % force channel 100 to be at 80 m 
dx = 2; % spacing between channels
x = channel * dx + xoffset; 

%% fetch growth
Xfetch = 9.8 * x ./ U10^2;
Hfetch = 0.002 * Xfetch.^0.5;
H = Hfetch .* U10^2 ./ 9.8; 

%% shoaling
for ci=1:length(channel)
    if depth(ci) > 0
        k(ci) = wavenumber( 1/Tp, depth(ci) );
    else 
        k(ci) = NaN;
    end
end

c = 9.8*Tp./(2*3.14) .* tanh(k.*depth');
cg = 0.5 * ( 1 + 2*k.*depth' ./ sinh( 2*k.*depth') );
shoaling = (cg ./cg(end) ).^0.5;


%% plots
figure(1), clf

subplot(2,1,1)
area(x./1000,-depth), hold on
plot(x./1000,-depth,'linewidth',10,'color',[.5 .3 .3])
set(gca,'FontSize',16,'fontweight','demi')
xlabel('Distance along cable [km]')
ylabel('Depth')

subplot(2,1,2)
plot(x./1000,H0*ones(1,length(x)),'K--','linewidth',4), hold on
plot(x./1000,H0+flipud(H),'-','linewidth',4), hold on
plot(x./1000,H0*shoaling,'-','linewidth',4)
set(gca,'FontSize',16,'fontweight','demi')
xlabel('Distance along cable [km]')
ylabel('Wave Height [m]')
legend('incident','Fetch growth','shoaling')
axis([0 40 0 3])

print -dpng shoaling_and_fetch.png
