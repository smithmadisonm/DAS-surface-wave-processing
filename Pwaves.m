function [ Hs, Tp, Hig, Tig, E, f ] = Pwaves(p,fs) 

% matlab function to read and process pressure measurements
%   to estimate significant wave height, peak period,
%   infragravity wave height and periodn, and spectra
%
% Inputs are pressure time series (in decibar or m) and sampling rate [Hz]
%
% Sampling rate must be at least 1 Hz and the same
% for all variables.  Additionaly,  time series data must 
% have at least 1024 points and all be the same size.
%
% Outputs will be '9999' for invalid results.
%
% Outputs can be supressed, in order, full usage is as follows:
%
%   [ Hs, Tp, Hig, Tig, E, f ] = Pwaves(p,fs) 

% J. Thomson, Oct 2020, adapted from PUVspectra.m (2002)
%
%#codegen
  

%% fixed parameters
wsecs = 900;   % window length in seconds, should make 2^N samples
merge = 1;      % freq bands to merge, must be odd?
maxf = 1;       % frequency cutoff for output
   

%% begin processing, if data sufficient
pts = length(p);       % record length in data points

if pts >= 2*wsecs & fs>=1,  % minimum length and quality for processing


%% break into windows (use 75 percent overlap)
w = round(fs * wsecs);  % window length in data points
if rem(w,2)~=0, w = w-1; else end  % make w an even number
windows = floor( 4*(pts/w - 1)+1 );   % number of windows, the 4 comes from a 75% overlap
dof = 2*windows*merge; % degrees of freedom
% loop to create a matrix of time series, where COLUMN = WINDOW 
pwindow = zeros(w,windows);
for q=1:windows, 
  	pwindow(:,q) = p(  (q-1)*(.25*w)+1  :  (q-1)*(.25*w)+w  );  
end

%% detrend individual windows (full series already detrended)
for q=1:windows
pwindow(:,q) = detrend(pwindow(:,q));
end

%% taper and rescale (to preserve variance)
% form taper matrix (columns of taper coef)
taper = sin ( (1:w) * pi/w )' * ones(1,windows); 
% taper each window
pwindowtaper = pwindow .* taper;
% now find the correction factor (comparing old/new variance)
factp = sqrt( var(pwindow) ./ var(pwindowtaper) );
% and correct for the change in variance
% (mult each window by it's variance ratio factor)
pwindowready = (ones(w,1)*factp).* pwindowtaper;


%% FFT
% calculate Fourier coefs
Pwindow = fft(pwindowready);
% second half of fft is redundant, so throw it out
Pwindow( (w/2+1):w, : ) = [];
% throw out the mean (first coef) and add a zero (to make it the right length)  
Pwindow(1,:)=[]; 
Pwindow(w/2,:)=0; 
% POWER SPECTRA (auto-spectra)
PPwindow = real ( Pwindow .* conj(Pwindow) );


%% merge neighboring freq bands (number of bands to merge is a fixed parameter)
% initialize
PPwindowmerged = zeros(floor(w/(2*merge)),windows);

for mi = merge:merge:(w/2) 
   	PPwindowmerged(mi/merge,:) = mean( PPwindow((mi-merge+1):mi , : ) );
end
% freq range and bandwidth
n = (w/2) / merge;                         % number of f bands
Nyquist = .5 * fs;                % highest spectral frequency 
bandwidth = Nyquist/n ;                    % freq (Hz) bandwitdh
% find middle of each freq band, ONLY WORKS WHEN MERGING ODD NUMBER OF BANDS!
f = 1/(wsecs) + bandwidth/2 + bandwidth.*(0:(n-1)) ; 


%% ensemble average windows together
% take the average of all windows at each freq-band
% and divide by N*samplerate to get power spectral density
% the two is b/c Matlab's fft output is the symmetric FFT, and we did not use the redundant half (so need to multiply the psd by 2)
E = mean( PPwindowmerged.' ) / (w/2 * fs  );


%% bulk parameters

fwaves = f>0.05 & f<maxf; % frequency cutoff for wave stats
fig = f>0.005 & f<0.05; % infragravity


% significant wave height
Hs  = 4*sqrt( sum( E(fwaves) ) * bandwidth);
Hig  = 4*sqrt( sum( E(fig) ) * bandwidth);


%  energy period 
fe = sum( f(fwaves).*E(fwaves) )./sum( E(fwaves) );
feig = sum( f(fig).*E(fig) )./sum( E(fig) );
Te = 1./fe;
Tig = 1./feig;
Tp = Te;

% peak period
%[~ , fpindex] = max(E);
%Tp = 1./f(fpindex);


%% prune high frequency results
E( f > maxf ) = [];
f( f > maxf ) = [];


else % if not enough points or sufficent sampling rate or data, give 9999
  
     Hs = 9999;
     Tp = 9999; 
     Hig = 9999;
     Tig = 9999; 
     E = 9999; 
     f = 9999;

end



