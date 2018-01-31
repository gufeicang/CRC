% This script draws figure 2 of paper [1] (see paper for more details)
% [1] A Novel Adaptive Approach to Mode Decomposition for Multicomponent
% Signals, Duong-Hung Pham, and Sylvain Meignen.
% Duong Hung PHAM 

% 2017, Jul 1st


clear all;
close all; %clc; 
set(0,'DefaultAxesFontSize',18);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Letters_IEEE_2017/Latex/figures';

%% Signal and parameters
N = 2^10; % signal length
Nx=N/2;
Nfft = 2*Nx; % number of FFT bins
Fs = N;

%R  = 1 ; % number of freq spacings
%df = linspace(.01,.01,R); 
a=1;

%Ks = 24 % resolution of submatrice 
SNR = 0;
cas = 6;
NumC = 1;

 if (cas == 1)
  s1 = 0.5*fmconst(N,0.15+rand*.25) ; 
  s2 = zeros(N,1);
  s2(N/2+round(30*rand)) = 15;
  numC = 2; 
 elseif (cas == 2)
  s1 = fmsin(N,0.05+rand*.05,.25);   
  s2 = fmlin(N,0.25+rand*.05,.45); 
  numC = 2; 
 elseif (cas == 3)
  s1 = 0.8*fmsin(N,.05+rand*0.05,.2+rand*0.15);
  s1(1:50)=0;
  s1(N/2-100:end)=0;
  s1(N/2+100:N-50)= 0.8*fmconst(N/2-99-50,.25);
  
  s2 = zeros(N,1);
  s2(101:N/2-100)= 0.8*fmlin(N/2-200,.42,.3);
  s2(N/2-round(20*rand)) = 15;
  numC = 2; 
  
 elseif (cas == 4)  
  s1 = fmlin(N,0.1+rand*.05,0.1+rand*.25);   
  s2 = fmlin(N,0.1+rand*.2,.4);
  numC = 2; 
 elseif (cas == 5)
  t = (0:N-1)/N;
  s = exp(2*pi*1i*(180*t+ 100*t.^2));
  s = s(:);
  numC = 1; 
 elseif (cas == 6)
  s = fmsin(N,.05+rand*0.1,.3+rand*0.15);
  numC = 1; 
 elseif (cas == 7)
  s = 0.5*fmconst(N,0.25) ;
  numC = 1; 
 else
  s = zeros(N,1);
  s(N/2-round(rand*75)) = 1;
  numC = 1; 
 end    
    
load 'noise_fig23.mat';
%noise = randn(N,1);
 
 s=real(s);
 sig=real(s);
 s = sigmerge(sig,noise,SNR) ; % noisy signal
 s = s(:);
 disp([num2str(round(20*log10(norm(sig)/norm(s-sig)))) ' dB for input SNR']);

prec = 10^(-6);
%[H L] = roundgauss(Nx, prec); 
L  = sqrt(Nfft/a);
l  = floor(sqrt(-Nfft*log(prec)/(a*pi)))+1;
w  = amgauss(2*l+1,l+1,L);
H  = w/norm(w);

%% Ridges and Basins
[q,reassign,dx] = get_tfrs(s,Nfft,H,L); % get TF representations
[~,~,dx_noround] = get_tfrs_noround(s,Nfft,H,L); % get TF representations

t = (0:N-1)/Fs;
f = (-0:Nfft/2-1)*Fs/N;

%% Compute the basins of attraction associated to the different ridges/modes
%[contours, rzeros, ang, dx_vec,basins, ang1] = get_contour_basins(q,dx,reassign,Ks);

%% Rounding effect
ang1 = angle(dx);
ang1 = mod(ang1,pi);

FigHandle(1) = figure();
imagesc(ang1);
set(gca,'YDir','normal'); xlabel('time');ylabel('frequency');
set(gca,'xtick',[],'ytick',[]);
%%%%%%%%%%%%%%%%%%%%%% print Figures
 export_fig(FigHandle(1), ... % figure handle
     sprintf('%s/image_with_rounding', chemin0),... % name of output file without extension
     '-painters', ...      % renderer
     '-transparent', ...   % renderer
     '-pdf', ...           % file format
     '-r500' );             % resolution in dpi

ang1(abs(real(dx_noround))<0.025|abs(imag(dx_noround))<0.025)=[];
FigHandle(1) = figure();
histogram(ang1,90); 

 %%%%%%%%%%%%%%%%%%%%%% print Figures
 export_fig(FigHandle(1), ... % figure handle
     sprintf('%s/histogram_with_rounding', chemin0),... % name of output file without extension
     '-painters', ...      % renderer
     '-transparent', ...   % renderer
     '-pdf', ...           % file format
     '-r500' );             % resolution in dpi

 