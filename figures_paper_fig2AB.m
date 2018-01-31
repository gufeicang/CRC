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
a=1;

Ks = 16; % resolution of submatrice 
SNR = 0;
cas = 8;
NumC = 1;

load 'noise_fig23.mat';
noise = randn(N,1);
s=noise;
s = s(:);
  
prec = 10^(-6);
%[H L] = roundgauss(Nx, prec); 
L  = sqrt(Nfft/a);
l  = floor(sqrt(-Nfft*log(prec)/(a*pi)))+1;
w  = amgauss(2*l+1,l+1,L);
H  = w/norm(w);

%% Ridges and Basins
[q_noround,reassign_noround,dx_noround] = get_tfrs_noround(s,Nfft,H,L); % get TF representations

[q,reassign,dx] = get_tfrs(s,Nfft,H,L); % get TF representations

t = (0:N-1)/Fs;
f = (-0:Nfft/2-1)*Fs/N;

%% Rounding effect

%dx_noround1 = dx_noround.*(1-((abs(real(dx_noround))<0.5).*(abs(imag(dx_noround))<0.5)));

 ang1 = angle(dx_noround);
% figure();
% histogram(ang1,180);

FigHandle(1) = figure();
ang1(abs(real(dx_noround))<0.025|abs(imag(dx_noround))<0.01)=[];
ang1 = mod(ang1,pi);
histogram(ang1,90); 

 
ang2 = angle(dx);
FigHandle(2) = figure();
ang2(abs(real(dx_noround))<0.025|abs(imag(dx_noround))<0.025)=[];
ang2 = mod(ang2,pi);
histogram(ang2,90); 

 %%%%%%%%%%%%%%%%%%%%%% print Figures
 export_fig(FigHandle(1), ... % figure handle
     sprintf('%s/dx_arg_without_rounding', chemin0),... 
     '-painters', ...      % renderer
     '-transparent', ...   % renderer
     '-pdf', ...           % file format
     '-r500' );             % resolution in dpi
 
  %%%%%%%%%%%%%%%%%%%%%% print Figures
 export_fig(FigHandle(2), ... % figure handle
     sprintf('%s/dx_arg_with_rounding', chemin0),... 
     '-painters', ...      % renderer
     '-transparent', ...   % renderer
     '-pdf', ...           % file format
     '-r500' );             % resolution in dpi
 