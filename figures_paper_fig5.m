% This script draws figure 5 of paper [1] (see paper for more details)
% [1] A Novel Adaptive Approach to Mode Decomposition for Multicomponent
% Signals, Duong-Hung Pham, and Sylvain Meignen.
% Duong Hung PHAM 

% 2017, Jul 1st

clear all;
close all; clc; 
set(0,'DefaultAxesFontSize',18);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/Letters_IEEE_2017/Latex/figures';

%% Signal and parameters

a=1;

Ks = 30
NumC = 3; % Number of modes

load -ascii batsig.txt
s0 = [batsig];
N = length(s0);

Nx=N/2;
Nfft = 2*Nx; % number of FFT bins
Fs = N;

noise = randn(N,1);
SNR = 5;

s = sigmerge(s0,noise,SNR);

disp([num2str(round(20*log10(norm(s0)/norm(s-s0)))) ' dB for input SNR']);

prec = 10^(-6);
%[H L] = roundgauss(Nx, prec); 
L  = sqrt(Nfft/a);
l  = floor(sqrt(-Nfft*log(prec)/(a*pi)))+1;
w  = amgauss(2*l+1,l+1,L);
H  = w/norm(w);

%% Ridges and Basins
[q,reassign,dx] = get_tfrs(s,Nfft,H,L); % get TF representations

t = (0:N-1)/Fs;
f = (-0:Nfft/2-1)*Fs/N;
t1= (0:N-1)/143000;

% Display spectrogram
FigHandle(6) =figure;
imagesc(t1,f,abs(q));%title('STFT');
set(gca,'YDir','normal');xlabel('time');ylabel('frequency'); ylim([0 200]); set(gca,'XTick',0:0.2*max(t1):max(t1)); xlim([0 max(t1)]);
xtickformat('%.2f')


%% Compute the basins of attraction associated to the different ridges/modes
[contours, rzeros, ang, dx_vec,basins, ang1] = get_contour_basins(q,dx,reassign,Ks);

%display angle dx
FigHandle(2) = figure();
subplot(2,2,[1,2])
imagesc(mod(ang,pi));
set(gca,'YDir','normal'); %title('anlges dx (modulo pi)')
subplot(2,2,3)
histogram(ang,90); %title('histogram of anlges dx')
subplot(2,2,4)
histogram(mod(ang,pi),90); %title('histogram of anlges dx (modulo pi)')

FigHandle(3) = figure();
imagesc(ang1); set(gca,'xtick',[],'ytick',[]);
set(gca,'YDir','normal');xlabel('time');ylabel('frequency'); %title('Local projecting angles (modulo pi)')

FigHandle(4) = figure();
histogram(ang1,90); %title('histogram of local projecting anlges (modulo pi)')

FigHandle(41) = figure();
subplot(2,1,1);
contour(dx_vec,[0 0]); %title('contours dx-vec')
subplot(2,1,2);
imagesc(t,f,10*log10(abs(dx_vec)));%title('signal');
set(gca,'YDir','normal');xlabel('time');ylabel('frequency');% title('amplitude of dx-vec')


% Display ridges, zeros and basins
FigHandle(5) = figure(); 
imagesc(t1,f,basins);set(gca,'clim',[-0 4],'ydir','normal');
hold on
for i=1:NumC
    plot(t1(contours{i}.x),f(contours{i}.y),'k-','linewidth',2)
end
set(gca,'YDir','normal');
tmp = find(~rzeros);[tmp1 tmp2] = ind2sub(size(rzeros),tmp);plot(t1(tmp2),f(tmp1),'r.','MarkerSize',15);
xlabel('time');ylabel('frequency');ylim([0 200]); set(gca,'XTick',0:0.2*max(t1):max(t1));xlim([0 max(t1)]);
xtickformat('%.2f')

%% Reconstruction
srec1 = get_resynth(q.*(basins==1),Nfft,H);
srec2 = get_resynth(q.*(basins==2),Nfft,H);
srec3 = get_resynth(q.*(basins==3),Nfft,H);
srec = srec1+srec2+srec3;    

% Display reconstructed signals
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.2 0.1], [0.1 0.05]);
if ~make_it_tight,  clear subplot;  end

FigHandle(7) = figure();
subplot(2,1,1)
plot(t1,srec); 
legend('total reconstruction','Location','southeast'); set(gca,'xtick',[]);ylim([-0.4 0.2]);xlim([0 max(t1)]);
subplot(2,1,2)
plot(t1,s0','r'); 
legend('noise-free','Location','southeast'); xlabel('time'); set(gca,'XTick',0:0.2*max(t1):max(t1));ylim([-0.4 0.2]);xlim([0 max(t1)]);
xtickformat('%.2f')

FigHandle(8) = figure();
%subplot(3,1,1)
plot(t1,srec1,'b',t1,srec2,'g',t1,srec3,'r'); 
legend('mode 1','mode 2','mode 3','Location','southeast'); %set(gca,'xtick',[],'ytick',[]);
xlabel('time');  set(gca,'XTick',0:0.2*max(t1):max(t1)); xlim([0 max(t1)]);
xtickformat('%.2f')

% Computes output SNRs
res = 20*log10(norm(s0)/norm(srec'-s0));
disp([num2str(res) ' dB for total output-snr']);

%%%%%%%%%%%%%%%%%%%%%% print Figures
for i = 5:8
 export_fig(FigHandle(i), ... % figure handle
     sprintf('%s/batsignal_%d', chemin0,i),... % name of output file without extension
     '-painters', ...      % renderer
     '-transparent', ...   % renderer
     '-pdf', ...           % file format
     '-r500' );             % resolution in dpi
end