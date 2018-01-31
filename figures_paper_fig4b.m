% This script draws figure 3 of paper [1] (see paper for more details)
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

Ks = 10
SNR = 0;
%cas = 3;
NumC = 3; % Number of modes

sig1 = 0.8*fmsin(N,.1,.4);
sig1(1:50)=0;
sig1(N/2-75:end)=0;

sig2 = zeros(N,1);
sig2(N/2-round(10*rand)) = 15;

sig3 = zeros(N,1);
sig3(N/2+75:N-50)= 0.8*fmconst(N/2-74-50,.32);

sig1 = real(sig1);
sig2 = real(sig2);
sig3 = real(sig3);

% save('sig1.mat','sig1');save('sig2.mat','sig2');save('sig3.mat','sig3');
% save('noise1.mat','noise1');save('noise2.mat','noise2');save('noise3.mat','noise3');
load sig1; load sig2; load sig3;
load noise1; load noise2; load noise3;


s0 = sig1+sig2+sig3;

% noise1 = randn(N,1);
% noise2 = randn(N,1);
% noise3 = randn(N,1);

s1 = sigmerge(sig1,noise1,SNR);
s2 = sigmerge(sig2,noise2,SNR);
s3 = sigmerge(sig3,noise3,SNR);

s=s1+s2+s3;

disp([num2str(round(20*log10(norm(sig1)/norm(s1-sig1)))) ' dB for input SNR mode 1']);
disp([num2str(round(20*log10(norm(sig2)/norm(s2-sig2)))) ' dB for input SNR mode 2']);
disp([num2str(round(20*log10(norm(sig3)/norm(s3-sig3)))) ' dB for input SNR mode 3']);

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

% Display spectrogram
FigHandle(1) =figure;
imagesc(t,f,abs(q));%title('STFT');
set(gca,'YDir','normal');xlabel('time');ylabel('frequency');
set(gca,'xtick',[],'ytick',[]);

%% Compute the basins of attraction old method
[contours1, rzeros1, basins1] = get_contour_basins_old(q,dx,reassign);

% Display ridges, zeros and basins
FigHandle(8) = figure(); 
imagesc(t,f,basins1);set(gca,'clim',[-0 4],'ydir','normal');
hold on
for i=1:10
    plot(t(contours1{i}.x),f(contours1{i}.y),'k-','linewidth',2)
end
set(gca,'YDir','normal');
tmp = find(~rzeros1);[tmp1 tmp2] = ind2sub(size(rzeros1),tmp);plot(t(tmp2),f(tmp1),'r.','MarkerSize',15);
xlabel('time');ylabel('frequency');set(gca,'xtick',[],'ytick',[]);

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
imagesc(t,f,basins);set(gca,'clim',[-0 4],'ydir','normal');
hold on
for i=1:10
    plot(t(contours{i}.x),f(contours{i}.y),'k-','linewidth',2)
end
set(gca,'YDir','normal');
tmp = find(~rzeros);[tmp1 tmp2] = ind2sub(size(rzeros),tmp);plot(t(tmp2),f(tmp1),'r.','MarkerSize',15);
xlabel('time');ylabel('frequency');xlim([0 1]);set(gca,'XTick',0:0.2:1)%set(gca,'xtick',[],'ytick',[]);

%% Reconstruction

q1= q.*(basins==1); q4(1,:,:)=q1;
q2 =q.*(basins==2); q4(2,:,:)=q2;
q3 =q.*(basins==3); q4(3,:,:)=q3;

srec(1,:) = get_resynth(q1,Nfft,H);  
srec(2,:) = get_resynth(q2,Nfft,H);
srec(3,:) = get_resynth(q3,Nfft,H);

srec1(1,:) = srec(1,:);  
srec1(2,:) = srec(2,:);
srec1(3,:) = srec(3,:);

sigg(1,:) = sig1;
sigg(2,:) = sig2;
sigg(3,:) = sig3;

S=zeros(3,3);

for i = 1:NumC
    for j = 1:NumC
        S(i,j) = log10(norm(sigg(j,:))/norm(srec(i,:)-sigg(j,:)))>0;
    end
end
S
for i = 1:NumC
    for j = 1:NumC
        if S(i,j)~=0
            srec(i,:) = srec1(j,:);
            q5(i,:,:)=q4(j,:,:);
        end
    end
end
srec1 = srec(3,:);  q1 = squeeze(q5(1,:,:));
srec2 = srec(1,:);  q2 = squeeze(q5(2,:,:));
srec3 = srec(2,:);  q3 = squeeze(q5(3,:,:));

srec = srec1+srec2+srec3;    

% Display reconstructed signals

FigHandle(6) = figure();
subplot(3,1,1)
imagesc(t,f,abs(q1));%title('basin mode 1');
set(gca,'YDir','normal');xlabel('time');ylabel('frequency');
subplot(3,1,2)
imagesc(t,f,abs(q2));%title('basin mode 2');
set(gca,'YDir','normal');xlabel('time');ylabel('frequency');
subplot(3,1,3)
imagesc(t,f,abs(q3));%title('basin mode 3');
set(gca,'YDir','normal');xlabel('time');ylabel('frequency');

% make_it_tight = true;
% subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.1 0.1], [0.1 0.01]);
% if ~make_it_tight,  clear subplot;  end

FigHandle(7) = figure();
%set(FigHandle(7),'units','normalized','outerposition',[0 0 1 1]);

%subplot(4,1,1)
plot(t,[srec;s0']); xlim([0.3 0.7])
legend('reconstructed','noise-free'); xlabel('time');% set(gca,'xtick',[]);

% subplot(4,1,2)
% plot(t,[srec1;sig1']); xlim([0.3 0.7])
% legend('reconstructed mode','noise-free mode 1'); set(gca,'xtick',[]);
% 
%  subplot(4,1,3)
%  plot(t,[srec2;sig2']);xlim([0.3 0.7])
%  legend('reconstructed mode','noise-free mode 2'); set(gca,'xtick',[]);
% 
%  subplot(4,1,4)
% plot(t,[srec3;sig3']);xlim([0.3 0.7])
% legend('reconstructed mode','noise-free mode 3'); set(gca,'xtick',[]);

% Computes output SNRs
res = 20*log10(norm(s0)/norm(srec'-s0));
disp([num2str(res) ' dB for total output-snr']);
disp([num2str(20*log10(norm(sig1)/norm(srec1'-sig1))) ' dB for mode 1']);
disp([num2str(20*log10(norm(sig2)/norm(srec2'-sig2))) ' dB for mode 2']);
disp([num2str(20*log10(norm(sig3)/norm(srec3'-sig3))) ' dB for mode 3']);

%%%%%%%%%%%%%%%%%%%%%% print Figures
for i = 5
 export_fig(FigHandle(i), ... % figure handle
     sprintf('%s/3modes_%d_failing', chemin0,i),... % name of output file without extension
     '-painters', ...      % renderer
     '-transparent', ...   % renderer
     '-pdf', ...           % file format
     '-r500' );             % resolution in dpi
end