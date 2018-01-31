% This script draws figure 3 of paper [1] (see paper for more details)
% [1] A Novel Adaptive Approach to Mode Decomposition for Multicomponent
% Signals, Duong-Hung Pham, and Sylvain Meignen.
% Duong Hung PHAM 

% 2017, Jul 1st


clear all;
close all; %clc; 
set(0,'DefaultAxesFontSize',18);
chemin0 = '~/Dropbox/Papers_PHAM_MEIGNEN/EUSIPCO_2017/Latex/figures';

N = 1024; % signal length

Nx=N/2;
Nfft = 2*Nx; % number of FFT bins
Fs = N;
a=1;

s1 = 0.8*fmsin(N,.1,.4);
s1(1:50)=0;
s1(N/2-75:end)=0;

s2 = zeros(N,1);
s2(N/2-round(10*rand)) = 15;

s3 = zeros(N,1);
s3(N/2+75:N-50)= 0.8*fmconst(N/2-74-50,.32);

numC = 3; 

%% windown
prec = 10^(-6);
%[H L] = roundgauss(Nx, prec); 
L  = sqrt(Nfft/a);
l  = floor(sqrt(-Nfft*log(prec)/(a*pi)))+1;
w  = amgauss(2*l+1,l+1,L);
H  = w/norm(w);

%% Sensibility
R = 1:5; %number of realizations
SNR = 0:5:10; %four different SNR tested
Ks    = 0:2:40; 

measure_int  = zeros(length(R),length(SNR),length(Ks));
z=0;
for r = 1:length(R)
    disp([num2str(R(r)) ' number of realizations']);
    for  k = 1:length(SNR)
         
        noise1 = randn(N,1);
        noise2 = randn(N,1);
        noise3 = randn(N,1);

        s1=real(s1); sn1 = sigmerge(s1,noise1,SNR(k));
        s2=real(s2); sn2 = sigmerge(s2,noise2,SNR(k));
        s3=real(s3); sn3 = sigmerge(s3,noise3,SNR(k));
        sn=sn1+sn2+sn3;
        s = sn(:); 
        
        res = 20*log10(norm(s1+s2+s3)/norm(s-(s1+s2+s3)));
        disp(['input-SNR=' num2str(num2str(round(res)))]);

        for p = 1:length(Ks)
            disp([' Ks=' num2str(Ks(p)) ]);
            z=z+1;
            disp([num2str(z) ' Iteration']);
            [q,reassign,dx] = get_tfrs(s,Nfft,H,L); % get TF representations
            %[contours,  ~, ~, ~, ~, ~] = get_contour_basins(q,dx,reassign,Ks(p));
            [contours, rzeros, ~, ~,basins, ~] = get_contour_basins(q,dx,reassign,Ks(p));

            for i = 1:numC
                measure_int(r,k,p) = measure_int(r,k,p) + (norm(q.*(basins==i),2)).^2;
            end
            
%            measure_int(r,k,p) = contours{1}.pwr+contours{2}.pwr+contours{3}.pwr;
            
            disp([' En=' num2str(measure_int(r,k,p))]);
%             clear q;
%             clear reassign;
%             clear dx; 
%             clear contours; 
            if 0
                % Display ridges, zeros and basins
                t = (0:N-1)/Fs;
                f = (-0:Nfft/2-1)*Fs/N;
                figure();
                subplot(2,1,1)
                imagesc(abs(q));title('STFT');
                set(gca,'YDir','normal');xlabel('time');ylabel('frequency'); set(gca,'xtick',[],'ytick',[]);
                subplot(2,1,2)
                imagesc(t,f,basins);set(gca,'clim',[-0 4],'ydir','normal');
                hold on
                for i=1:numC
                    plot(t(contours{i}.x),f(contours{i}.y),'k-','linewidth',2)
                end
                set(gca,'YDir','normal');
                tmp = find(~rzeros);[tmp1 tmp2] = ind2sub(size(rzeros),tmp);plot(t(tmp2),f(tmp1),'ro');
                xlabel('time');ylabel('frequency');set(gca,'xtick',[],'ytick',[]);
                
                FigHandle2 = figure();
                subplot(3,1,1)
                imagesc(t,f,abs(q.*(basins==1)));
                set(gca,'YDir','normal');
                subplot(3,1,2)
                imagesc(t,f,abs(q.*(basins==2)));
                set(gca,'YDir','normal');
                subplot(3,1,3)
                imagesc(t,f,abs(q.*(basins==3)));
                set(gca,'YDir','normal');
                pause; close all;
            end
        end
    end
    measure = mean(measure_int,1);
end

if length(SNR)==1
    FigHandle(1) = figure()
      measure1 = squeeze(measure(1,1,:));
      plot(Ks, measure1,'LineWidth',1.5); xlabel('Paremeter Ts');ylabel('Energy of first three BAs');
      legend('SNR = 0 dB','Location','Southeast');
elseif  length(SNR)==2
    FigHandle(1) = figure();
        measure1 = squeeze(measure(1,1,:));
        measure2 = squeeze(measure(1,2,:));
        plot(Ks,measure1,'b-d', Ks, measure2,'r-*', 'LineWidth',1.5); xlabel('Paremeter Ts');ylabel('Energy of first three BAs'); 
        legend('SNR = 0 dB','SNR = 5 dB','Location','Southeast');  
else
    FigHandle(1) = figure();
        measure1 = squeeze(measure(1,1,:));
        measure2 = squeeze(measure(1,2,:));
        measure3 = squeeze(measure(1,3,:));
        plot(Ks,measure1,'r-d', Ks, measure2,'b-o', Ks, measure3,'g-*','LineWidth',1.5); xlabel('Paremeter Ts');ylabel('Energy of first three BAs'); 
        legend('SNR = 0 dB','SNR = 5 dB','SNR = 10 dB','Location','Southeast');  
end
pause
%% print Figures
export_fig(FigHandle(1), ... % figure handle
sprintf('%s/sensitivity_Ts_1', chemin0),... % name of output file without extension
'-painters', ...      % renderer
'-transparent', ...   % renderer
'-pdf', ...           % file format
'-r500' );             % resolution in dpi

% export_fig(FigHandle(1), ... % figure handle
% sprintf('%s/sensitivity_Ts_%d_SNR_%d_R_%d_Ks_%d_v2', chemin0,N,SNR(end),R(end),length(Ks)),... % name of output file without extension
% '-painters', ...      % renderer
% '-transparent', ...   % renderer
% '-pdf', ...           % file format
% '-r500' );             % resolution in dpi
save('measure_norm_SNR_0510_R_60_v4.mat','measure_int')
%end

