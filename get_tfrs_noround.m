function [q,reassign,dx] = get_tfrs_noround(s,Nfft,H,L)
% Extracts a mode from the TF representation of a signal.
%
% INPUTS:
%   s: signal to analyze
%   Nfft: total number of bins when computing the STFT (but most of the
%       negative frequencies where removed, since one only deals with real
%       signals)
%   H: window
%   L: window size
%
% OUTPUTS
%   q: STFT
%   reassign: coordinates to where the coefficients are reallocated (see get_tfrs)
%   dx: reassignment vector (coded as a complex number, see get_tfrs)
%   Offset_Freq: Number of negative frequency bins kept in S
%   Total_Freq: Number of total frequency bins kept in S
%
% REFERENCES:
% [1] "Adaptive multimode signal reconstruction from time-frequency representations",
%      by Sylvain Meignen, Thomas Oberlin, Philippe Depalle, Patrick Flandrin, and Stephen McLaughlin,
%      submitted.
% [2] "Time-frequency ridge analysis based on reassignment vector", 
%      by Sylvain Meignen, Tim Gardner and Thomas Oberlin, 
%      in Proceedings of the 23st European Signal Processing Conference (EUSIPCO-15), 2015.
%
% Thomas Oberlin
% 2015, July 28th
%

n = length(s);

[q,t,f] = tfrstft(s,1:n,Nfft,H,0) ;
[~,~,dx] = tfrrsp_new(s,t,Nfft,H,0) ;

% remove negative frequencies
dx(isinf(dx)) = NaN;
Offset_Freq = 0;
Total_Freq = Nfft/2 + 2*Offset_Freq;

S = q;q = zeros(Total_Freq, n);q(1+Offset_Freq:Total_Freq,:) = S(1:Total_Freq - Offset_Freq,:);q(1:Offset_Freq, :) =  S(Nfft-Offset_Freq+1:Nfft, :);
S = dx;dx = zeros(Total_Freq, n);dx(1+Offset_Freq:Total_Freq,:) = S(1:Total_Freq - Offset_Freq,:);dx(1:Offset_Freq, :) =  S(Nfft-Offset_Freq+1:Nfft, :);
clear S;
reassign = zeros(2,Total_Freq,n);
reassign(2,:,:) = imag(dx);
reassign(1,:,:) = real(dx)+Offset_Freq;
dx = (real(dx)+Offset_Freq-(1:Total_Freq)'*ones(1,n))*1i + (imag(dx)-ones(Total_Freq,1)*t);

end

