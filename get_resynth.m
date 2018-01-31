function srec = get_resynth(q,Nfft,H)
% Extracts a mode from the TF representation of a signal.
%
% INPUTS:
%   : Thresholded STFT of the signal
%   Offset_Freq: Number of negative frequency bins in S
%   Nfft/2: Number of total frequency bins in S
%   Nfft: total number of bins when computing the STFT (but most of the
%       negative frequencies where removed, since one only deals with real
%       signals)
%   H: window
%
% OUTPUTS
%   srec: reconstructed mode
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

n = size(q,2);
tmp = zeros(Nfft,n);
tmp(1:Nfft/2,:) = q(1:Nfft/2,:);
%tmp(Nfft+1:Nfft, :) = q(1:0, :);
[srec1,~] = tfristft(tmp,1:n,H,0);
srec = 2*real(srec1');
end