function [contours, rzeros, ang, dx_vec,basins, ang1] = get_contour_basins(q,dx,reassign,Ks)
% Estimates the ridges of a given short-time Fourier transform, and the
% corresponding "basins of attraction", from the reassignment vector.
%
% INPUTS:
%   q: STFT of the signal
%   dx: reassignment vector (coded as a complex number, see get_tfrs)
%   reassign: coordinates to where the coefficients are reallocated (see get_tfrs)
%
% OUTPUTS
%   contours: cell which contains the contours, sorted by descending power. 
%       contours{i} is a 2 by n array, which contains the coordinates of the 
%       ith contour in q (of size n)
%   rzeros: matrix, same size as q, contains the mask of the zeros of q.
%   basins: matrix, same size as q, contains the labels of the basins
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

[fM,fN] = size(q);
absq = abs(q);

%% orientation modulo pi of the reassignment vector
ang = angle(dx);
ang1 = angle(dx);
ang1 = mod(ang1,pi);
ang1=round(180.*ang1./(pi))*pi/180;
ang2 = zeros(fM,fN);

for k =1:fM
    for p=1:fN
    A = ang1(max(1,k-Ks):min(fM,k+Ks),max(1,p-Ks):min(fN,p+Ks));
    ang2(k,p) = mode(A(:));   
    end
end

ang1 =ang2;
% for i = 1:Ks/2:floor(fM/Ks)*Ks
%   for j = 1:Ks:floor(fN/Ks)*Ks
%     tmp= ang1(i:(i+Ks/2-1), j:(j+Ks-1));  
%     a = unique(tmp);
%     a(isnan(a))=[];
%     tmp(isnan(tmp))=[];
%     out = [a,histc(tmp(:),a)];
%     %out=round(4.*out./(pi))*pi/4;
%     out = flipud(sortrows(out,2));
%     %figure();histogram(tmp,180); out(1,1)
%     ang1(i:(i+Ks/2-1), j:(j+Ks-1)) = out(1,1);
%     %pause;close;
%   end
% end

dx_vec = real(dx).*cos(ang1)+imag(dx).*sin(ang1);
dx_vec = round(dx_vec);

% figure()
% subplot(2,2,[1,2])
% imagesc(mod(angle(dx),pi));
% set(gca,'YDir','normal'); title('anlges dx (modulo pi)')
% subplot(2,2,3)
% histogram(ang,36); title('histogram of anlges dx')
% subplot(2,2,4)
% histogram(mod(ang,pi),36); title('histogram of anlges dx (modulo pi)')
% 
% figure();
% subplot(2,1,1)
% imagesc(ang1); title('signal');
% set(gca,'YDir','normal');xlabel('time');ylabel('frequency'); title('Local projecting angles (modulo pi)')
% 
% subplot(2,1,2)
% histogram(ang1,36); title('histogram of local projecting anlges (modulo pi)')
% 
% figure()
% subplot(2,1,1);
% contour(dx_vec,[0 0]); title('contours dx-vec')
% subplot(2,1,2);
% imagesc(abs(log10(abs(dx_vec))));title('signal');
% set(gca,'YDir','normal');xlabel('time');ylabel('frequency'); title('amplitude of dx-vec')


%% Create contour segmentation mask
%[rzeros] = ConSegMaskCreate(abs(dx),fM,fN);% Garnder's solution
% better: Patrick's zeros detection
[xZ,yZ] = extr2minth(absq,max(absq(:))/10^14) ; % get 0's coordinates
xZ(xZ<2)=[];yZ(yZ<2)=[];xZ(xZ>fM-1)=[];yZ(yZ>fN-1)=[];
rzeros = ones(size(q));
for i=1:length(xZ)
    rzeros(xZ(i)+(-1:1),yZ(i)+(-1:1)) = 0;
end


%% Contour extraction

% Extract 0 level-sets
cont = contourc(dx_vec,[0 0]);

% some inits
k = 1; contours = {};
cpt = 0;
p = length(cont);
while k<=p
    % get full contours
    if cont(1,k) ~=0
        error('should be 0 ... operation stopped');
    end
    l = cont(2,k);
    px = cont(1,k+1:k+l);
    py = cont(2,k+1:k+l);
    k = k+l+1;
    
    % segment contours according to the zeros of the spectrogramm
    idx = find(px<1 | px>fN | py<1 | py>fM);
    px(idx)=[];py(idx)=[];
    % zeros on the ridge
    idx = [1 find(~rzeros(sub2ind([fM fN],round(py),round(px)))) length(px)];
    
    % remove doublons
    idx(diff(idx)==1) = [];
    
    for j=1:length(idx)-1
        cpt = cpt+1;
        % one contour
        ppx = px(idx(j):idx(j+1));
        ppy = py(idx(j):idx(j+1));
        % compute length
        contours{cpt}.len = sqrt(sum(diff(ppx).^2 + diff(ppy).^2));
        % round
        ppx = round(ppx);ppy = round(ppy);
        contours{cpt}.x = ppx;
        contours{cpt}.y = ppy;
        % compute power
        tmp = absq(sub2ind([fM fN],ppy,ppx));
        contours{cpt}.pwr = sum(tmp.^2);
    end
end
numC = cpt;

% sort contours by descending power
pwr = zeros(numC,1);
for k=1:numC
    pwr(k) = contours{k}.pwr;
end
[~,idx] = sort(pwr,'descend');
contours = contours(idx);
NumC = length(contours);

%% Basin computation 
TempCon = NaN*zeros(fM,fN);

% Mask with contours labeled with integers, sorted by descending power
for i = 1:NumC
    pos = fM*(contours{i}.x-1)+contours{i}.y;
    TempCon(pos) = i;
end
reassign = round(reassign);
basins = zeros(fM,fN);
% run through TF points
for i = 1:fM,
    for j = 1:fN,
        k = reassign(1,i,j); %row index
        p = reassign(2,i,j); %column index
        % "Majority vote": the point is affected to the ridge with is most
        % present inside a 11x11 window
        tmp = TempCon(max(1,k-2):min(fM,k+2),max(1,p-2):min(fN,p+2));
        basins(i,j) = mode(tmp(:));
    end
end















