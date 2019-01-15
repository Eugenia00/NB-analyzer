function DataMatR = CorrectCellMovement2(DataMat,Tmax,Rmax,verifyEvery,Nframestoaverage)

% DataMat has to be an N1 X N1 X N2 matrix

Size3DMShort = floor(size(DataMat,3)/Nframestoaverage);

DataMatShort = reshape(DataMat(:,:,1:Size3DMShort*Nframestoaverage),[size(DataMat,1),size(DataMat,2),Nframestoaverage,Size3DMShort]);
DataMatShort = squeeze(sum(DataMatShort,3));

% translation matrixes
%Tmax = 2;
[Tx,Ty] = meshgrid(-Tmax:Tmax);
TranslationsX = nan(size(DataMat,1),size(DataMat,2),numel(Tx));
TranslationsY = nan(size(DataMat,1),size(DataMat,2),numel(Tx));
for i = 1:numel(Tx)
    TranslationsX(:,:,i) = Tx(i);
    TranslationsY(:,:,i) = Ty(i);
end
%%
% rotation matrixes
%Rmax = 2;
rotations = pi/500*(-Rmax:1:Rmax);  % pixel shift
RotationsX = nan(size(DataMat,1),size(DataMat,2),numel(rotations));
RotationsY = nan(size(DataMat,1),size(DataMat,2),numel(rotations));
%for i = 1:numel(rotations)

if rem(size(DataMat,1),2)==1 %odd -2 -1, 0 1 2
    xint = -(size(DataMat,1)-1)/2:(size(DataMat,1)-1)/2;
else %even -3 -2 -1 0 1 2
    xint = -size(DataMat,1)/2:size(DataMat,1)/2-1;
end

if rem(size(DataMat,2),2)==1 %odd -2 -1, 0 1 2
    yint = -(size(DataMat,2)-1)/2:(size(DataMat,2)-1)/2;
else %even -3 -2 -1 0 1 2
    yint = -size(DataMat,2)/2:size(DataMat,2)/2-1;
end

[ax,ay] = meshgrid(xint,yint);
dist = sqrt(ax.^2+ay.^2);
beta = acos(ax./dist);
beta(ay<0) = 2*pi - beta(ay<0);

for i = 1:numel(rotations)
    alpha = rotations(i); % rotation angle
    
    % coordinates rotated matrix
    axr = dist.*cos(alpha+beta);
    ayr = dist.*sin(alpha+beta);
    
    % distance from original pixel
    RotationsX(:,:,i) = round(axr-ax);
    RotationsY(:,:,i) = round(ayr-ay);
end

%%
% find best tranformations
%verifyEvery = 4;
verifyFrames = 1:verifyEvery:size(DataMatShort,3);

correlations = nan(numel(Tx),numel(rotations),numel(verifyFrames));
h = fspecial('gaussian',6,4);
Im_old = imfilter(DataMatShort(:,:,1),h,'same');
[xmesh0,ymesh0] = meshgrid(1:size(DataMat,1),1:size(DataMat,2));

outXT = zeros(1,numel(verifyFrames));
outYT = zeros(1,numel(verifyFrames));
outR = zeros(1,numel(verifyFrames));

for i = 1:numel(verifyFrames)
    
    Im_now = DataMatShort(:,:,verifyFrames(i));
    Im_now_filt = imfilter(Im_now,h,'same');
    
    % try all the possible transformations
    for j = 1:size(TranslationsX,3)
        
        for k = 1:size(RotationsX,3)
            
            Mx = TranslationsX(:,:,j)+RotationsX(:,:,k);
            My = TranslationsY(:,:,j)+RotationsY(:,:,k);
            
            subindy = ymesh0+My;
            subindx = xmesh0+Mx;
            subind_out = subindx<1 | subindy<1 | ...
                subindx>size(DataMat,2) | subindy>size(DataMat,1) | ...
                isnan(subindx) | isnan(subindy);
            subindy(subind_out) = ymesh0(subind_out);
            subindx(subind_out) = xmesh0(subind_out);
            linearInd = sub2ind([size(DataMat,1),size(DataMat,2)], subindy(:), subindx(:));
            
            Im_nowR = Im_now_filt(linearInd);
            Im_nowR(subind_out) = nan;
            Im_nowR = reshape(Im_nowR,size(Im_old));
            prod = Im_nowR.*Im_old;
            correlations(j,k,i) = nansum(prod(:));
        end
    end
    
    % find the best one
    cor = correlations(:,:,i);
    [~,indmax] = max(cor(:));
    [Tra,Rot] = ind2sub(size(cor),indmax);
    if i>1
        outXT(i) = outXT(i-1) + Tx(Tra);
        outYT(i) = outYT(i-1) + Ty(Tra);
        outR(i) = outR(i-1) + rotations(Rot);
    end
    fprintf(' --- FRAME %d --- \n translation: %f, %f \n rotation: %f \n',verifyFrames(i),outXT(i),outYT(i),outR(i))
    
    % ready for next loop
    Im_old = Im_now_filt;
    
end
%%
% apply transformation
DataMatR = nan(size(DataMat));

verifyFrames_full = (1:verifyEvery*Nframestoaverage:Size3DMShort*Nframestoaverage)';

xi = (1:size(DataMat,3))';
outRFull = interp1q(verifyFrames_full,outR',xi);
outXTFull = interp1q(verifyFrames_full,outXT',xi);
outYTFull = interp1q(verifyFrames_full,outYT',xi);

nan_full = find(isnan(outRFull));
outRFull(nan_full) = outR(end);
outXTFull(nan_full) = outXT(end);
outYTFull(nan_full) = outYT(end);


for i = 1:size(DataMat,3)
    
    
    % coordinates rotated matrix
    axr = dist.*cos(outRFull(i)+beta);
    ayr = dist.*sin(outRFull(i)+beta);
    
    % distance from original pixel
    TransformationOutX = round(axr-ax + outXTFull(i));
    TransformationOutY = round(ayr-ay + outYTFull(i));
    
    % store the transformed image
    subindy = ymesh0+TransformationOutY;
    subindx = xmesh0+TransformationOutX;
    subind_out = subindx<1 | subindy<1 | ...
        subindx>size(DataMat,2) | subindy>size(DataMat,1) | ...
        isnan(subindx) | isnan(subindy);
    subindy(subind_out) = ymesh0(subind_out);
    subindx(subind_out) = xmesh0(subind_out);
    linearInd = sub2ind([size(DataMat,1),size(DataMat,2)], subindy(:), subindx(:));
    
    
    Im_now = DataMat(:,:,i);
    Im_nowR = Im_now(linearInd);
    Im_nowR(subind_out) = nan;
    Im_nowR = reshape(Im_nowR,size(Im_old));
    
    DataMatR(:,:,i) = Im_nowR;
    
    
end



end












