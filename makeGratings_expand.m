
close all
clear all
tic
interpImgSize = 714;
backgroundSize = 1024;
interpDeg = 8.4;
imgScaling = 0.5;
width = backgroundSize*imgScaling;
pixPerDeg = interpImgSize*imgScaling/interpDeg;
totalDeg = width/pixPerDeg; %including gray background

numFreqs = 30;
numAngles = 32;
minCycles = 1;
maxCycles = backgroundSize*imgScaling/2;
minCpd = minCycles/totalDeg;
maxCpd = maxCycles/totalDeg;
cpds = logspace(log10(minCpd),log10(maxCpd),numFreqs);

angles = linspace(0,180,numAngles+1);
angles = angles(1:numAngles);







% from mglMakeGrating.m

phase = 0;
phase = pi*phase/180;



% get a grid of x and y coordinates that has
% the correct number of pixels
x = -(width-1)/2:(width-1)/2;
y = -(width-1)/2:(width-1)/2;
[xMesh,yMesh] = meshgrid(x,y);


%%
% cpd = 1;
% angle = 45;%0 to 180
% sf = cpd/pixPerDeg;
% totalCycles = sf*width;
% cpd = totalCycles/totalDeg;%(sf*width)/(width/pixPerDeg) = sf*pixPerDeg

freqs = cpds/pixPerDeg;
allGratings = zeros(length(freqs),length(angles),width,width);

for isf=1:length(freqs)
    sf = freqs(isf);
    for iangle=1:length(angles)
        angle = angles(iangle);
        
        % calculate orientation
        angle = pi*angle/180;
        
        a=cos(angle)*sf*2*pi;
        b=sin(angle)*sf*2*pi;
        % compute grating
        im = cos(a*xMesh+b*yMesh+phase);
        allGratings(isf,iangle,:,:) = im;
%         imagesc(im)
        
        
        
        %% pass image through steerable pyramid
        normResp=0;
        % construct quad frequency filters
        numOrientations = 8;
        bandwidth = 1;
        dims = [backgroundSize backgroundSize];
        dims = dims*imgScaling;
        numLevels = maxLevel(dims,bandwidth);
%         [freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);
        [freqResps, pind] = makeSteerFRs(dims, numLevels, numOrientations, bandwidth);
        
        %%
%         [pyr, pind] = buildQuadBands(im, freqRespsImag, freqRespsReal);
        [pyr, pind] = buildSteerBands(im, freqResps);
        sumOri = cell(numLevels,1);
        modelOri = cell(numLevels,1);
        for ilev = 1:numLevels
            % loop over levels and orientations of the pyramid
            % initialize output
            sumOri{ilev}(:,:) = zeros(dims(1), dims(2),'single');
            modelOri{ilev} = zeros(numOrientations, dims(1), dims(2),'single');
            for orientation = 1:numOrientations
                if normResp
                    nEnergies = normEnergies(pyr,pind,numOrientations,0.1);
                    thisBand = abs(accessSteerBand(nEnergies,pind,numOrientations,ilev,orientation));
                else
%                     thisBand = abs(accessSteerBand(pyr, pind, numOrientations,ilev, orientation)).^2;
                    thisBand = abs(accessSteerBand(pyr, pind, numOrientations,ilev, orientation));
                end
                sumOri{ilev}(:,:) = sumOri{ilev}(:,:) + thisBand;
                modelOri{ilev}(orientation,:,:) = thisBand;
                
                %sum energy
                sumOriEnergy(isf,iangle,ilev) = sum(sumOri{ilev}(:));
                modelOriEnergy(isf,iangle,ilev,orientation) = sum(thisBand(:));
            end
            %     sumOri{iLev}(:,:) = temp;
        end
        %ADD LOWEST SPATIAL FREQUENCY
        thisBand = abs(pyrLow(freqResps,pind));
        sumOri{numLevels+1} = thisBand;
        modelOri{numLevels+1} = thisBand;
        
        %ADD HIGHEST SPATIAL FREQUENCY
        thisBand = abs(pyrHi(freqResps,pind));
        sumOri{numLevels+2} = thisBand;
        modelOri{numLevels+2} = thisBand;
        
        
    end
end
toc
saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand');
if ~isfolder(saveFolder)
    saveFolder = '/misc/data18/rothzn/nsd/repDrift_expand/';
end

save(fullfile(saveFolder,'gratings','gratings.mat'),'cpds','angles','freqs','numOrientations','numLevels',...
    'sumOriEnergy','modelOriEnergy','normResp','backgroundSize','imgScaling');
save(fullfile(saveFolder,'gratings','allGratings.mat'),'cpds','angles','freqs','numOrientations','numLevels',...
    'sumOriEnergy','modelOriEnergy','normResp','backgroundSize','imgScaling','allGratings','-v7.3');