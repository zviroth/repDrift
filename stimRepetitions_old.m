function res = stimRepetitions(isub,visualRegions)
if ieNotDefined('visualRegions'), visualRegions = 1; end
if ieNotDefined('isub'), isub = 1; end

tic

nsdFolder = fullfile('~','misc','data18','rothzn','nsd');
if ~isfolder(nsdFolder)
    nsdFolder = '/misc/data18/rothzn/nsd/';
end
nimgs = 10000;
nsdDesignFilename = fullfile(nsdFolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);


% nsessionsSub = [40 40 32 30 40 32 40 30];
nsplits=30;


ntrials = length(nsdDesign.masterordering);
splitImgTrials = zeros(nsplits,ntrials);
maxReps=3;
sessRepTrials = zeros(nsplits,ntrials,maxReps);%0 or 1 for each trial
trialsPerSession = 10 * 75; %12 runs, 75 trials per run
trialSession = zeros(1,ntrials);

%% DIVIDE TRIALS INTO 1st, 2nd, 3rd repetition

stimRepNum = NaN(1,3*nimgs);%for each trial, what repetition is it
imgRepTrial = zeros(nimgs,3);%the 3 trials on which this image was presented
imgRepSess = zeros(nimgs,3);%the 3 sessions on which this image was presented
for i=1:nimgs%10000 images
    ind = find(nsdDesign.masterordering == i); %trials on which this image was presented. Should be 3 trials.
    imgRepTrial(i,:) = ind;
    imgRepSess(i,:) = ceil(ind/trialsPerSession);%session number for each repetition of each image
    for irep=1:length(ind)%should be 3
        stimRepNum(ind(irep)) = irep;
    end
end

%%

saveFolder = fullfile(nsdFolder,'stimRepetitions','/');
% nsdFolder = fullfile('~','NSD');

betasfolder = fullfile(nsdFolder,['sub' num2str(isub) '_betas_func1pt8mm/']);
roifolder = betasfolder;
visualRoisFile = fullfile(roifolder,'prf-visualrois.nii');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
visRoiData = visRoiData(:);

%% LOAD BETAS



itrial=1;
for isplit=1:nsplits
    trialSession(itrial:itrial+trialsPerSession-1) = isplit;
    splitImgTrials(isplit,itrial:itrial+trialsPerSession-1) = ones;
    for irep=1:maxReps
        sessRepTrials(isplit,itrial:itrial+trialsPerSession-1,irep) = stimRepNum(itrial:itrial+trialsPerSession-1)==irep;
    end
    itrial = itrial+trialsPerSession;
end


for visualRegion=visualRegions
    switch visualRegion
        case 1
            rois=1:2;%1:2;
        case 2
            rois=3:4;
        case 3
            rois=5:6;
        case 4
            rois=7;
    end

    ['visual region: ' num2str(visualRegion)]
    for roinum=1:length(rois); roiBetas{roinum}=[]; end
    for isession=1:nsplits
        betasfilename = fullfile(betasfolder,['betas_session' num2str(isession,'%02.f') '.nii']);
        betas = niftiread(betasfilename);
        betas = cast(betas,'double');
        betas = betas/300;
%         betas = betas*100;
        betas=reshape(betas,[],size(betas,4));
        for roinum=1:length(rois)
            iroi = rois(roinum);
            roiInd{roinum} = find(visRoiData==iroi);
            roiBetas{roinum} = [roiBetas{roinum} betas(visRoiData==iroi,:)];
        end
    end
    %     [imgTrials, imgNum] = ismember(subDesign, allImgs);%logical array
    %if less than 40 sessions, only use image trials that were actually presented
    %     imgTrials = imgTrials(1:size(roiBetas{roinum},2));
    splitImgTrials = splitImgTrials(:,1:size(roiBetas{roinum},2));
    sessRepTrials = sessRepTrials(:,1:size(roiBetas{roinum},2),:);

    %     imgNum = imgNum(1:size(roiBetas{roinum},2));

    clear nvox
    for roinum=1:length(rois)
        nvox(roinum) = size(roiBetas{roinum},1);
        L = nvox(roinum);
        sessBetasSplit = NaN(nsplits,L);
        sessStdBetasSplit = NaN(nsplits,L);
        sessRepBetasSplit = NaN(nsplits,L,maxReps);
        sessRepStdBetasSplit = NaN(nsplits,L,maxReps);
        for isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            sessBetasSplit(isplit,:) = nanmean(roiBetas{roinum}(:,imgTrials>0),2);
            sessStdBetasSplit(isplit,:) = std(roiBetas{roinum}(:,imgTrials>0),0,2);
            for irep=1:maxReps
                imgTrials = sessRepTrials(isplit,:,irep);
                sessRepBetasSplit(isplit,:,irep) = nanmean(roiBetas{roinum}(:,imgTrials>0),2);%mean over images
                sessRepStdBetasSplit(isplit,:,irep) = std(roiBetas{roinum}(:,imgTrials>0),0,2);%std over images
            end
        end
        sessBetas{roinum} = sessBetasSplit;
        sessStdBetas{roinum} = sessStdBetasSplit;
        sessRepBetas{roinum} = sessRepBetasSplit;
        sessRepStdBetas{roinum} = sessRepStdBetasSplit;
        %keep all responses, with info of repetition
    end

    if visualRegion<4
        roiMeanBetas{visualRegion} = [sessBetas{1} sessBetas{2}];
        roiStdBetas{visualRegion} = [sessStdBetas{1} sessStdBetas{2}];
        roiRepBetas{visualRegion} = [sessRepBetas{1} sessRepBetas{2}];
        roiRepStdBetas{visualRegion} = [sessRepStdBetas{1} sessRepStdBetas{2}];
        combinedRoiBetas{visualRegion} = [roiBetas{1}; roiBetas{2}];
    else
        roiMeanBetas{visualRegion} = sessBetas{1};
        roiStdBetas{visualRegion} = sessStdBetas{1};
        roiRepBetas{visualRegion} = sessRepBetas{1};
        roiRepStdBetas{visualRegion} = sessRepStdBetas{1};
        combinedRoiBetas{visualRegion} = [roiBetas{1}];
    end
    %get mean betas for all repetitions of each image
    nvox = size(combinedRoiBetas{visualRegion},1);
    imgBetas{visualRegion} = NaN(nimgs,maxReps,nvox);
    imgRoiMean{visualRegion} = NaN(nimgs,maxReps);
    imgRoiCorr{visualRegion} = NaN(nimgs,3);%correlations: 1-2, 1-3, 2-3
    imgRepDist = NaN(nimgs,3);%distance in #session between repeitions 1-2, 1-3, 2-3
    
    for itrial=1:trialsPerSession*nsplits %30000
        iimg = nsdDesign.masterordering(itrial);
        irep = stimRepNum(itrial);
        imgBetas{visualRegion}(iimg,irep,:) = combinedRoiBetas{visualRegion}(:,itrial);%1/4 of the values will be NaNs
        imgRoiMean{visualRegion}(iimg,irep) = nanmean(combinedRoiBetas{visualRegion}(:,itrial),1);%mean over voxels
    end
    for iimg=1:nimgs
        imgRoiCorr{visualRegion}(iimg,:) = 1-pdist(squeeze(imgBetas{visualRegion}(iimg,:,:)),'correlation');%1-2, 1-3, 2-3
        imgRepDist(iimg,:) = [imgRepSess(iimg,2)-imgRepSess(iimg,1); imgRepSess(iimg,3)-imgRepSess(iimg,1); imgRepSess(iimg,3)-imgRepSess(iimg,2)];
    end
end
% roiBetas = combinedRoiBetas;
save(fullfile(saveFolder,['stimRepMean_sub' num2str(isub) '.mat']),...
    'roiMeanBetas','roiStdBetas','roiRepBetas','roiRepStdBetas',...
    'visualRegions',...
    'stimRepNum','combinedRoiBetas','trialSession',...
    'imgBetas','imgRoiMean','imgRoiCorr',...
    'imgRepTrial','imgRepDist',...
    'imgRepSess','splitImgTrials','sessRepTrials');

