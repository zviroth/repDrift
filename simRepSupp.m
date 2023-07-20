function res = simRepSupp(isub,iroi,repSuppFactor,toZscore)

% simulate repetition suppression. Uses pRF and tuning for a single subject
% and ROI.
% Simulate responses based on the steerable pyramid outputs. For each
% repetition decrease the response by a fixed factor. Then run simulated
% responses through model, to see if we see representational drift.


%set parameters and constants:
tic
if ieNotDefined('toZscore'), toZscore = 0; end
if ieNotDefined('toZscore'), repSuppFactor = 0.7; end
noiseAmp = 0.1;%1/SNR. Signal = STD of single voxel responses. Noise = STD of gaussian noise added. 
% iroi = 1;
% isub=1;

% nsess=30;
nsdFolder = fullfile('~','misc','data18','rothzn','nsd');
if ~isfolder(nsdFolder)
    nsdFolder = '/misc/data18/rothzn/nsd/';
end

%% load model fit for each voxel (for first session)
zscoreStr='';
saveFolder = fullfile(nsdFolder,'repDrift_expand/');
% load(fullfile(saveFolder,['regressPrfSplit_session_v' num2str(visualRegion) '_sub' num2str(isub) zscoreStr '.mat']), ...
%         'nsd', ...%'synth',...
%         'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits');
load([saveFolder 'regressSessCombineROI_sub' num2str(isub) zscoreStr '.mat'],'sessBetas',...
    'voxOriCoef','voxCoef');
nsessions = size(sessBetas{1},1);
%% load coefficients for each voxel in the ROI, for all images
prffolder = '/misc/data18/rothzn/nsd/prfsample_expand/';
visualRegion=iroi;
load(fullfile(prffolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
    'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
    'roiPrf');
if iroi<4%ventral and dorsal ROIs
    prfSampleOri = cat(2,prfSampleLevOri{1},prfSampleLevOri{2});%img, vox, lev
    prfSample = cat(2,prfSampleLev{1},prfSampleLev{2});
else
    prfSampleOri = prfSampleLevOri{1} ;%img, vox, lev
    prfSample = prfSampleLev{1};
end

%% simulate responses for first repetition of each image
nsdDesignFilename = fullfile(nsdFolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
subDesign = nsdDesign.subjectim(isub,nsdDesign.masterordering);%for each of 30000 trials, what is corresponding image (out of 73000 images)
% allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
% allImgs = unique(allImgs);
nimgs=10000;
simImgs = 1:nimgs;

nvox = size(prfSample,2);

% voxSessResp = zeros(nvox,length(simImgs));
voxSessRespOri = zeros(nvox,length(simImgs));

for ivox=1:nvox
    voxPrfSample = squeeze(prfSample(simImgs,ivox,:));
    voxPrfSample(:,end+1) = ones;

    voxPrfSampleOri = squeeze(prfSampleOri(simImgs,ivox,:,:));
    voxPrfSampleOri = reshape(voxPrfSampleOri,[],numLevels*numOrientations);
    %add highest and lowest SF
    voxPrfSampleOri(:,end+1:end+2) = squeeze(prfSample(simImgs,ivox,numLevels+1:numLevels+2));
    voxPrfSampleOri(:,end+1) = ones;

    %     isess=1;%only using the fit for session 1
    %     voxSessCoef = squeeze(voxCoef{iroi}(isess,ivox,:));
    %     voxSessOriCoef = squeeze(voxOriCoef{iroi}(isess,ivox,:));

    %OR: use the mean across sessions:
    voxSessCoef = squeeze(mean(voxCoef{iroi}(:,ivox,:),1));
    voxSessOriCoef = squeeze(mean(voxOriCoef{iroi}(:,ivox,:),1));

    voxSessResp(ivox,:) =  voxPrfSample(:,:)*voxSessCoef(:);
    voxSessRespOri(ivox,:) =  voxPrfSampleOri(:,:)*voxSessOriCoef(:);
end

%% get image and repetition for each trial, from analyzeRepetitions.m
stimRepNum = NaN(1,3*nimgs);%for each trial, what repetition is it
imgRepTrial = zeros(nimgs,3);%the 3 trials on which this image was presented
for i=1:nimgs%10000 images
    ind = find(nsdDesign.masterordering == i); %trials on which this image was presented. Should be 3 trials.
    imgRepTrial(i,:) = ind;
    for irep=1:length(ind)%should be 3
        stimRepNum(ind(irep)) = irep;
    end
end

toc


%% scale pyramid output by each trial's repetition number
% repSuppFactor = 0.7;
scalingAmp = [1 repSuppFactor repSuppFactor];
trialScaling = scalingAmp(stimRepNum);%scaling for each of 30,000 trials
%need to find the original amplitude for each trial
%what is the img for each of 30,000 trials?
% subDesign - but this is out of 73,000 images. How do we get it out of the
% 10,000 images listed in

[imgTrials, imgNum] = ismember(subDesign, allImgs);%logical array
%imgNum is the image number for each trial, out of 10,000 images
unscaledAmp = voxSessRespOri(:,imgNum);%for each trial the original simulated amplitude
%multiple each column by that trial's scaling
scaledAmp = bsxfun(@times,trialScaling,unscaledAmp);

%% add noise
normalNoise = randn(size(scaledAmp));
voxelStd = std(scaledAmp,0,2);
scaledNoise = bsxfun(@times,voxelStd*noiseAmp,normalNoise);%scaled by voxel STD
roiBetas = scaledAmp + scaledNoise;

% keyboard

%% zscore simulated betas
nsplits = nsessions;
ntrials = length(nsdDesign.masterordering);
splitImgTrials = zeros(nsplits,ntrials);
trialsPerSession = 10 * 75; %12 runs, 75 trials per run
itrial = 1;
for isplit=1:nsplits
    splitImgTrials(isplit,itrial:itrial+trialsPerSession-1) = ones;
    itrial = itrial+trialsPerSession;
end

for isplit=1:nsplits
    imgTrials = splitImgTrials(isplit,:);
    if toZscore==1
        roiBetas(:,imgTrials>0) = zscore(roiBetas(:,imgTrials>0),0,2);
    elseif toZscore==2%only remove mean
        roiBetas(:,imgTrials>0) = roiBetas(:,imgTrials>0)-mean(roiBetas(:,imgTrials>0),2);
    elseif toZscore==3%only normalize variance (subtract mean, zscore, add back mean)
        roiBetas(:,imgTrials>0) = mean(roiBetas(:,imgTrials>0),2)+zscore(roiBetas(:,imgTrials>0)-mean(roiBetas(:,imgTrials>0),2),0,2);
    elseif toZscore==4%only remove ROI mean
        roiBetas(:,imgTrials>0) = roiBetas(:,imgTrials>0)-mean(mean(roiBetas(:,imgTrials>0),2));
    elseif toZscore==0
        %nothing
    end
end

%% train model on each session, test model on each session
voxCoefSplit = nan(nsplits,nvox,size(voxSessCoef,1));
voxOriCoefSplit =  nan(nsplits,nvox,size(voxSessOriCoef,1));
parfor isplit=1:nsplits
    imgTrials = splitImgTrials(isplit,:);
    numTrials = sum(imgTrials);
    sessBetasSplit(isplit,:) = mean(roiBetas(:,imgTrials>0),2);
    sessStdBetasSplit(isplit,:) = std(roiBetas(:,imgTrials>0),0,2);
    for ivox=1:nvox%nvox(roinum)
        voxBetas = roiBetas(ivox,imgTrials>0)';

        voxPrfSample = squeeze(prfSample(imgNum(imgTrials>0),ivox,:));%includes lowest and highest SF
        %add constant predictor
        voxPrfSample(:,end+1) = ones;
        %voxePrfSample*coeff = voxBetas;
        voxCoefSplit(isplit,ivox,:) = voxPrfSample\voxBetas;

        voxPrfOriSample = squeeze(prfSampleOri(imgNum(imgTrials>0),ivox,:,:));
        voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);

        %add lowest and highest SF
        voxPrfOriSample(:,end+1:end+2) = voxPrfSample(:,end-2:end-1);

        %add constant predictor
        voxPrfOriSample(:,end+1) = ones;

        %perform regression
        voxOriCoefSplit(isplit,ivox,:) = voxPrfOriSample\voxBetas;

        %compute within-split residuals
        voxOriResidual = voxBetas' - squeeze(voxOriCoefSplit(isplit,ivox,:))'*voxPrfOriSample';
        voxResidual = voxBetas' - squeeze(voxCoefSplit(isplit,ivox,:))'*voxPrfSample';

        %r2 within split
        r2WithinSplit(isplit,ivox) = rsquared(voxResidual(1:sum(splitImgTrials(isplit,:))), roiBetas(ivox,imgTrials>0));
        r2oriWithinSplit(isplit,ivox) = rsquared(voxOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas(ivox,imgTrials>0));
    end
end

%%cross-session generalization
parfor isplit=1:nsplits
    imgTrials = splitImgTrials(isplit,:);
    numTrials = sum(imgTrials);

    for ivox=1:nvox%nvox(roinum)

        voxBetas = roiBetas(ivox,imgTrials>0)';

        voxPrfSample = squeeze(prfSample(imgNum(imgTrials>0),ivox,:));%includes lowest and highest SF
        %add constant predictor
        voxPrfSample(:,end+1) = ones;
        %voxePrfSample*coeff = voxBetas;

        voxPrfOriSample = squeeze(prfSampleOri(imgNum(imgTrials>0),ivox,:,:));
        voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);

        %add lowest and highest SF
        voxPrfOriSample(:,end+1:end+2) = voxPrfSample(:,end-2:end-1);

        %add constant predictor
        voxPrfOriSample(:,end+1) = ones;

        %%
        %cross-split goodness-of-fit (r-squared and correlation)
        %taking coefficients from nextSplit, multiplying by pRF energy in isplit, and comparing to isplit betas
        for nextSplit=1:nsplits
            voxOriResidualSplit = voxBetas' - squeeze(voxOriCoefSplit(nextSplit,ivox,:))'*voxPrfOriSample';
            voxResidualSplit = voxBetas' - squeeze(voxCoefSplit(nextSplit,ivox,:))'*voxPrfSample';

            %r2 between splits
            r2splitVox(isplit,ivox,nextSplit) = rsquared(voxResidualSplit(1:sum(splitImgTrials(isplit,:))), roiBetas(ivox,imgTrials>0));
            r2OrisplitVox(isplit,ivox,nextSplit) = rsquared(voxOriResidualSplit(1:sum(splitImgTrials(isplit,:))), roiBetas(ivox,imgTrials>0));

            %root sum of squares between splits
            rssSplitVox(isplit,ivox,nextSplit) = sum(voxResidualSplit(1:sum(splitImgTrials(isplit,:))).^2);
            rssOrisplitVox(isplit,ivox,nextSplit) = sum(voxOriResidualSplit(1:sum(splitImgTrials(isplit,:))).^2);

            %corr between splits
            pearsonRoriVox(isplit,ivox,nextSplit) = corr(voxBetas,(squeeze(voxOriCoefSplit(nextSplit,ivox,:))'*voxPrfOriSample')');
            pearsonRVox(isplit,ivox,nextSplit) = corr(voxBetas,(squeeze(voxCoefSplit(nextSplit,ivox,:))'*voxPrfSample')');
        end

    end
end

%%

saveFolder = fullfile(nsdFolder,'stimRepetitions','/');

zscoreStr='';
if toZscore==1
    zscoreStr = '_zscore';
elseif toZscore==2
    zscoreStr = '_zeroMean';
elseif toZscore==3
    zscoreStr = '_equalStd';
elseif toZscore==4
    zscoreStr = '_zeroROImean';
end

repSuppFactorStr = num2str(repSuppFactor,'%01.1f');
repSuppFactorStr(2)=[];

save(fullfile(saveFolder,['simRepSupp_' repSuppFactorStr '_sub' num2str(isub) '_v' num2str(iroi) zscoreStr '.mat']),...
    'r2OrisplitVox','rssOrisplitVox','pearsonRoriVox',...
    'r2oriWithinSplit',...
    'trialScaling','scalingAmp','noiseAmp');
toc
end
%%
function r2 = rsquared(Xresid, Xorig)
r2 = 1 - sum((Xresid).^2)/sum((Xorig - mean(Xorig)).^2);
end