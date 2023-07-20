%use only trials of a single repetition (default: first presentation of
%each image).
%equalize number of trials across sessions by randomly downsampling trials.

%

function res = regressPrfSplit_1rep(isub,visualRegions,onlyRep,nsessions,toZscore,noBackground)
if ieNotDefined('onlyRep'), onlyRep = 1; end%only use trials of this repetition. 1,2,3
if ieNotDefined('nsessions'), nsessions = 30; end%will determine the minimum number of images for this repetition in each session
if ieNotDefined('visualRegions'), visualRegions = 1:4; end
if ieNotDefined('toZscore'), toZscore = 0; end
if ieNotDefined('noBackground'), noBackground = 0; end

maxReps=3;

%uses file from prfSampleModel_expand.m 
close all
% clear all
tic


bandpass = 1; bandMin = 1; bandMax = 7;
noBackgroundStr = '';
if noBackground
    bandMax=6;
    noBackgroundStr = '_noBackground';
end

toPerm = 0;



%%%%%%%%%%%%%%%%%%%%%%%%
% nsessionsSub = [40 40 32 30 40 32 40 30];
% nsessions=nsessionsSub(isub);
nsplits=nsessions;


boxfolder = ['/misc/data18/rothzn/nsd/prfsample_expand/'];
saveFolder = ['/misc/data18/rothzn/nsd/repDrift_expand/'];


betasfolder = ['/misc/data18/rothzn/nsd/sub' num2str(isub) '_betas_func1pt8mm/'];
% stimfilename = fullfile(folder,'nsdsynthetic_colorstimuli_subj01.hdf5');
nsdfolder = '/misc/data18/rothzn/nsd/';
roifolder = ['/misc/data18/rothzn/nsd/sub' num2str(isub) '_betas_func1pt8mm/'];
synthfolder = ['/misc/data18/rothzn/nsd/sub' num2str(isub) '_synth_func1pt8mm/'];
visualRoisFile = fullfile(roifolder,'prf-visualrois.nii');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
visRoiData = visRoiData(:);



nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
nsdDesign.masterordering;%for each of 30000 trials, what is corresponding image (out of 10000 images)
subDesign = nsdDesign.subjectim(isub,nsdDesign.masterordering);%for each of 30000 trials, what is corresponding image (out of 73000 images)


% sum([1 2 3]*sum(nsdDesign.basiccnt,2)); %30,000
% sum([1 2 3]*nsdDesign.basiccnt,1); %750

%% get image and repetition for each trial, from analyzeRepetitions.m
nimgs=10000;
trialsPerSession = 10 * 75; %12 runs, 75 trials per run
stimRepNum = NaN(1,maxReps*nimgs);%for each trial, what repetition is it
imgRepTrial = zeros(nimgs,maxReps);%the 3 trials on which this image was presented
% imgRepSess = zeros(nimgs,maxReps);%the 3 sessions on which this image was presented
for i=1:nimgs%10000 images
    ind = find(nsdDesign.masterordering == i); %trials on which this image was presented. Should be 3 trials.
    imgRepTrial(i,:) = ind;
%     imgRepSess(i,:) = ceil(ind/trialsPerSession);%session number for each repetition of each image
    for irep=1:length(ind)%should be 3
        stimRepNum(ind(irep)) = irep;
    end
end

%% divide trials into nsplits, and also get repetition for each trial
ntrials = length(nsdDesign.masterordering);
splitImgTrials = zeros(nsplits,ntrials);

sessRepTrials = zeros(nsessions,ntrials,maxReps);%0 or 1 for each trial
%     trialsPerSplit = trialsPerSession * 40;%floor or ceil?
itrial = 1;
for isplit=1:nsplits
    splitImgTrials(isplit,itrial:itrial+trialsPerSession) = ones;
    for irep=1:maxReps
        sessRepTrials(isplit,itrial:itrial+trialsPerSession-1,irep) = stimRepNum(itrial:itrial+trialsPerSession-1)==irep;
    end
    itrial = itrial+trialsPerSession;
end
%% keep only trials of the chosen repetition
splitImgTrials = splitImgTrials .* sessRepTrials(:,:,onlyRep); % sum(splitImgTrials')

%% find session with minimal number of this repetition
% h=histogram(imgRepSess(:,onlyRep));
% [minNumRep minSessRep] = min(h.Values(1:nsessions));
[minNumRep minSessRep] = min(sum(sessRepTrials(:,:,onlyRep)'));
% [minNumRep minSessRep] = min(sum(splitImgTrials'));

%% keep only the minimal number of trials from each session
for isess=1:nsessions
    temp = splitImgTrials(isess,(splitImgTrials(isess,:)==1));
    randTrialNum = randperm(size(temp,2));
    temp(randTrialNum>minNumRep) = 0;%remove excess trials
    splitImgTrials(isess,(splitImgTrials(isess,:)==1)) = temp;
end

[minNumRep minSessRep]
% rep 1, nsessions 30: 127
% rep 1, nsessions 25: 186 (sess 23)
% rep 1, nsessions 20: 230
% rep 1, nsessions 15: 285 (sess 14)
% rep 1, nsessions 12 (things): 358
% rep 1, nsessions 10: 363
% rep 2: 136 (sess 1)
% rep 3: 31 (sess 1)


% keyboard






for visualRegion=visualRegions
    ['visual region: ' num2str(visualRegion)]
    
    load(fullfile(boxfolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
        'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
        'roiPrf');
%     %if prf sampling was done with the nonlinear CSS prf, then we want to
%     %define the weights for the constrained model as a sum of the
%     %orientation model across orientations:
%     for roinum=1:length(rois)
%         prfSampleLev{roinum} = squeeze(sum(prfSampleLevOri{roinum},4));
%     end
    
    
    
    
%     if bandpass
%         for roinum=1:length(rois)
%             iroi = rois(roinum);
%             prfSampleLevOri{roinum} = prfSampleLevOri{roinum}(:,:,bandMin:bandMax,:);
%             prfSampleLev{roinum} = prfSampleLev{roinum}(:,:,bandMin:bandMax);
%         end
%         numLevels = bandMax-bandMin+1;
%     end
    
    for roinum=1:length(rois); iroi = rois(roinum); roiBetas{roinum}=[]; end
    for isession=1:nsessions
        betasfilename = fullfile(betasfolder,['betas_session' num2str(isession,'%02.f') '.nii']);
        betas = niftiread(betasfilename);
        betas = cast(betas,'double');
        betas = betas/300;
        betas=reshape(betas,[],size(betas,4));
        for roinum=1:length(rois)
            iroi = rois(roinum);
            if toZscore==1
                roiBetas{roinum} = [roiBetas{roinum} zscore(betas(visRoiData==iroi,:),0,2)];
            elseif toZscore==2%only remove mean
                roiBetas{roinum} = [roiBetas{roinum} betas(visRoiData==iroi,:)-mean(betas(visRoiData==iroi,:),2)];
            elseif toZscore==3%only normalize variance (subtract mean, zscore, add back mean)
                roiBetas{roinum} = [roiBetas{roinum} mean(betas(visRoiData==iroi,:),2)+zscore(betas(visRoiData==iroi,:)-mean(betas(visRoiData==iroi,:),2),0,2)];
            elseif toZscore==4%only remove ROI mean
                roiBetas{roinum} = [roiBetas{roinum} betas(visRoiData==iroi,:)-mean(mean(betas(visRoiData==iroi,:),2))];
            elseif toZscore==0
                roiBetas{roinum} = [roiBetas{roinum} betas(visRoiData==iroi,:)];
            end
            roiInd{roinum} = find(visRoiData==iroi);
        end
    end
    
    if toPerm
        for roinum=1:length(rois)
            %         roiBetas{roinum} = roiBetas{roinum}(:,randperm(size(roiBetas{roinum},2)));
            for ivox=1:size(roiBetas{roinum},1)
                roiBetas{roinum}(ivox,:) = roiBetas{roinum}(ivox,randperm(size(roiBetas{roinum},2)));
            end
        end
    end
    
    [imgTrials, imgNum] = ismember(subDesign, allImgs);%logical array
    %if less than 40 sessions, only use image trials that were actually presented
%     imgTrials = imgTrials(1:size(roiBetas{roinum},2));
    splitImgTrials = splitImgTrials(:,1:size(roiBetas{roinum},2));
    imgNum = imgNum(1:size(roiBetas{roinum},2));
%     imgTrialSum = cumsum(imgTrials);
%     splitImgTrials = repmat(imgTrials,nsplits,1);
    maxNumTrials = max(sum(splitImgTrials,2));
    
    
    
    r2 = cell(length(rois),1);
    r2ori = cell(length(rois),1);
    r2oriSplit = cell(length(rois),1);
    r2split = cell(length(rois),1);
    voxOriResidOriR2 = cell(length(rois),1);
    voxOriPredOriR2 = cell(length(rois),1);
    voxResidOriR2 = cell(length(rois),1);
    voxPredOriR2 = cell(length(rois),1);
    sessBetas = cell(length(rois),1);
    sessStdBetas = cell(length(rois),1);
    pearsonRori = cell(length(rois),1);
    pearsonR = cell(length(rois),1);
    clear nvox
    for roinum=1:length(rois)
        nvox(roinum) = size(roiBetas{roinum},1);
        voxOriResidual{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxResidual{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        %         voxOriResidualSplit{roinum} = NaN(nsplits, nsplits, nvox(roinum),maxNumTrials);
        %         voxResidualSplit{roinum} = NaN(nsplits, nsplits, nvox(roinum),maxNumTrials);
        voxOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels+1);
        voxPredOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxPredCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels+1);
        voxOriPredOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        
        
        voxResidOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxOriResidOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        
        %get model coefficients for each voxel, within each split
        L = nvox(roinum);
        voxCoefSplit = zeros(nsplits,L,numLevels+1+2);%including constant and highest & lowest SF
        voxOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1+2);%including constant and highest & lowest SF
        voxPredOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1+2);
        voxOriPredOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1+2);
        voxResidOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1+2);
        voxOriResidOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1+2);
        r2WithinSplit = zeros(nsplits,L);
        r2oriWithinSplit = zeros(nsplits,L);
        voxOriResidOriR2Split = zeros(nsplits,L);
        voxOriPredOriR2Split = zeros(nsplits,L);
        voxResidOriR2Split = zeros(nsplits,L);
        voxPredOriR2Split = zeros(nsplits,L);
        r2splitVox = zeros(nsplits,L,nsplits);
        r2OrisplitVox = zeros(nsplits,L,nsplits);
        rssOrisplitVox = zeros(nsplits,L,nsplits);
        pearsonRoriVox = zeros(nsplits,L,nsplits);
        pearsonRVox = zeros(nsplits,L,nsplits);
        sessBetasSplit = zeros(nsplits,L);
        sessStdBetasSplit = zeros(nsplits,L);
        parfor isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            sessBetasSplit(isplit,:) = mean(roiBetas{roinum}(:,imgTrials>0),2);
            sessStdBetasSplit(isplit,:) = std(roiBetas{roinum}(:,imgTrials>0),0,2);
            for ivox=1:L%nvox(roinum)
                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';
                
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));%includes lowest and highest SF
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                %voxePrfSample*coeff = voxBetas;
                voxCoefSplit(isplit,ivox,:) = voxPrfSample\voxBetas;
                
                voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
                
                %add lowest and highest SF
                voxPrfOriSample(:,end+1:end+2) = voxPrfSample(:,end-2:end-1);

                %add constant predictor
                voxPrfOriSample(:,end+1) = ones;

                %perform regression
                voxOriCoefSplit(isplit,ivox,:) = voxPrfOriSample\voxBetas;%check vox 144 in first ROI
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %regress vignetting predicted timecourse with orientation model
                voxPred = squeeze(voxCoefSplit(isplit,ivox,:))'*voxPrfSample';
                voxPredOriCoefSplit(isplit,ivox,:) = regress(voxPred',voxPrfOriSample);
                
                voxOriPred = squeeze(voxOriCoefSplit(isplit,ivox,:))'*voxPrfOriSample';
                voxOriPredOriCoefSplit(isplit,ivox,:) = regress(voxOriPred',voxPrfOriSample);
                
                %compute within-split residuals
                voxOriResidual = voxBetas' - squeeze(voxOriCoefSplit(isplit,ivox,:))'*voxPrfOriSample';
                voxResidual = voxBetas' - squeeze(voxCoefSplit(isplit,ivox,:))'*voxPrfSample';
                
                
                %regress residuals of vignetting model with full model
                voxResidOriCoefSplit(isplit,ivox,:) = regress(squeeze(voxResidual)',voxPrfOriSample);

                %             voxResidOriCoef{roinum}(isplit,ivox,:) = voxPrfOriSample\squeeze(voxResidualSplit{roinum}(isplit,ivox,:));
                
                %regress residuals of orientation model with orientation model
                voxOriResidOriCoefSplit(isplit,ivox,:) = regress(squeeze(voxOriResidual)',voxPrfOriSample);
                
                
                %r2 within split
                %                 r2{roinum}(isplit,ivox) = 1 - (rssq(voxResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))),3)'./rssq(roiBetas{roinum}(ivox,imgTrials>0),2)).^2;
                r2WithinSplit(isplit,ivox) = rsquared(voxResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                %                 r2ori{roinum}(isplit,ivox) = 1 - (rssq(voxOriResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))),3)'./rssq(roiBetas{roinum}(ivox,imgTrials>0),2)).^2;
                r2oriWithinSplit(isplit,ivox) = rsquared(voxOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                %resid r2 within split
                %residuals of full orientation model:
                oriResidOriResidual = voxBetas' - squeeze(voxOriResidOriCoefSplit(isplit,ivox,:))'*voxPrfOriSample';
                %                 voxOriResidOriR2{roinum}(isplit,ivox) = 1 - (rssq(oriResidOriResidual(1:sum(splitImgTrials(isplit,:))))./rssq(roiBetas{roinum}(ivox,imgTrials>0),2)).^2;
                voxOriResidOriR2Split(isplit,ivox) = rsquared(oriResidOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                
                %prediction of full orientation model:
                oriPredOriResidual = voxBetas' - squeeze(voxOriPredOriCoefSplit(isplit,ivox,:))'*voxPrfOriSample';
                %                 voxOriPredOriR2{roinum}(isplit,ivox) = 1 - (rssq(oriPredOriResidual(1:sum(splitImgTrials(isplit,:))))./rssq(roiBetas{roinum}(ivox,imgTrials>0),2)).^2;
                voxOriPredOriR2Split(isplit,ivox) = rsquared(oriPredOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                
                %residuals of constrained model:
                residOriResidual = voxBetas' - squeeze(voxResidOriCoefSplit(isplit,ivox,:))'*voxPrfOriSample';
                %                 voxResidOriR2{roinum}(isplit,ivox) = 1 - (rssq(residOriResidual(1:sum(splitImgTrials(isplit,:))))./rssq(roiBetas{roinum}(ivox,imgTrials>0),2)).^2;
                voxResidOriR2Split(isplit,ivox) = rsquared(residOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                
                %prediction of constrained model:
                predOriResidual = voxBetas' - squeeze(voxPredOriCoefSplit(isplit,ivox,:))'*voxPrfOriSample';
                %                 voxPredOriR2{roinum}(isplit,ivox) = 1 - (rssq(predOriResidual(1:sum(splitImgTrials(isplit,:))))./rssq(roiBetas{roinum}(ivox,imgTrials>0),2)).^2;
                voxPredOriR2Split(isplit,ivox) = rsquared(predOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));

            end
        end
        sessBetas{roinum} = sessBetasSplit;
        sessStdBetas{roinum} = sessStdBetasSplit;
        
        voxPredOriR2{roinum} = voxPredOriR2Split;
        voxResidOriR2{roinum} = voxResidOriR2Split;
        voxOriPredOriR2{roinum} = voxOriPredOriR2Split;
        voxOriResidOriR2{roinum} = voxOriResidOriR2Split;
        r2ori{roinum} = r2oriWithinSplit;
        r2{roinum} = r2WithinSplit;
        voxResidOriCoef{roinum} = voxResidOriCoefSplit;
        voxOriResidOriCoef{roinum} = voxOriResidOriCoefSplit;
        voxOriCoef{roinum} = voxOriCoefSplit;
        voxCoef{roinum} = voxCoefSplit;
        voxPredOriCoef{roinum} = voxPredOriCoefSplit;
        voxOriPredOriCoef{roinum} = voxOriPredOriCoefSplit;
        
        L = nvox(roinum);
        parfor isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            
            for ivox=1:L%nvox(roinum)
                
                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';
                
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));%includes lowest and highest SF
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                %voxePrfSample*coeff = voxBetas;
%                 voxCoef{roinum}(isplit,ivox,:) = voxPrfSample\voxBetas;
                
                voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
                
                %add lowest and highest SF
                voxPrfOriSample(:,end+1:end+2) = voxPrfSample(:,end-2:end-1);

                %add constant predictor
                voxPrfOriSample(:,end+1) = ones;

                %%
                %cross-split goodness-of-fit (r-squared and correlation)
                %taking coefficients from nextSplit, multiplying by pRF energy in isplit, and comparing to isplit betas
                for nextSplit=1:nsplits
                    voxOriResidualSplit = voxBetas' - squeeze(voxOriCoef{roinum}(nextSplit,ivox,:))'*voxPrfOriSample';
                    voxResidualSplit = voxBetas' - squeeze(voxCoef{roinum}(nextSplit,ivox,:))'*voxPrfSample';
                    
                    %r2 between splits
                    r2splitVox(isplit,ivox,nextSplit) = rsquared(voxResidualSplit(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                    r2OrisplitVox(isplit,ivox,nextSplit) = rsquared(voxOriResidualSplit(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                    
                    %root sum of squares between splits
                    rssSplitVox(isplit,ivox,nextSplit) = sum(voxResidualSplit(1:sum(splitImgTrials(isplit,:))).^2);
                    rssOrisplitVox(isplit,ivox,nextSplit) = sum(voxOriResidualSplit(1:sum(splitImgTrials(isplit,:))).^2);

                    %corr between splits
                    pearsonRoriVox(isplit,ivox,nextSplit) = corr(voxBetas,(squeeze(voxOriCoef{roinum}(nextSplit,ivox,:))'*voxPrfOriSample')');
                    pearsonRVox(isplit,ivox,nextSplit) = corr(voxBetas,(squeeze(voxCoef{roinum}(nextSplit,ivox,:))'*voxPrfSample')');
                end

            end
        end
        r2split{roinum} = r2splitVox;
        r2oriSplit{roinum} = r2OrisplitVox;
        rssOriSplit{roinum} = rssOrisplitVox;
        rssSplit{roinum} = rssSplitVox;
        pearsonRori{roinum} = pearsonRoriVox;
        pearsonR{roinum} = pearsonRVox;
        
        %correlate coefficients between splits
        %with or without constant?
        for ivox=1:L
            voxCoefCorr{roinum}(ivox,:,:) = corr(squeeze(voxCoef{roinum}(:,ivox,1:end-1))');
            voxOriCoefCorr{roinum}(ivox,:,:) = corr(squeeze(voxOriCoef{roinum}(:,ivox,1:end-1))');
            voxCoefCorrWithConst{roinum}(ivox,:,:) = corr(squeeze(voxCoef{roinum}(:,ivox,:))');
            voxOriCoefCorrWithConst{roinum}(ivox,:,:) = corr(squeeze(voxOriCoef{roinum}(:,ivox,:))');
        end
    end
    
    %%
    nsd.voxOriCoefCorrWithConst = voxOriCoefCorrWithConst;
    nsd.voxCoefCorrWithConst = voxCoefCorrWithConst;
    nsd.voxOriCoefCorr = voxOriCoefCorr;    
    nsd.voxCoefCorr = voxCoefCorr;
    % if nsdOrSynth==1b
%     nsd.voxResidual = voxResidual;
%     nsd.voxOriResidual = voxOriResidual;
%     nsd.voxResidualSplit = NaN;%voxResidualSplit; These would be too big
%     nsd.voxOriResidualSplit = NaN;%voxOriResidualSplit; These would be too big
    nsd.r2 = r2;
    nsd.r2ori = r2ori;
    nsd.r2split = r2split;
    nsd.r2oriSplit = r2oriSplit;
    nsd.rssOriSplit = rssOriSplit;
    nsd.rssSplit = rssSplit;
    nsd.pearsonRori = pearsonRori;
    nsd.pearsonR = pearsonR;
    
    % nsd.imgTrials = imgTrials;
    nsd.imgNum = imgNum;
    nsd.splitImgTrials = splitImgTrials;
    % nsd.midTrial = midTrial;
    
    nsd.voxCoef = voxCoef;
    nsd.voxOriCoef = voxOriCoef;
    
    nsd.voxPredOriCoef = voxPredOriCoef;
    nsd.voxOriPredOriCoef = voxOriPredOriCoef;
    
    nsd.voxResidOriCoef = voxResidOriCoef;
    nsd.voxOriResidOriCoef = voxOriResidOriCoef;
    
    nsd.voxPredOriR2 = voxPredOriR2;
    nsd.voxOriPredOriR2 = voxOriPredOriR2;
    nsd.voxResidOriR2 = voxResidOriR2;
    nsd.voxOriResidOriR2 = voxOriResidOriR2;
    
    nsd.sessBetas = sessBetas;
    nsd.sessStdBetas = sessStdBetas;
    
    nsd.roiInd = roiInd;
    

    
    %% SAVE RESULTS
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
    bandpassStr = '';
    if bandpass
        bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
    end
    permStr = '';
    if toPerm
        permStr= '_perm';
    end
    save(fullfile(saveFolder,['regressPrfSplit_' num2str(onlyRep) 'rep_' num2str(nsessions) 'sess_v' num2str(visualRegion) '_sub' num2str(isub) zscoreStr permStr '.mat']), ...
        'nsd', ...%'synth',...
        'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits',...
        'splitImgTrials','sessRepTrials','minNumRep','minSessRep');
    toc
end

end


function r2 = rsquared(Xresid, Xorig)
r2 = 1 - sum((Xresid).^2)/sum((Xorig - mean(Xorig)).^2);
end