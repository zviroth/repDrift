%each session is a separate split.

%data saved is used by plotVoxPref.m

function res = regressPrfSplit_session(isub,visualRegions,toZscore,noBackground)
if ieNotDefined('toZscore'), toZscore = 0; end
if ieNotDefined('noBackground'), noBackground = 0; end


%uses file from prfSampleModel_nsd2.m and from
%for synthetic stimuli, averages across all repetitions for each image
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
nsessionsSub = [40 40 32 30 40 32 40 30];


nsessions=nsessionsSub(isub);
nsplits=nsessions;
boxfolder = ['/misc/data18/rothzn/nsd/prfsample' noBackgroundStr '/'];
prffolder = boxfolder;%['~/NSD/prfsample' noBackgroundStr '/'];
saveFolder = ['/misc/data18/rothzn/nsd/repDrift/'];
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


sum([1 2 3]*sum(nsdDesign.basiccnt,2)); %30,000
sum([1 2 3]*nsdDesign.basiccnt,1); %750

%divide trials into nsplits
ntrials = length(nsdDesign.masterordering);
splitImgTrials = zeros(nsplits,ntrials);
trialsPerSession = 10 * 75; %12 runs, 75 trials per run
%     trialsPerSplit = trialsPerSession * 40;%floor or ceil?
itrial = 1;
for isplit=1:nsplits
    if isplit==nsplits
        splitImgTrials(isplit,itrial:end) = ones;
    else
        splitImgTrials(isplit,itrial:itrial+trialsPerSession) = ones;
        itrial = itrial+trialsPerSession;
    end
end



for visualRegion=visualRegions
    ['visual region: ' num2str(visualRegion)]
    
    load(fullfile(boxfolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
        'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
        'roiPrf');
    %if prf sampling was done with the nonlinear CSS prf, then we want to
    %define the weights for the constrained model as a sum of the
    %orientation model across orientations:
    for roinum=1:length(rois)
        prfSampleLev{roinum} = squeeze(sum(prfSampleLevOri{roinum},4));
    end
    
    
    
    
    if bandpass
        for roinum=1:length(rois)
            iroi = rois(roinum);
            prfSampleLevOri{roinum} = prfSampleLevOri{roinum}(:,:,bandMin:bandMax,:);
            prfSampleLev{roinum} = prfSampleLev{roinum}(:,:,bandMin:bandMax);
        end
        numLevels = bandMax-bandMin+1;
    end
    
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
        voxCoefSplit = zeros(nsplits,L,numLevels+1);
        voxOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1);
        voxPredOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1);
        voxOriPredOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1);
        voxResidOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1);
        voxOriResidOriCoefSplit = zeros(nsplits,L,numLevels*numOrientations+1);
        r2WithinSplit = zeros(nsplits,L);
        r2oriWithinSplit = zeros(nsplits,L);
        voxOriResidOriR2Split = zeros(nsplits,L);
        voxOriPredOriR2Split = zeros(nsplits,L);
        voxResidOriR2Split = zeros(nsplits,L);
        voxPredOriR2Split = zeros(nsplits,L);
        r2splitVox = zeros(nsplits,L,nsplits);
        r2OrisplitVox = zeros(nsplits,L,nsplits);
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
                
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                %voxePrfSample*coeff = voxBetas;
                voxCoefSplit(isplit,ivox,:) = voxPrfSample\voxBetas;
                
                voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
                
                %add constant predictor
                voxPrfOriSample(:,end+1) = ones;
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
                
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                %voxePrfSample*coeff = voxBetas;
%                 voxCoef{roinum}(isplit,ivox,:) = voxPrfSample\voxBetas;
                
                voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
                
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
                    
                    %corr between splits
                    pearsonRoriVox(isplit,ivox,nextSplit) = corr(voxBetas,(squeeze(voxOriCoef{roinum}(nextSplit,ivox,:))'*voxPrfOriSample')');
                    pearsonRVox(isplit,ivox,nextSplit) = corr(voxBetas,(squeeze(voxCoef{roinum}(nextSplit,ivox,:))'*voxPrfSample')');
                end

            end
        end
        r2split{roinum} = r2splitVox;
        r2oriSplit{roinum} = r2OrisplitVox;
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
    
    
%     %% SYNTHETIC STIMULI
%     load(fullfile(prffolder,['prfSampleSynth_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
%         'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg','imgScaling',...
%         'roiPrf');
%     if bandpass
%         for roinum=1:length(rois)
%             %             iroi = rois(roinum);
%             prfSampleLevOri{roinum} = prfSampleLevOri{roinum}(:,:,bandMin:bandMax,:);
%             prfSampleLev{roinum} = prfSampleLev{roinum}(:,:,bandMin:bandMax);
%         end
%         numLevels = bandMax-bandMin+1;
%     end
%     
%     allImgs = 105:216;%1:220; ONLY GRATINGS AND SPIRALS
%     
%     betasfilename = fullfile(synthfolder,'betas_nsdsynthetic.nii');
%     betas = niftiread(betasfilename);
%     betas = cast(betas,'double');
%     betas = betas/300;
%     
%     synthDesignFilename = fullfile(nsdfolder, 'nsdsynthetic_expdesign.mat');
%     synthDesign = load(synthDesignFilename);
%     synthDesign.masterordering;
%     
%     [imgTrials, imgNum] = ismember(synthDesign.masterordering, allImgs);%logical array
%     betas=reshape(betas,[],size(betas,4));
%     
%     %     unique(imgNum);%length 113 because it includes 0
%     nimgs=length(allImgs);
%     
%     
%     imgTrialSum = cumsum(imgTrials);
%     totalTrials = max(imgTrialSum);%246
%     
%     
%     
%     clear roiBetas pearsonRori pearsonR voxOriResidual voxResidual
%     for roinum=1:length(rois)
%         iroi = rois(roinum);
%         roiBetas{roinum} = betas(visRoiData==iroi,:);
%     end
%     r2 = cell(length(rois),1);
%     r2ori = cell(length(rois),1);
%     voxOriResidual=cell(length(rois),1);
%     voxResidual=cell(length(rois),1);
%     voxOriCoef = cell(length(rois),1);
%     voxCoef = cell(length(rois),1);
%     for roinum=1:length(rois)
%         nvox(roinum) = size(roiBetas{roinum},1);
%         voxOriResidual{roinum} = zeros(nsplits, nvox(roinum),nimgs);
%         voxResidual{roinum} = zeros(nsplits, nvox(roinum),nimgs);
%         %compute R^2 for both models, within splits, and between splits
%         for isplit=1:nsplits
%             
%             %AVERAGE ACROSS TRIALS OF EACH IMAGE
%             clear imgBetas imgPrfSample imgPrfOriSample
%             for iimg=1:nimgs
%                 imgBetas(:,iimg) = mean(roiBetas{roinum}(:,imgNum==iimg),2);
%             end
%             
%             
%             for ivox=1:nvox(roinum)
%                 voxBetas = imgBetas(ivox,:)';
%                 voxPrfSample = squeeze(prfSampleLev{roinum}(:,ivox,:));
%                 voxPrfSample(:,end+1) = ones;
%                 voxCoef{roinum}(ivox,:) = voxPrfSample\voxBetas;
%                 
%                 voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(:,ivox,:,:));
%                 voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
%                 voxPrfOriSample(:,end+1) = ones;
%                 
%                 voxOriCoef{roinum}(ivox,:) = voxPrfOriSample\voxBetas;%check vox 144 in first ROI
%                 
%                 
%                 voxOriResidual{roinum}(isplit,ivox,:) = voxBetas' - squeeze(nsd.voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
%                 voxResidual{roinum}(isplit,ivox,:) = voxBetas' - squeeze(nsd.voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample';
%                 
%                 pearsonRori{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(nsd.voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample')');
%                 pearsonR{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(nsd.voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample')');
%                 
%             end
%             
%         end
%     end
%     
%     
%     synth.roiBetas = roiBetas;
%     synth.prfSampleLevOri = prfSampleLevOri;
%     synth.prfSampleLev = prfSampleLev;
%     
%     synth.imgTrials = imgTrials;
%     synth.imgNum = imgNum;
%     
%     synth.pearsonR = pearsonR;
%     synth.pearsonRori = pearsonRori;
%     
%     
%     synth.voxOriCoef = voxOriCoef;
%     synth.voxCoef = voxCoef;
%     
%     
%     synth.voxResidual = voxResidual;
%     synth.voxOriResidual = voxOriResidual;

    
    %% SAVE RESULTS
    zscoreStr='';
    if toZscore==1
        zscoreStr = '_zscore';
    elseif toZscore==2
        zscoreStr = '_zeroMean';
    elseif toZscore==3
        zscoreStr = '_equalStd';
    end
    bandpassStr = '';
    if bandpass
        bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
    end
    permStr = '';
    if toPerm
        permStr= '_perm';
    end
    save(fullfile(saveFolder,['regressPrfSplit_session_v' num2str(visualRegion) '_sub' num2str(isub) zscoreStr permStr '.mat']), ...
        'nsd', ...%'synth',...
        'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits');
    toc
end

end


function r2 = rsquared(Xresid, Xorig)
r2 = 1 - sum((Xresid).^2)/sum((Xorig - mean(Xorig)).^2);
end