function res = simPopResponse(nrois,toZscore,nimgs)
%simPopResponse(4,0,100)
if ieNotDefined('nrois'), nrois = 4; end
if ieNotDefined('toZscore'), toZscore = 0; end%1=GLMdenoise, 2=without GLMdenoise, i.e. betas2.
if ieNotDefined('nimgs'), nimgs = 100; end

tic 

% toZscore = 0;%0, 1, 2==zeroMean, 3==equalStd
% nrois = 4;
subjects = 1:8;%[5:8];

% simImgs = [1:5];
simImgs = 1:nimgs;

% prffolder = '~/Box/NSD/prfsample/';
if ieNotDefined('version'), version = 1; end%1=GLMdenoise, 2=without GLMdenoise, i.e. betas2.
versionStr = '';
if version==2
    versionStr = '2';
end

saveFolder = fullfile('~','misc','data18','rothzn','nsd',['repDrift' versionStr]);
prffolder = fullfile('~','misc','data18','rothzn','nsd','prfsample');
if ~isfolder(saveFolder)
    saveFolder = '/misc/data18/rothzn/nsd/repDrift/';
    prffolder = '/misc/data18/rothzn/nsd/prfsample/';
end
zscoreStr='';
if toZscore==1
    zscoreStr = '_zscore';
elseif toZscore==2
    zscoreStr = '_zeroMean';
elseif toZscore==3
    zscoreStr = '_equalStd';
end


for isub=subjects
    load([saveFolder 'regressSessCombineROI_sub' num2str(isub) zscoreStr '.mat'],'allRoiPrf','roiLevVig','roiLevFull',...
        'visRoiData','roiNames','combinedRoiNames','roiInd','prefAnalysis','nsplits',...
        'sessBetas','sessStdBetas','voxConstCoef','voxConstOriCoef',...
        'voxOriCoef','voxCoef');
    
    nsessions = size(sessBetas{1},1);
    for iroi=1:nrois
        load(fullfile(prffolder,['prfSampleStim_v' num2str(iroi) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
            'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
            'roiPrf');
        if iroi<4%ventral and dorsal ROIs
            prfSampleOri = cat(2,prfSampleLevOri{1},prfSampleLevOri{2});%img, vox, lev
            prfSample = cat(2,prfSampleLev{1},prfSampleLev{2});
        else
            prfSampleOri = prfSampleLevOri{1} ;%img, vox, lev
            prfSample = prfSampleLev{1};
        end
%         for iimg=1:length(simImgs)
           %get each voxel's prf sampled weight for each filter
           %multiply weights by each session's coefficients to get
           %simulated response for each session.
           
           nvox = size(prfSample,2);
           voxSessResp = zeros(nvox,nsessions,length(simImgs));
           voxSessRespOri = zeros(nvox,nsessions,length(simImgs));
           
           for ivox=1:nvox
%                ivox
%                voxPrfSample = squeeze(prfSample(iimg,ivox,:));
%                voxPrfSample(end+1) = ones;
               voxPrfSample = squeeze(prfSample(simImgs,ivox,:));
               voxPrfSample(:,end+1) = ones;
               
               voxPrfSampleOri = squeeze(prfSampleOri(simImgs,ivox,:,:));
               voxPrfSampleOri = reshape(voxPrfSampleOri,[],numLevels*numOrientations);
               voxPrfSampleOri(:,end+1) = ones;
               
               for isess=1:nsessions
                   voxSessCoef = squeeze(voxCoef{iroi}(isess,ivox,:));
                   voxSessOriCoef = squeeze(voxOriCoef{iroi}(isess,ivox,:));
%                     for iimg=1:length(simImgs)
%                        voxSessResp(ivox,isess,iimg) =  voxPrfSample(iimg,:)*voxSessCoef(:)';
%                        voxSessRespOri(ivox,isess,iimg) =  voxPrfSampleOri(iimg,:)*voxSessOriCoef(:)';
%                     end
                    voxSessResp(ivox,isess,:) =  voxPrfSample(:,:)*voxSessCoef(:);
                    voxSessRespOri(ivox,isess,:) =  voxPrfSampleOri(:,:)*voxSessOriCoef(:);
               end
               
           end
        save(fullfile(saveFolder,['simPopResp_v' num2str(iroi) '_sub' num2str(isub) zscoreStr '.mat']),...
            'voxSessResp','voxSessRespOri', 'simImgs','nsessions');
        toc
    end
   
end
