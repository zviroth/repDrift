function res = regressSessionCombineRoi_1rep(isub,numregions,onlyRep,nsessions,toZscore)
if ieNotDefined('onlyRep'), onlyRep = 1; end%only use trials of this repetition. 1,2,3
if ieNotDefined('nsessions'), nsessions = 30; end%will determine the minimum number of images for this repetition in each session
if ieNotDefined('numregions'), numregions = 4; end
if ieNotDefined('toZscore'), toZscore = 0; end

%uses data saved by regressPrfSplit_v2_session.m or regressPrfSplit_session.m


saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand');
if ~isfolder(saveFolder)
    saveFolder = '/misc/data18/rothzn/nsd/repDrift_expand/';
end

global interpSz;
global backgroundSz;
global degPerPix;
global prefAnalysis;
prefAnalysis = 3;
imgScaling = 0.5;
interpSz= 714*imgScaling;
backgroundSz= 1024*imgScaling;
pixPerDeg = imgScaling*714/8.4;%=85
degPerPix = 8.4/(714*imgScaling);

if ieNotDefined('toZscore'), toZscore = 0; end
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
for iregion=1:numregions
    
    visualRegion = iregion;%V1,V2,V3,V4
    
    load(fullfile(saveFolder,['regressPrfSplit_' num2str(onlyRep) 'rep_' num2str(nsessions) 'sess_v' num2str(visualRegion) '_sub' num2str(isub) zscoreStr  '.mat']), ...
        'nsd', ...%'synth',...
        'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits');
    if length(rois)>1
        oldNsd = nsd;
        clear nsd
        
        nsd.sessBetas{1} = [];
        nsd.sessStdBetas{1} = [];
        
        nsd.voxOriCoefCorrWithConst{1} = [];
        nsd.voxCoefCorrWithConst{1} = [];
        nsd.voxOriCoefCorr{1} = [];
        nsd.voxCoefCorr{1} = [];
    
        nsd.r2{1} = [];
        nsd.r2ori{1} = [];
        nsd.r2split{1} = [];
        nsd.r2oriSplit{1} = [];
        nsd.rssSplit{1} = [];
        nsd.rssOriSplit{1} = [];
        nsd.voxCoef{1} = [];
        nsd.voxOriCoef{1} = [];
        nsd.voxPredOriCoef{1} = [];
        nsd.voxOriPredOriCoef{1} = [];
        nsd.voxResidOriCoef{1} = [];
        nsd.voxOriResidOriCoef{1} = [];
        
        nsd.voxPredOriR2{1} = [];
        nsd.voxOriPredOriR2{1} = [];
        nsd.voxResidOriR2{1} = [];
        nsd.voxOriResidOriR2{1} = [];
        
        nsd.pearsonRori{1} = [];
        nsd.pearsonR{1} = [];
        
        nsd.roiInd{1} = [];

        
        for iroi=1:length(rois)%rois=2
            
            nsd.sessBetas{1} = cat(2,nsd.sessBetas{1}, oldNsd.sessBetas{iroi});
            nsd.sessStdBetas{1} = cat(2,nsd.sessStdBetas{1}, oldNsd.sessStdBetas{iroi});
            
            nsd.roiInd{1} = cat(1,nsd.roiInd{1}, oldNsd.roiInd{iroi});
            nsd.pearsonRori{1} = cat(2,nsd.pearsonRori{1},oldNsd.pearsonRori{iroi});
            nsd.pearsonR{1} = cat(2,nsd.pearsonR{1},oldNsd.pearsonR{iroi});
            nsd.r2split{1} = cat(2,nsd.r2split{1},oldNsd.r2split{iroi});
            nsd.r2oriSplit{1} = cat(2,nsd.r2oriSplit{1},oldNsd.r2oriSplit{iroi});
            nsd.rssSplit{1} = cat(2,nsd.rssSplit{1},oldNsd.rssSplit{iroi});
            nsd.rssOriSplit{1} = cat(2,nsd.rssOriSplit{1},oldNsd.rssOriSplit{iroi});
            nsd.voxCoef{1} = cat(2,nsd.voxCoef{1},oldNsd.voxCoef{iroi});
            nsd.voxOriCoef{1} = cat(2,nsd.voxOriCoef{1},oldNsd.voxOriCoef{iroi});
            nsd.voxPredOriCoef{1} = cat(2,nsd.voxPredOriCoef{1},oldNsd.voxPredOriCoef{iroi});
            nsd.voxOriPredOriCoef{1} = cat(2,nsd.voxOriPredOriCoef{1},oldNsd.voxOriPredOriCoef{iroi});
            nsd.voxResidOriCoef{1} = cat(2,nsd.voxResidOriCoef{1},oldNsd.voxResidOriCoef{iroi});
            nsd.voxOriResidOriCoef{1} = cat(2,nsd.voxOriResidOriCoef{1},oldNsd.voxOriResidOriCoef{iroi});
            
            nsd.r2{1} = cat(2,nsd.r2{1},oldNsd.r2{iroi});
            nsd.r2ori{1} = cat(2,nsd.r2ori{1},oldNsd.r2ori{iroi});
            
            nsd.voxOriCoefCorrWithConst{1} = cat(1,nsd.voxOriCoefCorrWithConst{1},oldNsd.voxOriCoefCorrWithConst{iroi});
            nsd.voxCoefCorrWithConst{1} = cat(1,nsd.voxCoefCorrWithConst{1},oldNsd.voxCoefCorrWithConst{iroi});
            nsd.voxOriCoefCorr{1} = cat(1,nsd.voxOriCoefCorr{1},oldNsd.voxOriCoefCorr{iroi});
            nsd.voxCoefCorr{1} = cat(1,nsd.voxCoefCorr{1},oldNsd.voxCoefCorr{iroi});
        
            nsd.voxPredOriR2{1} = cat(2,nsd.voxPredOriR2{1},oldNsd.voxPredOriR2{iroi});
            nsd.voxOriPredOriR2{1} = cat(2,nsd.voxOriPredOriR2{1},oldNsd.voxOriPredOriR2{iroi});
            nsd.voxResidOriR2{1} = cat(2,nsd.voxResidOriR2{1},oldNsd.voxResidOriR2{iroi});
            nsd.voxOriResidOriR2{1} = cat(2,nsd.voxOriResidOriR2{1},oldNsd.voxOriResidOriR2{iroi});
            
            
            
%             synth.voxResidual{1} = cat(2,synth.voxResidual{1},oldSynth.voxResidual{iroi});
%             synth.voxOriResidual{1} = cat(2,synth.voxOriResidual{1},oldSynth.voxOriResidual{iroi});
%             synth.pearsonRori{1} = cat(2,synth.pearsonRori{1},oldSynth.pearsonRori{iroi});
%             synth.pearsonR{1} = cat(2,synth.pearsonR{1},oldSynth.pearsonR{iroi});
%             synth.voxCoef{1} = cat(1,synth.voxCoef{1},oldSynth.voxCoef{iroi});
%             synth.voxOriCoef{1} = cat(1,synth.voxOriCoef{1},oldSynth.voxOriCoef{iroi});
            
        end
        oldPrf = roiPrf; clear roiPrf;
        roiPrf{1}.ecc=[];
        roiPrf{1}.ang=[];
        roiPrf{1}.sz=[];
        roiPrf{1}.r2=[];
        roiPrf{1}.x=[];
        roiPrf{1}.y=[];
        for iroi=1:length(rois)
            roiPrf{1}.ecc = cat(1,roiPrf{1}.ecc,oldPrf{iroi}.ecc);
            roiPrf{1}.ang = cat(1,roiPrf{1}.ang,oldPrf{iroi}.ang);
            roiPrf{1}.sz = cat(1,roiPrf{1}.sz,oldPrf{iroi}.sz);
            roiPrf{1}.r2 = cat(1,roiPrf{1}.r2,oldPrf{iroi}.r2);
            roiPrf{1}.x = cat(1,roiPrf{1}.x,oldPrf{iroi}.x);
            roiPrf{1}.y = cat(1,roiPrf{1}.y,oldPrf{iroi}.y);
        end
        rois = 1;
    end
    iroi=1;
    
    %% AVERAGE SPLITS
    
    nsd.voxOriCoefCorrWithConst{1}(nsplits+1,:,:) = mean(nsd.voxOriCoefCorrWithConst{1},1);
    nsd.voxCoefCorrWithConst{1}(nsplits+1,:,:) = mean(nsd.voxCoefCorrWithConst{1},1);
    nsd.voxOriCoefCorr{1}(nsplits+1,:,:) = mean(nsd.voxOriCoefCorr{1},1);
    nsd.voxCoefCorr{1}(nsplits+1,:,:) = mean(nsd.voxCoefCorr{1},1);
            
    nsd.voxCoef{1}(nsplits+1,:,:) = mean(nsd.voxCoef{1},1);
    nsd.voxOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriCoef{1},1);
    nsd.voxPredOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxPredOriCoef{1},1);
    nsd.voxOriPredOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriPredOriCoef{1},1);
    nsd.voxResidOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxResidOriCoef{1},1);
    nsd.voxOriResidOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriResidOriCoef{1},1);
    
    nsd.pearsonRori{1}(nsplits+1,:,:) = mean(nsd.pearsonRori{1},1);
    nsd.pearsonR{1}(nsplits+1,:,:) = mean(nsd.pearsonR{1},1);
    nsd.r2{1}(nsplits+1,:) = mean(nsd.r2{1},1);
    nsd.r2ori{1}(nsplits+1,:) = mean(nsd.r2ori{1},1);
    nsd.r2split{1}(nsplits+1,:,:) = mean(nsd.r2split{1},1);
    nsd.r2oriSplit{1}(nsplits+1,:,:) = mean(nsd.r2oriSplit{1},1);
    nsd.rssSplit{1}(nsplits+1,:,:) = mean(nsd.rssSplit{1},1);
    nsd.rssOriSplit{1}(nsplits+1,:,:) = mean(nsd.rssOriSplit{1},1);
%     synth.voxResidual{1}(nsplits+1,:,:) = mean(synth.voxResidual{1},1);
%     synth.voxOriResidual{1}(nsplits+1,:,:) = mean(synth.voxOriResidual{1},1);
%     synth.pearsonRori{1}(nsplits+1,:) = mean(synth.pearsonRori{1},1);
%     synth.pearsonR{1}(nsplits+1,:) = mean(synth.pearsonR{1},1);
    
    
    nsplits = nsplits+1;
    
   
    %% COMPUTE PREFERRED ORIENTATION - CIRCULAR CENTER OF MASS
    
    clear fullPrefOri residPrefOri residOriPrefOri predCoefOri predOriCoefOri oriDeviation vertDeviation cardDeviation
    clear fullOriModul fullPrefAmp fullAntiAmp
   
    %save model coefficients
    voxOriCoef{iregion} = nsd.voxOriCoef{iroi};
    voxCoef{iregion} = nsd.voxCoef{iroi};
    
    %save preferred orientation and level for this ROI
    allRoiPrf{iregion} = roiPrf{iroi};
%     roiLevVig{iregion} = vigPrefLevel{iroi};
%     roiLevFull{iregion} = fullPrefLevel{iroi};
    %     roiOri{iregion} = prefOri{iroi};
%     roiOri{iregion} = fullPrefOri{iroi};
%     residOri{iregion} = residPrefOri{iroi};
%     residOriOri{iregion} = residOriPrefOri{iroi};
%     predOri{iregion} = predCoefOri{iroi};
%     predOriOri{iregion} = predOriCoefOri{iroi};
%     oriModulation{iregion} = fullOriModul{iroi};
%     oriPrefAmp{iregion} = fullPrefAmp{iroi};
%     oriAntiAmp{iregion} = fullAntiAmp{iroi};
    
    roiInd{iregion} = nsd.roiInd{iroi};
    roiNsdCorr{iregion} = nsd.pearsonR{iroi};
    roiNsdOriCorr{iregion} = nsd.pearsonRori{iroi};
    roiNsdOriR2{iregion} = nsd.r2oriSplit{iroi};
    roiNsdR2{iregion} = nsd.r2split{iroi};
    roiNsdOriRss{iregion} = nsd.rssOriSplit{iroi};
    roiNsdRss{iregion} = nsd.rssSplit{iroi};
    roiNsdR2within{iregion} = nsd.r2{iroi};
    roiNsdOriR2within{iregion} = nsd.r2ori{iroi};
    
    %constant conefficent
    voxConstCoef{iregion} = squeeze(nsd.voxCoef{iroi}(:,:,end));
    voxConstOriCoef{iregion} = squeeze(nsd.voxOriCoef{iroi}(:,:,end));
    
    voxOriCoefCorrWithConst{iregion} = nsd.voxOriCoefCorrWithConst{iroi};
    voxCoefCorrWithConst{iregion} = nsd.voxCoefCorrWithConst{iroi};
    voxOriCoefCorr{iregion} = nsd.voxOriCoefCorr{iroi};
    voxCoefCorr{iregion} = nsd.voxCoefCorr{iroi};
    
    sessBetas{iregion} = nsd.sessBetas{iroi};
    sessStdBetas{iregion} = nsd.sessStdBetas{iroi};
    
    roiNsdOriPredOriR2{iregion} = nsd.voxOriPredOriR2{iroi};
    roiNsdOriResidOriR2{iregion} = nsd.voxOriResidOriR2{iroi};
    roiNsdPredOriR2{iregion} = nsd.voxPredOriR2{iroi};
    roiNsdResidOriR2{iregion} = nsd.voxResidOriR2{iroi};
    
%     roiSynthCorr{iregion} = synth.pearsonR{iroi};
%     roiSynthOriCorr{iregion} = synth.pearsonRori{iroi};
%     roiSynthOri{iregion} = synthFullPrefOri{iroi};
%     roiSynthLevVig{iregion} = synthVigPrefLevel{iroi};
%     roiSynthLevFull{iregion} = synthFullPrefLevel{iroi};
%     roiOriDeviation{iregion} = oriDeviation{iroi};
%     roiVertDeviation{iregion} = vertDeviation{iroi};
%     roiCardDeviation{iregion} = cardDeviation{iroi};

end




%save all ROIs to create overlay
roifolder = ['~/misc/data18/rothzn/nsd/sub' num2str(isub) '_betas_func1pt8mm/'];
if ~isfolder(roifolder)
    roifolder = ['/misc/data18/rothzn/nsd/sub' num2str(isub) '_betas_func1pt8mm/'];
end
visualRoisFile = fullfile(roifolder,'prf-visualrois.nii');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
combinedRoiNames = {'V1','V2','V3','hV4'};
%If we analyzed all 4 ROIs - save the data

% if numregions==4
    save(fullfile(saveFolder, ['regressSessCombineROI_' num2str(onlyRep) 'rep_' num2str(nsessions) 'sess_sub' num2str(isub) zscoreStr '.mat']),'allRoiPrf',...
        'roiNsdCorr','roiNsdOriCorr','roiNsdOriR2','roiNsdR2',...
        'roiNsdOriRss','roiNsdRss','numregions',...
        'roiNsdResidOriR2','roiNsdOriResidOriR2','roiNsdPredOriR2','roiNsdOriPredOriR2',...
        'visRoiData','roiNames','combinedRoiNames','roiInd','prefAnalysis','nsplits',...
        'roiNsdOriR2within','roiNsdR2within',...
        'voxOriCoefCorrWithConst', 'voxCoefCorrWithConst','voxOriCoefCorr','voxCoefCorr',...
        'sessBetas','sessStdBetas','voxConstCoef','voxConstOriCoef',...
        'voxOriCoef','voxCoef');
% end


end
