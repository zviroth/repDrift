% function res = saveBetasVar(isub,visualRegions)
if ieNotDefined('visualRegions'), visualRegions = 1:4; end
subjects=1:8;
% if ieNotDefined('isub'), isub = 1; end

% nsess=30;
subSess = [40 40 32 30 40 32 40 30];
tic

nsdFolder = fullfile('~','misc','data18','rothzn','nsd');
if ~isfolder(nsdFolder)
    nsdFolder = '/misc/data18/rothzn/nsd/';
end

saveFolder = fullfile(nsdFolder,'repDrift_expand','/');
% nsdFolder = fullfile('~','NSD');
for isub=subjects
    betasfolder = fullfile(nsdFolder,['sub' num2str(isub) '_betas_func1pt8mm/']);
    roifolder = betasfolder;
    visualRoisFile = fullfile(roifolder,'prf-visualrois.nii');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
    visRoiData = niftiread(visualRoisFile);
    roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
    visRoiData = visRoiData(:);

    %% LOAD BETAS
clear combinedRoiBetas
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
        for isession=1:subSess(isub)
            betasfilename = fullfile(betasfolder,['betas_session' num2str(isession,'%02.f') '.nii']);
            betas = niftiread(betasfilename);
            betas = cast(betas,'double');
            betas = betas/300;
            %         betas = betas*100;
            betas=reshape(betas,[],size(betas,4));
            for roinum=1:length(rois)
                iroi = rois(roinum);
                %             roiInd{roinum} = find(visRoiData==iroi);
                roiBetas{roinum} = [roiBetas{roinum} betas(visRoiData==iroi,:)];
            end
        end

        clear nvox

        if visualRegion<4
            combinedRoiBetas = [roiBetas{1}; roiBetas{2}];
        else

            combinedRoiBetas = [roiBetas{1}];
        end

        roiBetasVar{isub,visualRegion} = var(combinedRoiBetas,0,2);
        roiBetasVar30{isub,visualRegion} = var(combinedRoiBetas(:,1:30*750),0,2);
        roiBetasVar20{isub,visualRegion} = var(combinedRoiBetas(:,1:20*750),0,2);
    end
end
% roiBetas = combinedRoiBetas;
save(fullfile(saveFolder,['betasVar.mat']),...
    'subjects','visualRegions','subSess','roiBetasVar','roiBetasVar30','roiBetasVar20');

