% function res = simPopResponse_expand(nrois,toZscore,nimgs)

nrois=1;
toZscore = 0;%0, 1, 2==zeroMean, 3==equalStd
nimgs=100;
% nrois = 4;
subjects = 1:8;%[5:8];

simImgs = 1:nimgs;


saveFolder = fullfile('/','misc','data18','rothzn','nsd','repDrift_expand','/');
% prffolder = fullfile('/','misc','data18','rothzn','nsd','prfsample_expand','/');

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


for isub=subjects
        load(fullfile(saveFolder, ['regressSessCombineROI_sub' num2str(isub) zscoreStr '.mat']),'allRoiPrf',...
        'sessBetas','sessStdBetas','voxConstCoef','voxConstOriCoef',...
        'voxOriCoef','voxCoef');
        keyboard
    for iroi=1:nrois
        betaDiff = diff(sessBetas{iroi},1,1);
        betaDiffCorr(isub,iroi,:,:) = corr(betaDiff');
        load(fullfile(saveFolder,['simPopResp_v' num2str(iroi) '_sub' num2str(isub) zscoreStr '.mat']),...
            'voxSessResp','voxSessRespOri', 'simImgs','nsessions');
        popDiff = diff(squeeze(mean(voxSessRespOri,3)),1,2);%vox X session X image
        popDiffCorr(isub,iroi,:,:) = corr(popDiff);%correlation between drift vectors between adjacent sessions.
keyboard
%         popImgDiff = diff(voxSessRespOri,1,2);%vox X session X image
%         popDiff = squeeze(mean(popImgDiff,3));%mean over images, gives approximation of drift between adjacent sessions. vox X session
        
    end

end
