%uses data saved by regressSessionCombineRoi.m

close all
clear all
saveFigs=0;


% saveFolder = fullfile('~','misc','data18','rothzn','nsd',['repDrift' versionStr]);
% if ~isfolder(saveFolder)
%     saveFolder = '/misc/data18/rothzn/nsd/repDrift/';
% end
saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand','/');
if ~isfolder(saveFolder)
    saveFolder = ['/misc/data18/rothzn/nsd/repDrift_expand/'];
end
figsFolder = fullfile(saveFolder,'figs');

tic
rng(1)
% rois = [1];
iroi=1;
nbins=30;
% nperms = 2;

toZscore = 0;%0, 1, 2==zeroMean, 3==equalStd
% toNormalize = 0;% (r2distant-r2initial)/abs(r2distant+r2initial)
singleSubject=0;

useMedian=1;

subjects = [1:8];
% subjects = [1 2 5 7];


colormapName = 'parula';%'cool';

% percentileThresh = 28;%r2 within percentile
r2thresh = 00;%33;%50;

r2threshStr = '';
if r2thresh>0
    r2threshStr = ['r2thresh' num2str(r2thresh,'%4.2f')];
end

subColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840],...
    [0.4660 0.6740 0.1880]/2+[0.8500 0.3250 0.0980]/2};
% subColor = {[255 179 186]/255, [255 223 186]/255, [255 255 186]/255, [186 255 201]/255, ...
%     [186 225 255]/255, [255 179 255]/255, [255 186 225]/255, [201 186 225]/255};
linewidth=2;
nsubjects=length(subjects);
nrois=4;
nsessions=40;
maxSessions = nsessions;
nsplits = nsessions+1;
errorbarColor = [0.2 0.2 0.2];
surfaceAlpha = 0.1;

distMatrix = toeplitz(0:nsessions-1);%nsplits = nsessions+1. from 0 to 39.
lagMatrix = toeplitz(0:-1:1-nsessions,0:nsessions-1);%from -39 to +39
subSessions = zeros(nsubjects,nsessions);

zscoreStr='';
if toZscore==1
    zscoreStr = '_zscore';
elseif toZscore==2
    zscoreStr = '_zeroMean';
elseif toZscore==3
    zscoreStr = '_equalStd';
end
ifig=1;
autocorrLength = 21;


f=figure(ifig); clf; ifig=ifig+1;
rows=1;
cols=4;
for isub=1:nsubjects
    isub
    subjNum = subjects(isub);
    load(fullfile(saveFolder, ['regressSessCombineROI_sub' num2str(subjNum) zscoreStr '.mat']),'allRoiPrf','roiLevVig','roiLevFull',...
        'roiOri','roiNsdCorr','roiNsdOriCorr','roiNsdOriR2','roiNsdR2',...
        'visRoiData','roiNames','combinedRoiNames','roiInd','prefAnalysis','nsplits',...
        'roiNsdOriR2within','roiNsdR2within',...
        'voxOriCoefCorrWithConst', 'voxCoefCorrWithConst','voxOriCoefCorr','voxCoefCorr',...
        'sessBetas','sessStdBetas','voxConstCoef','voxConstOriCoef',...
        'voxOriCoef','voxCoef');
    
    nsessions = nsplits-1; %maxnumber of intervening sessions. nsplits is <40 for some subjects
    subSessions(isub,1:nsessions) = ones;
    
    %if r2thresh>0 then use only voxels above threshold
    goodVox = allRoiPrf{iroi}.r2>r2thresh;
    r2GoodVox(iroi,isub) = r2thresh;
    numGoodVox(iroi,isub) = sum(goodVox);
    
    for i=1:cols
        switch i
            case 1
                data = sessBetas{iroi};
                titleStr = 'mean beta';
            case 2
                data = sessStdBetas{iroi};
                titleStr = 'std betas';
            case 3
                data = voxConstCoef{iroi}(1:nsessions,:);
                titleStr = 'const. constant coef';
            case 4
                data = voxConstOriCoef{iroi}(1:nsessions,:);
                titleStr = 'full constant coef';
        end
        titles{i} = titleStr;
        figure(1)
        subplot(rows,cols,i)
        crossCorr{isub,i} = corr(data(:,goodVox));
        ecc{isub} = allRoiPrf{iroi}.ecc(goodVox);
        ang{isub} = allRoiPrf{iroi}.ang(goodVox);
        r2{isub} = allRoiPrf{iroi}.r2(goodVox);
        %         histogram(crossCorr(triu(ones(size(crossCorr)),1)==1),'NumBins',nbins,'FaceColor',subColor{isub},'FaceAlpha',0.2, 'Normalization','probability'); hold on;
        histogram(crossCorr{isub,i}(triu(ones(size(crossCorr{isub,i})),1)==1),'NumBins',nbins,'FaceColor',subColor{isub},'FaceAlpha',0.2); hold on;
        if isub==1
            title(titleStr)
        end
        xlabel('correlation (r)');
        ylabel('# voxels');
        axis square
        
    end
    
end

save(fullfile(saveFolder,['driftVoxCorr' zscoreStr '.mat']),'crossCorr','ecc','ang','r2','titles');

set(gcf,'position',[100 100 600 200]);
if saveFigs
    savepdf(f,fullfile(figsFolder,['crosscorr' zscoreStr '.pdf']));
end
%%
eccBorder=0.3;
figure(ifig); clf; ifig=ifig+1;
for isub=1:nsubjects
    for i=1:cols
        subplot(rows,cols,i)
        
        
        corrWithinFov = crossCorr{isub,i}(ecc{isub}<eccBorder,ecc{isub}<eccBorder);
        corrWithinPer = crossCorr{isub,i}(ecc{isub}>=eccBorder,ecc{isub}>=eccBorder);
        corrWithin = [corrWithinFov(triu(ones(size(corrWithinFov)),1)==1); corrWithinPer(triu(ones(size(corrWithinPer)),1)==1)];
        corrFovToPer = crossCorr{isub,i}(ecc{isub}<eccBorder,ecc{isub}>=eccBorder);
        corrPerToFov = crossCorr{isub,i}(ecc{isub}>=eccBorder,ecc{isub}<eccBorder);
        corrBetween = [corrFovToPer(triu(ones(size(corrFovToPer)),1)==1); corrPerToFov(triu(ones(size(corrPerToFov)),1)==1)];
        histogram(corrWithin,'NumBins',nbins,'FaceColor',[0.1 0.1 0.9],'FaceAlpha',0.2, 'Normalization','probability'); hold on;
        histogram(corrBetween,'NumBins',nbins,'FaceColor',[0.9 0.1 0.1],'FaceAlpha',0.2, 'Normalization','probability'); hold on;
        separationInd(isub,i) = (mean(corrWithin) - mean(corrBetween))/(mean(corrWithin) + mean(corrBetween));
        axis square
    end
end
set(gcf,'position',[100 100 1000 400]);

% %% finding the best eccentricity to separate fovea from periphery
% borders = linspace(0.1,8,20);
% tic
% 
% figure(ifig); clf; ifig=ifig+1;
% for i=1:cols
%     for isub=1:nsubjects
%         for iborder=1:length(borders)
%             
%             eccBorder=borders(iborder);
%             corrWithinFov = crossCorr{isub,i}(ecc{isub}<eccBorder,ecc{isub}<eccBorder);
%             corrWithinPer = crossCorr{isub,i}(ecc{isub}>=eccBorder,ecc{isub}>=eccBorder);
%             corrWithin = [corrWithinFov(triu(ones(size(corrWithinFov)),1)==1); corrWithinPer(triu(ones(size(corrWithinPer)),1)==1)];
%             corrFovToPer = crossCorr{isub,i}(ecc{isub}<eccBorder,ecc{isub}>=eccBorder);
%             corrPerToFov = crossCorr{isub,i}(ecc{isub}>=eccBorder,ecc{isub}<eccBorder);
%             corrBetween = [corrFovToPer(triu(ones(size(corrFovToPer)),1)==1); corrPerToFov(triu(ones(size(corrPerToFov)),1)==1)];
%             %             separationInd(i,isub,iborder) = (mean(corrWithin) - mean(corrBetween))/(mean(corrWithin) + mean(corrBetween));
%             separatEccInd(i,isub,iborder) = mean(corrWithin) - mean(corrBetween);
%             
%         end
%     end
%     subplot(rows,cols,i)
%     plot(borders,squeeze(separatEccInd(i,:,:))','.-');
%     axis square
% end
% toc


%% finding the best polar-angle
tic
allPrfParams{1} = ecc;
allPrfParams{2} = ang;
allPrfParams{3} = r2;

borders{1} = linspace(0.1,8,20);
borders{2} = linspace(0,2*pi,20);
borders{3} = linspace(0,80,20);

prfxlabel{1} = 'eccentricity (deg)';
prfxlabel{2} = 'angle (deg)';
prfxlabel{3} = 'R2';


    figure(ifig); clf; ifig=ifig+1;
    rows = length(borders);
for j=1:length(borders)
    prfParam = allPrfParams{j};

    for i=1:cols
        for isub=1:nsubjects
            for iborder=1:length(borders{j})
                
                thisBorder=borders{j}(iborder);
                corrWithinFov = crossCorr{isub,i}(prfParam{isub}<thisBorder,prfParam{isub}<thisBorder);
                corrWithinPer = crossCorr{isub,i}(prfParam{isub}>=eccBorder,prfParam{isub}>=thisBorder);
                corrWithin = [corrWithinFov(triu(ones(size(corrWithinFov)),1)==1); corrWithinPer(triu(ones(size(corrWithinPer)),1)==1)];
                corrFovToPer = crossCorr{isub,i}(prfParam{isub}<thisBorder,prfParam{isub}>=thisBorder);
                corrPerToFov = crossCorr{isub,i}(prfParam{isub}>=thisBorder,prfParam{isub}<thisBorder);
                corrBetween = [corrFovToPer(triu(ones(size(corrFovToPer)),1)==1); corrPerToFov(triu(ones(size(corrPerToFov)),1)==1)];
                %             separationInd(i,isub,iborder) = (mean(corrWithin) - mean(corrBetween))/(mean(corrWithin) + mean(corrBetween));
                separationInd(j,i,isub,iborder) = mean(corrWithin) - mean(corrBetween);
                
            end
        end
        subplot(rows,cols,i+(j-1)*cols)
        plot(borders{j},squeeze(separationInd(j,i,:,:))','.-');
        ylabel('\Deltacorr(within - between)');
        xlabel(prfxlabel{j});
        axis square
    end
end
toc
