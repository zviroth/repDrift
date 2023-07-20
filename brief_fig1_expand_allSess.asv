close all
clear all
%uses data saved by savePerms_expand.m
tic

saveFigs=0;

addColorbars = 0;
fixedFirst=0;
toZscore=0;%0=none, 1=zscore, 2=zero mean, 3=normalized std, 4=zero ROI mean
r2thresh = 0;%33;%50;

singleSubject=1;

coreyVersion = 2;
addScatter = 0;

% version=2;
% if ieNotDefined('version'), version = 1; end%1=GLMdenoise, 2=without GLMdenoise, i.e. betas2.
subjects = [1:8];

% subjects = [1 2 5 7];
nperms=1000;
histBins = 20;
scatterSize=5;
scatterFill = 'filled';

fixedFirstStr='';
if fixedFirst
    fixedFirstStr = '_fixedFirst_';
end


toNormalize = 0;

figRoi=1;
linewidthWide=2;
linewidthNarrow = 1;
errorbarColor = [0.2 0.2 0.2];
surfaceAlpha = 0.1;
subColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840],...
    [0.4660 0.6740 0.1880]/2+[0.8500 0.3250 0.0980]/2};
graySubColor = 0.7*ones(1,3);
%histogram color
faceColor = graySubColor;
edgeColor = 'none';

rois = [1];
versionStr = '';
if version==2
    versionStr = '2';
end
saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand');
figsFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand_figs');
if ~isfolder(saveFolder)
    saveFolder = ['/misc/data18/rothzn/nsd/repDrift_expand/'];
    figsFolder = ['/misc/data18/rothzn/nsd/repDrift_expand_figs/'];
end

colormapName = 'parula';%'cool';
% colormapName = 'cool';
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
normalizeStr = '';
if toNormalize
    normalizeStr = '_normalized';
end

r2threshStr = '';
if r2thresh>0
    r2threshStr = ['r2thresh' num2str(r2thresh,'%4.0f')];
end
load(fullfile(saveFolder, ['perm_allSess' fixedFirstStr num2str(nperms) zscoreStr normalizeStr r2threshStr '.mat']),...
    'toNormalize','toZscore', 'useMedian','r2thresh','nrois','rois', ...
    'permOrders', 'subSessions',...
    'r2split','r2oriSplit','pearsonRori','pearsonR',...
    'r2Dist','r2OriDist','pearsonDist','pearsonOriDist',...
    'r2DistPerm','r2OriDistPerm','pearsonDistPerm','pearsonOriDistPerm',....
    'r2DistSess','r2OriDistSess','pearsonDistSess','pearsonOriDistSess',...
    'r2DistSessPerm','r2OriDistSessPerm','pearsonDistSessPerm','pearsonOriDistSessPerm',...
    'meanAutocorrBetas', 'autocorrMeanBetas','meanAutocorrStdBetas','autocorrMeanStdBetas',...
    'meanAutocorrConstant','meanAutocorrConstantOri','autocorrMeanConstant','autocorrMeanConstantOri',...
    'autocorrMeanCoef','autocorrMeanCoefOri','meanAutocorrCoef','meanAutocorrCoefOri',...
    'meanAutocorrBetasPerm','autocorrMeanBetasPerm','meanBetas',...
    'meanAutocorrStdBetasPerm','autocorrMeanStdBetasPerm',...
    'meanAutocorrConstantPerm','meanAutocorrConstantOriPerm',...
    'autocorrMeanConstantPerm','autocorrMeanConstantOriPerm',...
    'autocorrMeanCoefPerm','autocorrMeanCoefOriPerm',...
    'meanAutocorrCoefPerm','meanAutocorrCoefOriPerm',...
    'subRoiPrf','numGoodVox',...
    'r2perm','r2oriPerm','pearsonPerm','pearsonOriPerm');
nsubjects = length(subjects);
ifig=0;

%% use data for selected subjects
subSessions = subSessions(subjects,:);
numSessions = sum(subSessions');

%% distance matrix
% sess = 1:minSessions;
% sessDiff = sess - sess';
nsessions = 40;
distMatrix = toeplitz(0:nsessions-1);
sessDistVec = distMatrix(1:nsessions,1:nsessions);
%%
ifig=ifig+1;

iroi=figRoi;%nrois
f=figure(ifig); ifig=ifig+1;  isubplot=0;
rows=4;
cols=3;


for i=1:2
    switch i
        case 1
            similarityData = r2oriSplit{iroi};
            distSessData = r2OriDistSess;
            distData = r2OriDist{iroi};
            titleStr = 'cvR^2';
            ylabelStr = 'cvR^2';
            distSessDataPerm = r2OriDistSessPerm;
            distDataPerm = r2OriDistPerm{iroi};
            similarityPerm = r2oriPerm{iroi};
        case 2
            similarityData = pearsonRori{iroi};
            distSessData = pearsonOriDistSess;
            distData = pearsonOriDist{iroi};
            titleStr = 'Pearson''s r';
            ylabelStr = 'Pearson''s r';
            distSessDataPerm = pearsonOriDistSessPerm;
            distDataPerm = pearsonOriDistPerm{iroi};
            similarityPerm = pearsonOriPerm{iroi};

    end
    %keep data for chosen subjects only
    similarityData = similarityData(subjects,:,:);
    distData = distData(subjects,:);
    distDataPerm = distDataPerm(subjects,:,:);

    for isub=1:length(subjects)
        %distance matrix defined by session number
        sessSubDist(isub,:,:) =  sessDistVec;
    end

    %%
    subSessDistVec = sessSubDist;

    %keep data for chosen subjects only
    similarityData = similarityData(subjects,:,:);

    %similarity as function of distance, averaged over all possible
    %train-test pairs
    subplot(rows,cols,3*cols+i); axis square; hold all

    %compute correlation per subject and then average
    % scatter(subDistMatrix(subDistMatrix(:)>0),r2oriSplit{iroi}(subDistMatrix(:)>0))
    for isub=1:length(subjects)
        distVec = subSessDistVec(isub,1:numSessions(isub),1:numSessions(isub));
        %correlate this factor with goodness-of-fit
        subSimilarityData = similarityData(isub,1:numSessions(isub),1:numSessions(isub));
        [subDistCorr(isub,i) pSubDistCorr(isub,i)] = corr(subSimilarityData(distVec>0),distVec(distVec>0));

        permDistMat = squeeze(similarityPerm(isub,:,:,:));
        for iperm=1:nperms
            permDistVec = permDistMat(iperm,1:numSessions(isub),1:numSessions(isub));
            [subDistCorrPerm(i,isub,iperm), pSubDistCorrPerm(i,isub,iperm)] = corr(distVec(distVec>0),permDistVec(distVec>0));
        end

        %         pSubPermDistCorr(isub,i) = 1-sum(subDistCorrPerm(i,isub,:)>=subDistCorr(isub,i))/nperms;
        pSubPermDistCorr(isub,i) = sum(subDistCorrPerm(i,isub,:)<=subDistCorr(isub,i))/nperms;

        %PLOT
        for idist=1:numSessions(isub)-1
            subMeanDistData{isub}(idist) = mean(subSimilarityData(distVec(:)==idist));
        end
        plot(subMeanDistData{isub},'linewidth', linewidthNarrow,'color',subColor{isub}); hold all
    end
    axis square
    %     xlabel('distance (sessions)');
    xlabel('\Delta session');
    ylabel(ylabelStr);
    xlim([1 size(distData,2)-1]);
    xticks([1 size(distData,2)-1]);
    xticklabels([1 size(distData,2)-1]);
    box on
    keyboard
end
if saveFigs
        savepdf(f,fullfile(figsFolder,['brief_fig1_allSess' fixedFirstStr zscoreStr '.pdf']));
end

%%
%p-values fpr single subjects
pSubPermDistCorr
