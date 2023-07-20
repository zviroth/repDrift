close all
clear all
tic

saveFigs=0;
imgFormat = '.jpg';

% coreyVersion = 1;
toZscore=0;%0=none, 1=zscore, 2=zero mean, 3=normalized std

zscoreYlim(1,:) = [-0.22 -0.12];
zscoreYlim(2,:) = [-0.15 -0.05];

nperms=1000;
r2thresh = 00;%33;%50;
fixedFirst=0;
subjects = [1:8];


% saveFolder = fullfile('~','misc','data18','rothzn','nsd',['repDrift' versionStr]);
% if ~isfolder(saveFolder)
%     saveFolder = '/misc/data18/rothzn/nsd/repDrift/';
% end
saveFolder =  '/misc/data18/rothzn/nsd/repDrift_expand/';
figsFolder = '/misc/data18/rothzn/nsd/repDrift_expand_figs/';


fixedFirstStr='';
if fixedFirst
    fixedFirstStr = '_fixedFirst_';
end


% subjects = [1 2 5 7];

histBins = 20;

scatterSize=5;
scatterFill = 'filled';
toNormalize = 0;
singleSubject=1;
figRoi=1;
linewidthWide=2;
linewidthNarrow = 1;
errorbarColor = [0.2 0.2 0.2];
surfaceAlpha = 0.1;
subColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840],...
    [0.4660 0.6740 0.1880]/2+[0.8500 0.3250 0.0980]/2};

colormapName = 'parula';%'cool';

zscoreStr='';
if toZscore==1
    zscoreStr = '_zscore';
elseif toZscore==2
    zscoreStr = '_zeroMean';
elseif toZscore==3
    zscoreStr = '_equalStd';
end
normalizeStr = '';
if toNormalize
    normalizeStr = '_normalized';
end

r2threshStr = '';
if r2thresh>0
    r2threshStr = ['r2thresh' num2str(r2thresh,'%4.0f')];
end
load(fullfile(saveFolder, ['perm' fixedFirstStr num2str(nperms) zscoreStr normalizeStr r2threshStr '.mat']),...
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
    'meanAutocorrBetasPerm','autocorrMeanBetasPerm',...
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
minSessions = min(sum(subSessions,2));

%% distance matrix for correlating R2 with inter-session distance
behaviorFolder = '/misc/data18/rothzn/nsd/behavior/';
load(fullfile(behaviorFolder,'subSessBehavior.mat'),'subjects','nsessions','ntrials',...
    'sessTime','sessOld','sessCurrentOld','sessRespond','sessRespondOld',...
    'sessCorrect','sessIncorrect','sessHits','sessMisses','sessFA','sessCR',...
    'sessEasyRespondOld','sessHardRespondOld','sessNovelRespondOld',...
    'sessRT','sessCorrectRT','sessIncorrectRT','sessRespondOldRT','sessRespondNewRT',...
    'sessEasyRespondOldRT','sessEasyRespondNewRT','sessHardRespondOldRT',...
    'sessHardRespondNewRT','sessNovelRespondOldRT','sessNovelRespondNewRT','sessBias','sessDprime',...
    'hitRate','faRate','sessCorrect','sessHour','adjHitRate');
sessCorrect = sessCorrect./750;
nsessions = 40;
distMatrix = toeplitz(0:nsessions-1);
% subDistMatrix = repmat(distMatrix,1,1,length(subjects));
% subDistMatrix = permute(subDistMatrix,[3 1 2]);
% sessDistVec = subDistMatrix(:);
sessDistVec = distMatrix(1:minSessions,1:minSessions);



for isub=1:length(subjects)
    %distance matrix defined by session number
    sessSubDist(isub,:,:) =  sessDistVec;
    distStr{1} = 'session';
    distUnit{1} = 'session';
    %distance matrix defined by session time/date
    timeSubDist(isub,:,:) = abs(sessTime(isub,1:minSessions) - sessTime(isub,1:minSessions)');
    distStr{2} = 'date';
    distUnit{2} = 'day';
    %distance matrix defined by session hour
    hourSubDist(isub,:,:) = abs(sessHour(isub,1:minSessions) - sessHour(isub,1:minSessions)');
    distStr{3} = 'time-of-day';
    distUnit{3} = 'hour';
    %distance matrix defined by decision criterion (=bias)
    biasSubDist(isub,:,:) = abs(sessBias(isub,1:minSessions) - sessBias(isub,1:minSessions)');
    distStr{4} = 'bias';
    distUnit{4} = 'criterion';
    %distance matrix defined by d prime (a measure of performance)
    dprimeSubDist(isub,:,:) = abs(sessDprime(isub,1:minSessions) - sessDprime(isub,1:minSessions)');
    distStr{5} = 'd-prime';
    distUnit{5} = 'd''';
        %distance matrix defined by adjusted hit rate (hit rate minus false
        %alarms)
    adjhitSubDist(isub,:,:) = abs(adjHitRate(isub,1:minSessions) - adjHitRate(isub,1:minSessions)');
    distStr{6} = 'adj. hit rate';
    distUnit{6} = 'prop. hits';
    %distance matrix defined by correctness
    correctSubDist(isub,:,:) = abs(sessCorrect(isub,1:minSessions) - sessCorrect(isub,1:minSessions)');
    distStr{7} = 'correctness';
    distUnit{7} = 'prop. correct';
end

%%
ifig=ifig+1;

iroi=figRoi;%nrois

f=figure(ifig); ifig=ifig+1;  isubplot=0;
rows=5;
cols=1*length(distStr);
for i=1:1%2
    switch i
        %         case 1
        %             similarityData = r2split{iroi};
        %             distSessData = r2DistSess;
        %             distData = r2Dist{iroi};
        %             titleStr = 'constrained R^2';
        %             ylabelStr = 'R^2';
        %             distSessDataPerm = r2DistSessPerm;
        %             distDataPerm = r2DistPerm{iroi};%isub,iperm,isession
        %             similarityPerm = r2perm{iroi};
        case 1
            similarityData = r2oriSplit{iroi};
            titleStr = 'cvR^2';
            ylabelStr = 'cvR^2';
            similarityPerm = r2oriPerm{iroi};
            %         case 3
            %             similarityData = pearsonR{iroi};
            %             distSessData = pearsonDistSess;
            %             distData = pearsonDist{iroi};
            %             titleStr = 'constrained r';
            %             ylabelStr = 'correlation (r)';
            %             distSessDataPerm = pearsonDistSessPerm;
            %             distDataPerm = pearsonDistPerm{iroi};
            %             similarityPerm = pearsonPerm{iroi};
        case 2
            similarityData = pearsonRori{iroi};
            titleStr = 'correlation';
            ylabelStr = 'correlation (r)';
            similarityPerm = pearsonOriPerm{iroi};
    end
    %keep data for chosen subjects only
    similarityData = similarityData(subjects,1:minSessions,1:minSessions);
    
    %matrix: trained on each session and tested on every session
    isubplot=isubplot+1; subplot(rows,cols,(i-1)*length(distStr)+1);
    temp = squeeze(mean(similarityData));

    imagesc(temp); axis square%mean across subjects
    title(titleStr);
    ylabel('trained session');
    xlabel('tested session');
    colormap(colormapName);
    
    for idistType = 1:length(distStr)
        switch idistType
            case 1
                subSessDistVec = sessSubDist;
            case 2
                subSessDistVec = timeSubDist;
                subData = sessTime(:,1:minSessions);
            case 3
                subSessDistVec = hourSubDist;
                subData = sessHour(:,1:minSessions);
            case 4
                subSessDistVec = biasSubDist;
                subData = sessBias(:,1:minSessions);
            case 5
                subSessDistVec = dprimeSubDist;
                subData = sessDprime(:,1:minSessions);
            case 6
                subSessDistVec = adjhitSubDist;
                subData = adjHitRate(:,1:minSessions);
            case 7
                subSessDistVec = correctSubDist;
                subData = sessCorrect(:,1:minSessions);
        end
        %first, just plot this measure as a function of session number
        if idistType>1
            subplot(rows,cols,(i-1)*length(distStr)+idistType);
            for isub=1:length(subjects)
                plot(subData(isub,:),'color',subColor{isub},'linewidth',linewidthNarrow); hold on
            end
            axis square
            title(distStr{idistType});
            ylabel(distUnit{idistType});
            xticks([1 30]);
            xlabel('session');
        end
        %similarity as function of distance, averaged over all possible
        %train-test pairs
        subplot(rows,cols,(i-1)*length(distStr)+cols+idistType); axis square; hold all
        
        %compute correlation per subject and then average
        % scatter(subDistMatrix(subDistMatrix(:)>0),r2oriSplit{iroi}(subDistMatrix(:)>0))
        for isub=1:length(subjects)
            distVec = subSessDistVec(isub,1:minSessions,1:minSessions);
            %distance in #sessions
            sessDistVec = squeeze(sessSubDist(isub,:,:));
            
            %correlate this factor with goodness-of-fit
            subSimilarityData = similarityData(isub,:,:);
            [subDistCorr(isub,i) pSubDistCorr(isub,i)] = corr(subSimilarityData(sessDistVec>0),distVec(sessDistVec>0));
            
            %partially correlate after accounting for this factor
%             temp = subSessDistVec(isub,:,:);
            [subPartCorrSess(isub,idistType,i) pSubPartCorrSess(isub,idistType,i)] = partialcorr(subSimilarityData(sessDistVec>0),sessDistVec(sessDistVec>0),distVec(sessDistVec>0));
            
            %partially correlate this factor after accounting for session
            [subPartCorrBehav(isub,idistType,i) pSubPartCorrBehav(isub,idistType,i)] = partialcorr(subSimilarityData(sessDistVec>0),distVec(sessDistVec>0),sessDistVec(sessDistVec>0));
            
%             corr(sessDistVec(sessDistVec>0),subData(sessDistVec>0));
            
            permDistMat = squeeze(similarityPerm(isub,:,:,:));
            for iperm=1:nperms
                permDistVec = permDistMat(iperm,1:minSessions,1:minSessions);
                [subDistCorrPerm(i,isub,iperm), pSubDistCorrPerm(i,isub,iperm)] = corr(sessDistVec(sessDistVec>0),permDistVec(sessDistVec>0));
            end
            %plot
            scatter1 = scatter(distVec(sessDistVec(:)>0),subSimilarityData(sessDistVec>0),scatterSize,repmat(subColor{isub},sum(sessDistVec(:)>0),1),scatterFill);
            scatter1.MarkerFaceAlpha = surfaceAlpha;
            scatter1.MarkerEdgeAlpha = surfaceAlpha;
            hold on
            %average this factor according to session distance
            for idist=1:minSessions-1
                subMeanDistData(isub,idist) = mean(distVec(sessDistVec(:)==idist));
            end
            plot(subMeanDistData(isub,:),'linewidth', linewidthNarrow,'color',subColor{isub}); hold all
            pSubPermDistCorr(isub,i) = 1-sum(subDistCorrPerm(i,isub,:)>=subDistCorr(isub,i))/nperms;
        end
        distCorr(i) = mean(subDistCorr(:,i),1);%mean over subjects
        distCorrPerm(i,:) = mean(subDistCorrPerm(i,:,:),2);%mean over subjects
        pPermDistCorr(i) = 1-sum(distCorrPerm(i,:)>=distCorr(i))/nperms;
        plot(mean(subMeanDistData,1),'linewidth', linewidthWide,'color','k'); hold all
        
        %     title(['r=' num2str(distCorr(i),'%4.2f') ' p=' num2str(pDistCorr(i),'%4.3f')]);
        if ~saveFigs
            title(['r=' num2str(distCorr(i),'%4.2f') ' p=' num2str(pPermDistCorr(i),'%4.3f')]);
        end
        
        axis square
        %     xlabel('distance (sessions)');
        xlabel(['\Delta ' distStr{idistType}]);
        ylabel(ylabelStr);
        %         xlim([1 size(distData,2)]);
        %         xticks([1 size(distData,2)]);
        %         xticklabels([1 size(distData,2)]);
        ylim([-0.5 0]);
        if i==1
            %             ylim([-0.15 -0.05]); %toZscore=1
            %             ylim([-0.23 -0.13]); %toZscore=0
        end
        
        %permutation histogram
        subplot(rows,cols,2*cols+(i-1)*length(distStr)+idistType)
        h=histogram(distCorrPerm(i,:),histBins); hold all
        h.FaceColor = 0.7*[1 1 1];
        h.Normalization = 'probability';
        axis square
        ylmt = get(gca,'ylim');
        plot([distCorr(i) distCorr(i)], [ylmt(1) ylmt(2)],'k','linewidth',linewidthWide);
        xlim([-1 1]);
        xlabel('correlation(r)');
        ylabel('permutations prob.');
        %     vline(distCorr(i));
        
        
        %generalization to adjacent sessions
        subplot(rows,cols,3*cols+(i-1)*length(distStr)+idistType); hold all
        
        idist=1;
        
        %compute correlation per subject and then average
        % scatter(subDistMatrix(subDistMatrix(:)>0),r2oriSplit{iroi}(subDistMatrix(:)>0))
        for isub=1:length(subjects)
            subData = similarityData(isub,1:minSessions,1:minSessions);
            dataDistSess = subData(distMatrix(1:minSessions,1:minSessions)==idist);
            [trainSess testSess] = ind2sub(size(distMatrix(1:minSessions,1:minSessions)),find(distMatrix(1:minSessions,1:minSessions)==idist));%make sure this is correct
            
            switch idistType
                case 1
                    trainSess = trainSess;
                case 2
                    trainSess = sessTime(isub,trainSess)';
                case 3
                    trainSess = sessHour(isub,trainSess)';
                case 4
                    trainSess = sessBias(isub,trainSess)';
                case 5
                    trainSess = sessDprime(isub,trainSess)';
                case 6
                    trainSess = adjHitRate(isub,trainSess)';
                case 7
                    trainSess = sessCorrect(isub,trainSess)';
            end
            
            [subDistSessCorr(isub,i) pSubDistSessCorr(isub,i)] = corr(trainSess,dataDistSess);
            
            %                 [subDistSessCorr(isub,i) pSubDistSessCorr(isub,i)] = corr(sessDistVec(sessDistVec>0),subData(sessDistVec>0));
            permDistMat = squeeze(similarityPerm(isub,:,:,:));
            for iperm=1:nperms
                tempPermMat = squeeze(permDistMat(iperm,1:minSessions,1:minSessions));
                permDistVec = tempPermMat(distMatrix(1:minSessions,1:minSessions)==idist);
                [subDistSessCorrPerm(i,isub,iperm), pSubDistSessCorrPerm(i,isub,iperm)] = corr(trainSess,permDistVec);
            end
            %plot
            scatter1=scatter(trainSess,dataDistSess,scatterSize,repmat(subColor{isub},length(trainSess),1),scatterFill);
            scatter1.MarkerFaceAlpha = 1;%4*surfaceAlpha;
            scatter1.MarkerEdgeAlpha = 1;%4*surfaceAlpha;
            hold on
            for ifirst=1:minSessions-1
                subMeanDistSessData(isub,ifirst) = mean(dataDistSess(trainSess(:)==ifirst));
            end
            plot(subMeanDistSessData(isub,:),'linewidth', linewidthNarrow,'color',subColor{isub}); hold all
            
            distSessCorr(i) = mean(subDistSessCorr(:,i),1);%mean over subjects
            distSessCorrPerm(i,:) = mean(subDistSessCorrPerm(i,:,:),2);%mean over subjects
            pPermDistSessCorr(i) = 1-sum(distSessCorrPerm(i,:)>=distSessCorr(i))/nperms;
            plot(mean(subMeanDistSessData,1),'linewidth', linewidthWide,'color','k'); hold all
            
            axis square
            xlabel(distStr{idistType});
            ylabel(ylabelStr);
            %         xlim([1 minSessions-idist]);
            %         xticks([1 minSessions-idist]);
            %         xticklabels([1 minSessions-idist]);
            
        end
        if ~saveFigs
            title(['r=' num2str(distSessCorr(i),'%4.2f') ' p=' num2str(pPermDistSessCorr(i),'%4.3f')]);
        end
        
        %adjacent permutation histogram
        subplot(rows,cols,4*cols+(i-1)*length(distStr)+idistType)
        h=histogram(distSessCorrPerm(i,:),histBins); hold all
        h.FaceColor = 0.7*[1 1 1];
        h.Normalization = 'probability';
        axis square
        ylmt = get(gca,'ylim');
        plot([distSessCorr(i) distSessCorr(i)], [ylmt(1) ylmt(2)],'k','linewidth',linewidthWide);
        xlim([-1 1]);
        xlabel('correlation(r)');
        ylabel('permutations prob.');
    end
end

for isubplot=1:rows*cols
    subplot(rows,cols,isubplot)
    %     set(ax2, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', []);
    set(gca, 'box', 'on', 'Visible', 'on');
end
set(gcf,'position',[250+90 200+20 160*cols 700]);


if saveFigs
%     savepdf(f,fullfile(figsFolder,['brief2_behavior' zscoreStr normalizeStr '.pdf']));
    saveas(h, [fullfile(figsFolder,['brief2_behavior' zscoreStr normalizeStr]) imgFormat]);
end

pSubPermDistCorr;
mean(subPartCorrBehav)
mean(subPartCorrSess)
%% plot the various factors
f=figure(ifig); ifig=ifig+1;  isubplot=0;
for isub=1:length(subjects)
    
    subCorrMat(isub,:,:) = corr([1:minSessions; sessTime(isub,1:minSessions); sessHour(isub,1:minSessions);...
        sessBias(isub,1:minSessions); sessDprime(isub,1:minSessions); ... 
        adjHitRate(isub,1:minSessions); sessCorrect(isub,1:minSessions)]');
end
imagesc(squeeze(mean(subCorrMat,1)));
xticklabels(distStr)
yticklabels(distStr)
colorbar
%% correlations between the various factors
f=figure(ifig); ifig=ifig+1;  isubplot=0;
cols = length(distStr)-1;
rows=1;
for isub=1:length(subjects)
    for idistType=2:length(distStr)
        switch idistType
            case 2
                subData = sessTime(isub,1:minSessions)';
            case 3
                subData = sessHour(isub,1:minSessions)';
            case 4
                subData = sessBias(isub,1:minSessions)';
            case 5
                subData = sessDprime(isub,1:minSessions)';
            case 6
                subData = adjHitRate(isub,1:minSessions)';
            case 7
                subData = sessCorrect(isub,1:minSessions)';
        end
        subplot(rows,cols,idistType-1)
        plot(subData,'color',subColor{isub},'linewidth',linewidthNarrow); hold on
        axis square
        title(distStr{idistType});
        xticks([1 30]);
        xlabel('session');
    end
end
set(gcf,'position',[250+90 200+20 160*cols 200]);