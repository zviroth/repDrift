close all
clear all

saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand');
figsFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand_figs');
if ~isfolder(saveFolder)
    saveFolder = ['/misc/data18/rothzn/nsd/repDrift_expand/'];
    figsFolder = ['/misc/data18/rothzn/nsd/repDrift_expand_figs/'];
end

saveFigs=0;
ifig=0;
subColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840],...
    [0.4660 0.6740 0.1880]/2+[0.8500 0.3250 0.0980]/2};
graySubColor = 0.7*ones(1,3);
colorFractions = [0.1 0.5 1];
for i=1:length(colorFractions)
    plotColors{i} = subColor{1}*colorFractions(i)*1.1 + subColor{2}*(1-colorFractions(i)*1.1);
    plotColors{i}(plotColors{i}>1) = ones;
    plotColors{i}(plotColors{i}<0) = zeros;
end
plotColors = {[255,179,186]./255, 0.5*([255,179,186]+[186,225,255])./255, [186,225,255]./255};
barFacecolor = 0.8*ones(1,3);
barEdgecolor = 0.2*ones(1,3);
zeroColor= 0.6*[1 1 1];
linewidth=2;


cols = 2;
rows=2;


ifig=ifig+1; f=figure(ifig); clf;

%% LOAD AUTOCORR DATA AND PERMUTATIONS
if ieNotDefined('version'), version = 1; end%1=GLMdenoise, 2=without GLMdenoise, i.e. betas2.
subjects = [1:8];

% subjects = [1 2 5 7];
nperms=1000;

singleSubject=1;
figRoi=1;
linewidthWide=2;
linewidthNarrow = 1;
errorbarColor = [0.2 0.2 0.2];
surfaceAlpha = 0.1;

iroi = [1];
% saveFolder = fullfile('~','NSD',['repDrift' versionStr]);

r2thresh = 0;%33;%50;
r2threshStr = '';
if r2thresh>0
    r2threshStr = ['r2thresh' num2str(r2thresh,'%4.2f')];
end
toZscore=0;
zscoreStr='';
if toZscore==1
    zscoreStr = '_zscore';
elseif toZscore==2
    zscoreStr = '_zeroMean';
elseif toZscore==3
    zscoreStr = '_equalStd';
end
load(fullfile(saveFolder, ['perm' num2str(nperms) zscoreStr  r2threshStr '.mat']),...
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
    'subRoiPrf','numGoodVox');
nsubjects = length(subjects);
% use data for selected subjects
subSessions = subSessions(subjects,:);
minSessions = min(sum(subSessions,2));

%% AUTOCORRELATION - Mean betas, and STD betas
autocorrMin = -0.15;
autocorrMax = 0.25;
starHeight = 0.23; starSize=20; starColor = [0 0 0];

for i=1:cols
    switch i
        
        %         case 1
        %             data = meanAutocorrConstant;
        %             titleStr = 'constrained constant';
        %             permData = meanAutocorrConstantPerm;
        case 1
            data = meanAutocorrConstantOri;
            titleStr = 'constant';
            permData = meanAutocorrConstantOriPerm;
            %         case 3
            %             data = meanAutocorrCoef;
            %             titleStr = 'constrained std coef';
            %             permData = meanAutocorrCoefPerm;
        case 2
            data = meanAutocorrCoefOri;
            titleStr = 'std coef';
            permData = meanAutocorrCoefOriPerm;
    end
    %keep only data for chosen subjects
    data = data(:,subjects,:);
    %     if i<3 | i>4
    permData = permData(:,subjects,:,:);
    %     end
    
    
    subplot(rows,cols,i)
    %     plot(squeeze(data(iroi,:,2:end))'); hold all
    for isub=1:length(subjects)
        plot(squeeze(data(iroi,isub,2:end))','color',subColor{subjects(isub)},'linewidth',linewidthNarrow); hold all
        %         plot(squeeze(data(iroi,isub,2:end))','color',graySubColor,'linewidth',linewidthNarrow); hold all
    end
    empiricalMean = squeeze(mean(data(iroi,:,2:end),2));
    plot(empiricalMean','linewidth',linewidthWide,'color',[0 0 0]);
    ylim([autocorrMin autocorrMax]);
    if ~saveFigs
        title(titleStr);
    end
    axis square
    xlabel('lag (sessions)');
    if i==1
        ylabel('correlation (r)');
    end
    plot(0:size(data,3),zeros(1,1+size(data,3)),'color',0.5*[1 1 1]);
    xlim([1 size(data,3)-1]);
    xticks([1 size(data,3)-1]);
    xticklabels([1 size(data,3)-1]);
    
    permMean = squeeze(mean(permData(iroi,:,:,2:end),2));
    for t=1:length(empiricalMean)
        pPermDistAutocorr(i,t) = sum(permMean(:,t)>=empiricalMean(t))/nperms;
    end
    
    subplot(rows,cols,cols+i)
    plot(squeeze(mean(permData(iroi,:,:,2:end),2))','color',[0.8 0.8 0.8]); hold all%mean over subjects, plot all permutations
    plot(squeeze(mean(squeeze(mean(permData(iroi,:,:,2:end),2)),1))','linewidth',linewidthWide,'color',[0 0 0]); hold all%mean over subjects, then over permutations
    ylim([autocorrMin autocorrMax]);
    %         title(titleStr)
    axis square
    xlabel('lag (sessions)');
    if i==1
        ylabel('correlation (r)');
    end
    plot(0:size(data,3),zeros(1,1+size(data,3)),'color',0.5*[1 1 1]);
    xlim([1 size(data,3)-1]);
    xticks([1 size(data,3)-1]);
    xticklabels([1 size(data,3)-1]);
    
    %add significance stars
    subplot(rows,cols,i)
    for t=1:length(empiricalMean)
        if pPermDistAutocorr(i,t)<0.05
            scatter(t,starHeight,starSize,starColor,'.');
        end
        if pPermDistAutocorr(i,t)<0.01
            scatter(t,starHeight-0.02,starSize,starColor,'.');
        end
    end
end
% set(gcf,'position',[500 300 800 350]);
set(gcf,'position',[50 300 160*cols 300]);
if saveFigs
    %    savepdf(f,fullfile(figsFolder,['fig2' zscoreStr normalizeStr '.pdf']));
    print('-painters','-dpdf',fullfile(figsFolder,['brief_figS5.pdf']));
end