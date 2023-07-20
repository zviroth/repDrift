close all
clear all
rng(2);

saveFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift');
figsFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_figs');
if ~isfolder(saveFolder)
    saveFolder = ['/misc/data18/rothzn/nsd/repDrift/'];
    figsFolder = ['/misc/data18/rothzn/nsd/repDrift_figs/'];
end

% figsFolder = fullfile('/Users','rothzn','NSD',['repDrift_expand_figs']);
saveFigs=0;

linewidth=1.5;
zeroColor = 0.6*[1 1 1];
sessColor = {[180,243,192]./255,[255,223,186]./255};
graySubColor = 0.7*ones(1,3);
nsessions = 2;
T = 10;
gain=1;
origResp = gain*rand(T,nsessions)-gain/2;
% origResp = randn(T,nsessions);
for i=1:nsessions
    origResp(:,i) = zscore(origResp(:,i));
end
origResp(:,1) = origResp(:,1)*0.6 + 0.5;
origResp(:,2) = origResp(:,2)*1.5 - 0.5;
% origResp(:,3) = origResp(:,3)*2;
% origResp(:,4) = origResp(:,4)*0.4;
zeromean = origResp - mean(origResp);
normstd = (origResp - mean(origResp))./std(origResp) + mean(origResp);
f=figure(1);
rows=3;
cols=3;
subplot(rows,cols,1)
% plot(origResp)
plot(1:T,zeros(1,T),'color',zeroColor,'linewidth',1); hold on
for isess=1:nsessions
    plot(origResp(:,isess),'color',sessColor{isess},'linewidth',linewidth); hold all
end
axis square
subplot(rows,cols,2)
plot(1:T,zeros(1,T),'color',zeroColor,'linewidth',1); hold on
for isess=1:nsessions
    plot(zeromean(:,isess),'color',sessColor{isess},'linewidth',linewidth); hold all
end
axis square
subplot(rows,cols,3)
plot(1:T,zeros(1,T),'color',zeroColor,'linewidth',1); hold on
for isess=1:nsessions
    plot(normstd(:,isess),'color',sessColor{isess},'linewidth',linewidth); hold all
end
axis square
for isubplot=1:3
    subplot(rows,cols,isubplot)
    axis square
    xlim([1 T]);
    ylim([-3 3]);
    xlabel('trial #');
       if isubplot==1
    ylabel('response amp.');
       end
end


%% LOAD ACTIAL DATA AND PERMUTATIONS
if ieNotDefined('version'), version = 1; end%1=GLMdenoise, 2=without GLMdenoise, i.e. betas2.
subjects = [1:8];
% saveFigs=0;
% subjects = [1 2 5 7];
nperms=1000;

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
iroi = [1];

% saveFolder = fullfile('~','NSD',['repDrift_expand']);

normalizeStr = '';
if toNormalize
    normalizeStr = '_normalized';
end
r2thresh = 0;%33;%50;
r2threshStr = '';
if r2thresh>0
    r2threshStr = ['r2thresh' num2str(r2thresh,'%4.2f')];
end
zscoreVals=[0 2 3];
for izscore=1:length(zscoreVals)
    
    toZscore=zscoreVals(izscore);
    zscoreStr='';
    if toZscore==1
        zscoreStr = '_zscore';
    elseif toZscore==2
        zscoreStr = '_zeroMean';
    elseif toZscore==3
        zscoreStr = '_equalStd';
    end
    
    load(fullfile(saveFolder, ['perm' num2str(nperms) zscoreStr normalizeStr r2threshStr '.mat']),...
        'toNormalize','toZscore', 'useMedian','r2thresh','version','nrois','rois', ...
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
    
    %% use data for selected subjects
    subSessions = subSessions(subjects,:);
    minSessions = min(sum(subSessions,2));
    
    %% AUTOCORRELATION - Mean betas, and STD betas
    autocorrMin = -0.15;
    autocorrMax = 0.25;
    starHeight = 0.23; starSize=20; starColor = [0 0 0];
    
    
    data = meanAutocorrBetas;
%     titleStr = 'voxel mean betas';
    permData = meanAutocorrBetasPerm;
    
    %keep only data for chosen subjects
    data = data(:,subjects,:);
    permData = permData(:,subjects,:,:);
    
    subplot(rows,cols,cols+izscore);
    %     plot(squeeze(data(iroi,:,2:end))'); hold all
    for isub=1:length(subjects)
        plot(squeeze(data(iroi,isub,2:end))','color',subColor{subjects(isub)},'linewidth',linewidthNarrow); hold all
    % plot(squeeze(data(iroi,isub,2:end))','color',graySubColor,'linewidth',linewidthNarrow); hold all
    end
    empiricalMean = squeeze(mean(data(iroi,:,2:end),2));
    plot(empiricalMean','linewidth',linewidthWide,'color',[0 0 0]);
    ylim([autocorrMin autocorrMax]);
%     title(titleStr)
    axis square
    xlabel('lag (sessions)');
    if izscore==1
    ylabel('correlation (r)');
    end
    plot(0:size(data,3),zeros(1,1+size(data,3)),'color',0.5*[1 1 1]);
    xlim([1 size(data,3)-1]);
    xticks([1 size(data,3)-1]);
    xticklabels([1 size(data,3)-1]);
    
    permMean = squeeze(mean(permData(iroi,:,:,2:end),2));
    for t=1:length(empiricalMean)
        pPermDistAutocorr(izscore,t) = sum(permMean(:,t)>=empiricalMean(t))/nperms;
    end
    
    subplot(rows,cols,2*cols+izscore);
    plot(squeeze(mean(permData(iroi,:,:,2:end),2))','color',[0.8 0.8 0.8]); hold all%mean over subjects, plot all permutations
    plot(squeeze(mean(squeeze(mean(permData(iroi,:,:,2:end),2)),1))','linewidth',linewidthWide,'color',[0 0 0]); hold all%mean over subjects, then over permutations
    ylim([autocorrMin autocorrMax]);
    %         title(titleStr)
    axis square
    xlabel('lag (sessions)');
    if izscore==1
        ylabel('correlation (r)');
    end
    plot(0:size(data,3),zeros(1,1+size(data,3)),'color',0.5*[1 1 1]);
    xlim([1 size(data,3)-1]);
    xticks([1 size(data,3)-1]);
    xticklabels([1 size(data,3)-1]);
    
    %add significance stars
    subplot(rows,cols,cols+izscore);
    for t=1:length(empiricalMean)
        if pPermDistAutocorr(izscore,t)<=0.05
            scatter(t,starHeight,starSize,starColor,'.');
        end
        if pPermDistAutocorr(izscore,t)<=0.01
            scatter(t,starHeight-0.015,starSize,starColor,'.');
        end
    end
    
end


set(gcf,'position',[200 200 160*cols 500])
if saveFigs
%    savepdf(f,fullfile(figsFolder,['rep_fig2_left.pdf']));
   print('-painters','-dpdf',fullfile(figsFolder,['brief_figS3.pdf']));
%         saveas(h, [figFolder 'rep_fig2_left.' imgFormat]);
end
