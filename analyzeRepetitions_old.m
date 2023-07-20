close all
clear all
tic
saveFigs = 0;
nsdFolder = fullfile('~','misc','data18','rothzn','nsd');
if ~isfolder(nsdFolder)
    nsdFolder = '/misc/data18/rothzn/nsd/';
    figsFolder = '/misc/data18/rothzn/nsd/repDrift_expand_figs/';
end
% subjects=1:8;
subjects = [1:2];
subjects = [3 4];
% subjects=[4 8];
% nsessionsSub = [40 40 32 30 40 32 40 30];
numSubs = length(subjects);
totalSubs = 2;
saveFolder = fullfile(nsdFolder,'stimRepetitions','/');
subColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840],...
    [0.4660 0.6740 0.1880]/2+[0.8500 0.3250 0.0980]/2};
linewidthNarrow=1;
linewidthWide=3;

nreps=3;
nsess = 30;%30;
numROIs = 1;
subRoiSessMean = NaN(totalSubs,numROIs,nsess,nreps);
subRoiDistCorr = NaN(totalSubs,numROIs,nsess,nreps);
subRoiSessMeanStd = NaN(totalSubs,numROIs,nsess,nreps);
subRoiDistCorrStd = NaN(totalSubs,numROIs,nsess,nreps);

roiRepMean = NaN(totalSubs,numROIs,nreps);
normRoiRepMean = NaN(totalSubs,numROIs,nreps);
subRoiCorr = NaN(totalSubs,numROIs,nreps);
normRepCorr = NaN(totalSubs,numROIs,nreps);

allRepBetas = NaN(totalSubs,numROIs,nsess,nreps);

for isub=subjects
   load(fullfile(saveFolder,['stimRepMean_sub' num2str(isub) '.mat']),...
    'roiMeanBetas','roiStdBetas','roiRepBetas','roiRepStdBetas',...
    'visualRegions',...
    'stimRepNum','combinedRoiBetas','trialSession',...
    'imgBetas','imgRoiMean','imgRoiCorr',...
    'imgRepTrial','imgRepDist',...
    'imgRepSess','splitImgTrials','sessRepTrials');
   ntrials = size(combinedRoiBetas{1},2);
   stimRepNum = stimRepNum(1:ntrials);%rep 1: 9209, rep 2: 7846, rep 3: 5445
   trialSession = trialSession(1:ntrials);
   for iroi=1:numROIs
       allRepBetas(isub,iroi,:,:) = squeeze(nanmean(roiRepBetas{iroi},2));%mean over voxels
       for irep=1:3
            voxRepMean{isub,iroi}(:,irep) = nanmean(combinedRoiBetas{iroi}(:,stimRepNum==irep),2);%mean over images
            
            roimean = mean(combinedRoiBetas{iroi}(:,stimRepNum==irep),1);%mean over voxels
            roiRepMean(isub,iroi,irep) = nanmean(roimean); %mean over images and voxels
            roiRepStd(isub,iroi,irep) = std(roimean); %std over images, mean over voxels
            subRoiCorr(isub,iroi,irep) = nanmean(imgRoiCorr{iroi}(:,irep),1);%mean over images
            for isess=1:nsess
                temp = imgRepSess(:,irep)==isess;%images whose irep repetition is in this session
                subRoiSessMean(isub,iroi,isess,irep) = nanmean(imgRoiMean{iroi}(temp),1);%mean over images
                subRoiSessMeanStd(isub,iroi,isess,irep) = std(imgRoiMean{iroi}(temp),0,1);%mean over images
            end
            for idist=0:nsess-1
                temp = imgRepDist(:,irep)==idist;
                subRoiDistCorr(isub,iroi,idist+1,irep) = nanmean(imgRoiCorr{iroi}(temp),1);%mean over images
                subRoiDistCorrStd(isub,iroi,idist+1,irep) = nanstd(imgRoiCorr{iroi}(temp),0,1);%mean over images
            end
            
       end
   end    
end

iroi=1;
subMeanBetas = nanmean(allRepBetas(:,:,:,:),1);
subStdBetas = nanstd(allRepBetas(:,:,:,:),1);

%% PER SESSION
repColors = {[133,149,225]/255, [224,123,145]/255, [141,213,147]/255, [239,151,8]/255};
surfaceAlpha = 0.1;
linewidth=2;
figure(1)
clf
% plot(squeeze(subMeanBetas(iroi,:,:)))

%for legend:
for irep=1:3    
    plot(squeeze(subMeanBetas(:,iroi,:,irep)),'color',repColors{irep},'linewidth',linewidth); hold on
end
for irep=1:3   
    dsErrorsurface(1:nsess,squeeze(subMeanBetas(:,iroi,:,irep)),squeeze(subStdBetas(:,iroi,:,irep))./sqrt(numSubs), repColors{irep}, surfaceAlpha); hold on
end
for irep=1:3
    plot(squeeze(subMeanBetas(:,iroi,:,irep)),'color',repColors{irep},'linewidth',linewidth); hold on
    [repCorr(irep) repPval(irep)] = corr(squeeze(subMeanBetas(:,iroi,:,irep)),(1:nsess)');
end
legend('1','2','3');


%% ACROSS ALL SESSIONS
markersize = 10;
f=figure(2); clf
rows=1;
cols=4;
subplot(rows,cols,1)
for isub=subjects%1:length(subjects)
%     normRepMean(isub,:) = roiRepMean(isub,iroi,:)./roiRepMean(isub,iroi,1);
    plot(squeeze(roiRepMean(isub,iroi,:)),'.-','color',subColor{isub},'linewidth',linewidthNarrow,'markersize',markersize); hold all
end
plot(squeeze(nanmean(roiRepMean(:,iroi,:),1)),'.-','linewidth', linewidthWide,'color','k','markersize',markersize); hold all
xlabel('repetition');
ylabel('beta (% signal change)');

subplot(rows,cols,2)
for isub=subjects%1:length(subjects)
%     normRepMean(isub,:) = roiRepMean(isub,iroi,:)./roiRepMean(isub,iroi,1);
    normRoiRepMean(isub,iroi,:) = roiRepMean(isub,iroi,:)./roiRepMean(isub,iroi,1);
    plot(squeeze(normRoiRepMean(isub,iroi,:)),'.-','color',subColor{isub},'linewidth',linewidthNarrow,'markersize',markersize); hold all
end
plot(squeeze(nanmean(normRoiRepMean(:,iroi,:),1)),'.-','linewidth', linewidthWide,'color','k','markersize',markersize); hold all
xlabel('repetition');
ylabel('normalized beta');

subplot(rows,cols,3)
for isub=subjects%length(subjects)
%     normRepMean(isub,:) = roiRepMean(isub,iroi,:)./roiRepMean(isub,iroi,1);
    plot(squeeze(subRoiCorr(isub,iroi,:)),'.-','color',subColor{isub},'linewidth',linewidthNarrow,'markersize',markersize); hold all
end
plot(squeeze(nanmean(subRoiCorr(:,iroi,:),1)),'.-','linewidth', linewidthWide,'color','k','markersize',markersize); hold all
xlabel('corr');
ylabel('correlation');

subplot(rows,cols,4)
for isub=subjects%1:length(subjects)
    normRepCorr(isub,iroi,:) = subRoiCorr(isub,iroi,:)./subRoiCorr(isub,iroi,1);
    plot(squeeze(normRepCorr(isub,iroi,:)),'.-','color',subColor{isub},'linewidth',linewidthNarrow,'markersize',markersize); hold all
end
plot(squeeze(nanmean(normRepCorr(:,iroi,:),1)),'.-','linewidth', linewidthWide,'color','k','markersize',markersize); hold all
xlabel('repetition');
ylabel('normalized corr');

set(gcf,'position',[200 200 400 250]);

if saveFigs
    savepdf(f,fullfile(figsFolder,['repetitionBetas.pdf']));
end
%% ACROSS ALL SESSIONS, AS FUNCTION OF DISTANCE
f=figure(3); clf
rows=1;
cols=2;
subplot(rows,cols,1)
for irep=1:3   
    dsErrorsurface(1:nsess,squeeze(nanmean(subRoiSessMean(:,iroi,:,irep),1)),squeeze(nanstd(subRoiSessMean(:,iroi,:,irep)))./sqrt(numSubs), repColors{irep}, surfaceAlpha); hold on
end
for irep=1:3 
        plot(squeeze(nanmean(subRoiSessMean(:,iroi,:,irep),1)),'color',repColors{irep},'linewidth',linewidth); hold on
%     [repCorr(irep) repPval(irep)] = corr(squeeze(subRoiSessMean(iroi,:,irep))',(1:30)');
end
% plot(squeeze(mean(subRoiSessMean(:,iroi,:,:))));

xlabel('session')
xticklabels([1 10:10:nsess])
xticks([1 10:10:nsess])
axis square

subplot(rows,cols,2)
plot(squeeze(nanmean(subRoiDistCorr(:,iroi,:,:))),'linewidth',linewidth);
xlabel('\Deltasession')
xticklabels([0:10:20 nsess-1])
xticks([1:10:21 nsess])
axis square

set(gcf,'position',[300 300 800 300]);

toc