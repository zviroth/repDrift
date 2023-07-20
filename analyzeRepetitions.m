close all hidden

clear all
tic

saveFigs = 0;

toZscore = 0;


nsdFolder = fullfile('~','misc','data18','rothzn','nsd');
if ~isfolder(nsdFolder)
    nsdFolder = '/misc/data18/rothzn/nsd/';
    figsFolder = '/misc/data18/rothzn/nsd/repDrift_expand_figs/';
end
subjects=1:8;
% subjects = [1:2];
% subjects = [3 4];
% subjects=[4 8];
% nsessionsSub = [40 40 32 30 40 32 40 30];
numSubs = length(subjects);
totalSubs = 8;
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
subRoiBetaDiff = NaN(totalSubs,numROIs,nsess,nreps);
subRoiBetaAbsDiff = NaN(totalSubs,numROIs,nsess,nreps);
subRoiSessMeanStd = NaN(totalSubs,numROIs,nsess,nreps);
subRoiDistCorrStd = NaN(totalSubs,numROIs,nsess,nreps);

roiRepMean = NaN(totalSubs,numROIs,nreps);
normRoiRepMean = NaN(totalSubs,numROIs,nreps);
subRoiCorr = NaN(totalSubs,numROIs,nreps);
normRepCorr = NaN(totalSubs,numROIs,nreps);

allRepBetas = NaN(totalSubs,numROIs,nsess,nreps);

nimgs = 10000;
nsdDesignFilename = fullfile(nsdFolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
% nsessionsSub = [40 40 32 30 40 32 40 30];
nsess=30;
ntrials = length(nsdDesign.masterordering);
splitImgTrials = zeros(nsess,ntrials);
maxReps=3;
sessRepTrials = zeros(nsess,ntrials,maxReps);%0 or 1 for each trial
trialsPerSession = 10 * 75; %12 runs, 75 trials per run
trialSession = zeros(1,ntrials);

stimRepNum = NaN(1,3*nimgs);%for each trial, what repetition is it
imgRepTrial = zeros(nimgs,3);%the 3 trials on which this image was presented
imgRepSess = zeros(nimgs,3);%the 3 sessions on which this image was presented
for i=1:nimgs%10000 images
    ind = find(nsdDesign.masterordering == i); %trials on which this image was presented. Should be 3 trials.
    imgRepTrial(i,:) = ind;
    imgRepSess(i,:) = ceil(ind/trialsPerSession);%session number for each repetition of each image
    for irep=1:length(ind)%should be 3
        stimRepNum(ind(irep)) = irep;
    end
end

itrial=1;
for isplit=1:nsess
    trialSession(itrial:itrial+trialsPerSession-1) = isplit;
    splitImgTrials(isplit,itrial:itrial+trialsPerSession-1) = ones;
    for irep=1:maxReps
        sessRepTrials(isplit,itrial:itrial+trialsPerSession-1,irep) = stimRepNum(itrial:itrial+trialsPerSession-1)==irep;
    end
    itrial = itrial+trialsPerSession;
end

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
   load(fullfile(saveFolder,['stimRepMean_sub' num2str(isub) zscoreStr '.mat']),...
    'visualRegions','combinedRoiBetas','nsess');
   ntrials = size(combinedRoiBetas{1},2);
   
   stimRepNum = stimRepNum(1:ntrials);%rep 1: 9209, rep 2: 7846, rep 3: 5445
%    trialSession = trialSession(1:ntrials);
   for iroi=1:numROIs
        nvox = size(combinedRoiBetas{iroi},1);

       %FIG 1
       roiRepBetas{iroi} = NaN(nsess,nvox,nreps);
       for isess=1:nsess
           for irep=1:nreps
               imgTrials = sessRepTrials(isess,:,irep);
               roiRepBetas{iroi}(isess,:,irep) = nanmean(combinedRoiBetas{iroi}(:,imgTrials>0),2);%mean over images
           end
       end
       allRepBetas(isub,iroi,:,:) = squeeze(nanmean(roiRepBetas{iroi},2));%mean over voxels
       

       %FIG 2
       imgBetas{iroi} = NaN(nimgs,maxReps,nvox);
       imgRoiMean{iroi} = NaN(nimgs,maxReps);
       imgRoiCorr{iroi} = NaN(nimgs,3);%correlations: 1-2, 1-3, 2-3
       imgRoiBetaDiff{iroi} = NaN(nimgs,3);%ROI mean beta differences: 1-2, 1-3, 2-3
       imgRepDist = NaN(nimgs,3);%distance in #session between repeitions 1-2, 1-3, 2-3
       for itrial=1:trialsPerSession*nsess %30000
           iimg = nsdDesign.masterordering(itrial);
           irep = stimRepNum(itrial);
           imgBetas{iroi}(iimg,irep,:) = combinedRoiBetas{iroi}(:,itrial);%1/4 of the values will be NaNs
           imgRoiMean{iroi}(iimg,irep) = nanmean(combinedRoiBetas{iroi}(:,itrial),1);%mean over voxels
       end
       for iimg=1:nimgs
           %corr
           imgRoiCorr{iroi}(iimg,:) = 1-pdist(squeeze(imgBetas{iroi}(iimg,:,:)),'correlation');%1-2, 1-3, 2-3
           imgRepDist(iimg,:) = [imgRepSess(iimg,2)-imgRepSess(iimg,1); imgRepSess(iimg,3)-imgRepSess(iimg,1); imgRepSess(iimg,3)-imgRepSess(iimg,2)];
           %beta diff
            imgRoiBetaDiff{iroi}(iimg,:) = [imgRoiMean{iroi}(iimg,1)-imgRoiMean{iroi}(iimg,2)...
                                            imgRoiMean{iroi}(iimg,1)-imgRoiMean{iroi}(iimg,3)...
                                            imgRoiMean{iroi}(iimg,2)-imgRoiMean{iroi}(iimg,3)];%1-2, 1-3, 2-3
       end

       for irep=1:3
%             voxRepMean{isub,iroi}(:,irep) = nanmean(combinedRoiBetas{iroi}(:,stimRepNum==irep),2);%mean over images
            roimean = mean(combinedRoiBetas{iroi}(:,stimRepNum==irep),1);%mean over voxels
            roiRepMean(isub,iroi,irep) = nanmean(roimean); %mean over images and voxels
            roiRepStd(isub,iroi,irep) = std(roimean); %std over images, mean over voxels
            subRoiCorr(isub,iroi,irep) = nanmean(imgRoiCorr{iroi}(:,irep),1);%mean over images
            for isess=1:nsess
                temp = imgRepSess(:,irep)==isess;%images whose irep repetition is in this session
                subRoiSessMean(isub,iroi,isess,irep) = nanmean(imgRoiMean{iroi}(temp,irep),1);%mean over images
                subRoiSessMeanStd(isub,iroi,isess,irep) = std(imgRoiMean{iroi}(temp,irep),0,1);%mean over images
            end
            for idist=0:nsess-1
                temp = imgRepDist(:,irep)==idist;
                %correlations
                subRoiDistCorr(isub,iroi,idist+1,irep) = nanmean(imgRoiCorr{iroi}(temp,irep),1);%mean over images
                subRoiDistCorrStd(isub,iroi,idist+1,irep) = nanstd(imgRoiCorr{iroi}(temp,irep),0,1);%mean over images
                %beta difference
                subRoiBetaDiff(isub,iroi,idist+1,irep) = nanmean(imgRoiBetaDiff{iroi}(temp,irep),1);%mean over images
                subRoiBetaAbsDiff(isub,iroi,idist+1,irep) = nanmean(abs(imgRoiBetaDiff{iroi}(temp,irep)),1);%mean over images
            end
       end
   end    
end

iroi=1;
subMeanBetas = nanmean(allRepBetas(:,:,:,:),1);
subStdBetas = nanstd(allRepBetas(:,:,:,:),1);

%% FIG 1
repColors = {[133,149,225]/255, [224,123,145]/255, [141,213,147]/255, [239,151,8]/255};
surfaceAlpha = 0.1;
linewidth=2;
ifig=1;
figure(ifig); ifig=ifig+1;
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

%% FIG 2
f=figure(ifig); ifig=ifig+1; clf
rows=1;
cols=4;
subplot(rows,cols,1)
for irep=1:3 
        plot(squeeze(nanmean(subRoiSessMean(:,iroi,:,irep),1)),'color',repColors{irep},'linewidth',linewidth); hold on
end
for irep=1:3   
    dsErrorsurface(1:nsess,squeeze(nanmean(subRoiSessMean(:,iroi,:,irep),1)),squeeze(nanstd(subRoiSessMean(:,iroi,:,irep)))./sqrt(numSubs), repColors{irep}, surfaceAlpha); hold on
end
for irep=1:3 
        plot(squeeze(nanmean(subRoiSessMean(:,iroi,:,irep),1)),'color',repColors{irep},'linewidth',linewidth); hold on
%     [repCorr(irep) repPval(irep)] = corr(squeeze(subRoiSessMean(iroi,:,irep))',(1:30)');
end
% plot(squeeze(mean(subRoiSessMean(:,iroi,:,:))));
legend('1','2','3');
xlabel('session')
ylabel('mean beta');
xticklabels([1 10:10:nsess])
xticks([1 10:10:nsess])
axis square
%t-tests
ireps1 = [1 1 2];
ireps2 = [2 3 3];
pvalCorr = NaN(nreps,nsess);
%within each session, t-tests between repetitions
for irep=1:3
    irep1=ireps1(irep);
    irep2 = ireps2(irep);
    for isess=1:nsess
        [h pvalCorr(irep,isess)] = ttest(squeeze(subRoiSessMean(:,iroi,isess,irep1)),squeeze(subRoiSessMean(:,iroi,isess,irep2)));
    end
end
%within each repetition, correlation with session number
for irep=1:3
    [repCorr(irep), pvalRepCorr(irep)] = corr(squeeze(mean(subRoiSessMean(:,iroi,:,irep),1)),[1:nsess]');
end

subplot(rows,cols,2)
for irep=1:3 
    plot(squeeze(nanmean(subRoiDistCorr(:,iroi,:,irep))),'linewidth',linewidth,'color',repColors{irep}); hold on
end
for irep=1:3 
%     plot(squeeze(nanmean(subRoiDistCorr(:,iroi,:,irep))),'linewidth',linewidth,'color',repColors{irep});
    dsErrorsurface(1:nsess,squeeze(nanmean(subRoiDistCorr(:,iroi,:,irep))),squeeze(nanstd(subRoiDistCorr(:,iroi,:,irep)))./sqrt(numSubs), repColors{irep}, surfaceAlpha); hold on
end
for irep=1:3 
    plot(squeeze(nanmean(subRoiDistCorr(:,iroi,:,irep))),'linewidth',linewidth,'color',repColors{irep});
end
xlabel('\Deltasession')
ylabel('correlation')
xticklabels([0:10:20 nsess-1])
xticks([1:10:21 nsess])
axis square
legend('1-2','1-3','2-3');


subplot(rows,cols,3)
plot(squeeze(nanmean(subRoiBetaDiff(:,iroi,:,:))),'linewidth',linewidth);
xlabel('\Deltasession')
xticklabels([0:10:20 nsess-1])
xticks([1:10:21 nsess])
ylabel('\Delta\beta');
axis square
legend('1-2','1-3','2-3');

subplot(rows,cols,4)
plot(squeeze(nanmean(subRoiBetaAbsDiff(:,iroi,:,:))),'linewidth',linewidth);
xlabel('\Deltasession')
xticklabels([0:10:20 nsess-1])
xticks([1:10:21 nsess])
ylabel('\Deltaabs(\beta)');
axis square

set(gcf,'position',[300 300 1100 200]);

if saveFigs
    savepdf(f,fullfile(figsFolder,['repetitionDrift_sess' zscoreStr '.pdf']));
end

%% ACROSS ALL SESSIONS
markersize = 10;
f=figure(ifig); ifig=ifig+1; clf
rows=1;
cols=4;
subplot(rows,cols,1)
for isub=subjects%1:length(subjects)
%     normRepMean(isub,:) = roiRepMean(isub,iroi,:)./roiRepMean(isub,iroi,1);
    plot(squeeze(roiRepMean(isub,iroi,:)),'.-','color',subColor{isub},'linewidth',linewidthNarrow,'markersize',markersize); hold all
end
plot(squeeze(nanmean(roiRepMean(:,iroi,:),1)),'.-','linewidth', linewidthWide,'color','k','markersize',markersize); hold all
xlabel('repetition');
ylabel('beta');


subplot(rows,cols,2)
for isub=subjects%length(subjects)
%     normRepMean(isub,:) = roiRepMean(isub,iroi,:)./roiRepMean(isub,iroi,1);
    plot(squeeze(subRoiCorr(isub,iroi,:)),'.-','color',subColor{isub},'linewidth',linewidthNarrow,'markersize',markersize); hold all
end
plot(squeeze(nanmean(subRoiCorr(:,iroi,:),1)),'.-','linewidth', linewidthWide,'color','k','markersize',markersize); hold all
xlabel('rep. pair');
ylabel('correlation');
xticklabels({'1-2','1-3','2-3'});

subplot(rows,cols,3)
for isub=subjects%1:length(subjects)
%     normRepMean(isub,:) = roiRepMean(isub,iroi,:)./roiRepMean(isub,iroi,1);
    normRoiRepMean(isub,iroi,:) = roiRepMean(isub,iroi,:)./roiRepMean(isub,iroi,1);
    plot(squeeze(normRoiRepMean(isub,iroi,:)),'.-','color',subColor{isub},'linewidth',linewidthNarrow,'markersize',markersize); hold all
end
plot(squeeze(nanmean(normRoiRepMean(:,iroi,:),1)),'.-','linewidth', linewidthWide,'color','k','markersize',markersize); hold all
xlabel('repetition');
ylabel('normalized beta');

subplot(rows,cols,4)
for isub=subjects%1:length(subjects)
    normRepCorr(isub,iroi,:) = subRoiCorr(isub,iroi,:)./subRoiCorr(isub,iroi,1);
    plot(squeeze(normRepCorr(isub,iroi,:)),'.-','color',subColor{isub},'linewidth',linewidthNarrow,'markersize',markersize); hold all
end
plot(squeeze(nanmean(normRepCorr(:,iroi,:),1)),'.-','linewidth', linewidthWide,'color','k','markersize',markersize); hold all
xlabel('rep. pair');
ylabel('normalized corr');
xticklabels({'1-2','1-3','2-3'});

set(gcf,'position',[200 200 500 200]);

if saveFigs
    savepdf(f,fullfile(figsFolder,['repetitionDrift' zscoreStr '.pdf']));
end


%t-tests, difference between repetitions
for i=1:3
    [h pvalMeanBeta(i)] = ttest(roiRepMean(:,iroi,ireps1(i)),roiRepMean(:,iroi,ireps2(i)));
    [h pvalMeanCorr(i)] = ttest(subRoiCorr(:,iroi,ireps1(i)),subRoiCorr(:,iroi,ireps2(i)));
end
pvalMeanBeta
pvalMeanCorr

%% FIND IMAGES FOR AN RDM

for iimg=1:10000
    imgSess = imgRepSess(iimg,:);
    sessDiff = abs(imgRepSess - imgSess);
    sumDiff = sum(sessDiff,2);
    temp = find(sumDiff==0);
    sameSessImg{iimg} = temp;
    sameSessNum(iimg) = length(temp);
    imgNumDiffSess(iimg) = length(unique(imgSess));
    imgUniqueSess(iimg) = length(unique(imgSess));
end
[maxImagesSameSessions imageNum] = max(sameSessNum(imgUniqueSess==nreps));
% hist(sameSessNum(imgUniqueSess==nreps))

rdmImages = find(sameSessNum==maxImagesSameSessions & imgUniqueSess==nreps);
%there are 2 sets of 6 images each:
set1 = [2376        2439        6172        6838        8732        9178];%sameSessImg{2376}'. No
set2 = [2896        3657        3921        5656        6937        8670];%sameSessImg{2896}'
%there are 2 sets of 5 images each:
set1 = [678        2115        3024        5739        6976];%WORKS. sameSessImg{678}.
set2 = [1613        2000        2875        5999        8579];%NO. sameSessImg{1613}

%imgRepSess(set1,:) == 13, 14, 27
for irep=1:nreps
    repRDM(irep,:) = pdist(combinedRoiBetas{iroi}(:,imgRepTrial(set1,irep))','correlation');
end
repRDMsim = 1-pdist(repRDM);%1-2, 1-3, 2-3

% combinedRoiBetas{iroi}(:,imgRepTrial(set2,:));


%% 

figure(ifig); ifig=ifig+1;
% histogram(imgRepSess,40);
for i=1:3; histogram(imgRepSess(:,i),40,'FaceColor',subColor{i}); hold on; end
xlabel('session');
ylabel('#images');
legend('1st','2nd','3rd');
axis square
set(gcf,'position',[200 220 200 200]);

figure(ifig); ifig=ifig+1;
for i=1:2; histogram(imgRepSess(:,i),40,'FaceColor',subColor{i}); hold on; end
xlabel('session');
ylabel('#images');
legend('1st','2nd','3rd');
axis square
set(gcf,'position',[250 220 200 200]);

figure(ifig); ifig=ifig+1;
for i=1:1; histogram(imgRepSess(:,i),40,'FaceColor',subColor{i}); hold on; end
xlabel('session');
ylabel('#images');
legend('1st','2nd','3rd');
axis square
set(gcf,'position',[300 220 200 200]);

%only 2nd presentation
figure(ifig); ifig=ifig+1;
for i=2:2; histogram(imgRepSess(:,i),40,'FaceColor',subColor{i}); hold on; end
xlabel('session');
ylabel('#images');
% legend('1st','2nd','3rd');
legend('2nd');
axis square
set(gcf,'position',[350 220 200 200]);

%only 3rd presentation
figure(ifig); ifig=ifig+1;
for i=3:3; histogram(imgRepSess(:,i),40,'FaceColor',subColor{i}); hold on; end
xlabel('session');
ylabel('#images');
% legend('1st','2nd','3rd');
legend('3rd');
axis square
set(gcf,'position',[400 220 200 200]);


%% find number of unique images shown in each session. Do high repeats around session 25 cause the band?
%get the image for each trial: nsdDesign.masterordering
%get the session for each trial
%within each session count unique images
trialSess = ceil([1:ntrials]/trialsPerSession);%session number for each trial
trialImg = nsdDesign.masterordering(1:ntrials);
for isess=1:nsess
    temp = unique(trialImg(trialSess==isess));
    uniqueImgs(isess) = length(temp);
end
figure(ifig); ifig=ifig+1;
plot(uniqueImgs,'linewidth',linewidth)
xlabel('session')
ylabel('#unique images');
axis square
set(gcf,'position',[400 220 200 200]);

toc