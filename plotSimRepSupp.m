close all
clear all

repSuppFactor=0.9;
toZscore=4;
iroi=1;
nsessions = 30;
minSessions=30;
subjects=1:8;

distMatrix = toeplitz(0:nsessions-1);%nsplits = nsessions+1. from 0 to 39.
lagMatrix = toeplitz(0:-1:1-nsessions,0:nsessions-1);%from -39 to +39
nsdFolder = fullfile('~','misc','data18','rothzn','nsd');
if ~isfolder(nsdFolder)
    nsdFolder = '/misc/data18/rothzn/nsd/';
end
saveFolder = fullfile(nsdFolder,'stimRepetitions','/');

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
repSuppFactorStr = num2str(repSuppFactor,'%01.1f');
repSuppFactorStr(2)=[];

for isub=subjects
    load(fullfile(saveFolder,['simRepSupp_' repSuppFactorStr '_sub' num2str(isub) '_v' num2str(iroi) zscoreStr '.mat']),...
        'r2OrisplitVox','rssOrisplitVox','pearsonRoriVox',...
        'scalingAmp','noiseAmp');

    %take median over voxels
    rssOriSplit(isub,1:nsessions,1:nsessions) = squeeze(nanmedian(rssOrisplitVox(1:nsessions,:,1:nsessions),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
    r2oriSplit(isub,1:nsessions,1:nsessions) = squeeze(nanmedian(r2OrisplitVox(1:nsessions,:,1:nsessions),2));%nsplits x nsplits matrix, isplit,vox,nextSplit
    pearsonRori(isub,1:nsessions,1:nsessions) = squeeze(nanmedian(pearsonRoriVox(1:nsessions,:,1:nsessions),2));%nsplits x nsplits matrix, isplit,vox,nextSplit

    for idist=1:nsessions
        temp = squeeze(rssOriSplit(isub,1:minSessions,1:minSessions));
        rssOriDistSess{idist}(isub,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
        rssOriDist(isub,idist) = nanmean(temp(distMatrix==idist));

        temp = squeeze(r2oriSplit(isub,1:minSessions,1:minSessions));
        r2OriDistSess{idist}(isub,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
        r2OriDist(isub,idist) = nanmean(temp(distMatrix==idist));

        temp = squeeze(pearsonRori(isub,1:minSessions,1:minSessions));
        pearsonOriDistSess{idist}(isub,:) = 0.5*(temp(lagMatrix==idist)+temp(lagMatrix==-idist));
        pearsonOriDist(isub,idist) = nanmean(temp(distMatrix==idist));

    end

end
%%
rows=1;cols=2;
figure(1)
subplot(rows,cols,1)
plot(r2OriDist'); axis square
hold on
plot(mean(r2OriDist),'linewidth',2,'color','k')
xlabel('\Deltasession')
ylabel('cvR^2');
subplot(rows,cols,2)
plot(pearsonOriDist'); axis square
hold on
plot(mean(pearsonOriDist),'linewidth',2,'color','k')
xlabel('\Deltasession')
ylabel('correlation (r)');

set(gcf,'position',[300 400 400 150]);