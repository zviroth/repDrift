close all
clear all

figsFolder = fullfile('~','misc','data18','rothzn','nsd','repDrift_expand_figs');
if ~isfolder(figsFolder)
    figsFolder = ['/misc/data18/rothzn/nsd/repDrift_expand_figs/'];
end

saveFigs=0;

dx=0.01;
xmax = 5;
x = -xmax:dx:xmax;
sigma = 1;
x0=0;
amp=1;
orig = amp*(1/(sigma*sqrt(2*pi)))*exp(-0.5*((x-x0)/sigma).^2);
ifig = 0;
col=0;

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

%baseline 
col=col+1;
vals(col,:) = [-0.04 0 0.04];
for j=1:size(vals,2)
    data(col,j,:) = vals(col,j)+orig;
end
titleStr{col} = 'baseline';

%gain 
col=col+1;
gain=0.2;
vals(col,:) = [1-gain 1 1+gain];
for j=1:size(vals,2)
    data(col,j,:) = vals(col,j)*orig;
end
titleStr{col} = 'gain';

%tuning width
col=col+1;
vals(col,:) = sigma*[0.8 1 1.2];
for j=1:size(vals,2)
    data(col,j,:) = amp*(1/(vals(col,j)*sqrt(2*pi)))*exp(-0.5*((x-x0)/vals(col,j)).^2);
end
titleStr{col} = 'width';

%shift
col=col+1;
vals(col,:) = round(length(x)*[-0.06 0 0.06]);
for j=1:size(vals,2)
    data(col,j,:) = circshift(orig,vals(col,j));
end
titleStr{col} = 'preference';

%randomly sample the x dimension
% numSamples = 100;
numSamples = length(x);%NO RANDOM SAMPLING!!

randOrder = randperm(length(x));
randSample = randOrder(1:numSamples);

meanData = mean(data(:,:,randSample),3);
stdData = std(data(:,:,randSample),0,3);

cols = col;
rows=3;


ifig=ifig+1; f=figure(ifig); clf;
%PLOT DATA
padding = 0.05;
for i=1:cols
    subplot(rows,cols,i)
    plot(1:length(x),zeros(1,length(x)),'color',zeroColor,'linewidth',1); hold on
    for j=1:size(vals,2)
        plot(squeeze(data(i,j,:)),'color',plotColors{j},'linewidth',linewidth); hold on
    end
    axis square
    ymin = min(data(:));
    ymax = max(data(:));
    ylim([ymin-padding ymax+padding]);
    title(titleStr{i});
    set(gca,'xticklabels',[]); 
    set(gca,'yticklabels',[]); 
    xlabel('stimulus');
%     xlabel('stimulus parameter');
    if i==1
       ylabel('response'); 
    end
end
% l=legend({'1', '2', '3'},'location','northeast');
l=legend({'', '1', '2', '3'},'location','northeast');


%PLOT MEAN
for i=1:cols
    subplot(rows,cols,cols+i)
    bar(squeeze(meanData(i,:)),'facecolor',barFacecolor,'edgecolor',0.2*ones(1,3));
    axis square
    ymin = min(meanData(:));
    ymax = max(meanData(:));
    ylim([ymin-padding ymax+padding]);
    yticks([]);
%     title('mean');
    if i==1
       ylabel('response mean'); 
    end
    xlabel('condition');
end


%PLOT STD
for i=1:cols
    subplot(rows,cols,2*cols+i)
    bar(squeeze(stdData(i,:)),'facecolor',0.8*ones(1,3),'edgecolor',0.2*ones(1,3));
    axis square
    ymin = min(stdData(:));
    ymax = max(stdData(:));
    ylim([ymin-padding ymax+padding]);
    yticks([]);
%     title('std');
    if i==1
       ylabel('response std'); 
    end
    xlabel('condition');
end


% set(gcf,'position',[500 300 800 350]);
set(gcf,'position',[50 300 160*cols 400]);
if saveFigs
%    savepdf(f,fullfile(figsFolder,['fig2' zscoreStr normalizeStr '.pdf']));
    print('-painters','-dpdf',fullfile(figsFolder,['brief_figS4.pdf']));
end