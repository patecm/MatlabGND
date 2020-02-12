%% Various types of plots for the GND Code %%
% to add grain boundaries to a plot: hold on; plot(grains.boundary,'linewidth',1); hold off
%rotate plots with: plotx2west; plotzOutOfPlane; 

%% --- Grain orientation and boundaries plot ---
ipfKey = ipfHSVKey(ebsd);
ipfKey.inversePoleFigureDirection = yvector;
plot(ebsd,ipfKey.orientation2color(ebsd.orientations),'micronBar','off')
hold on
plot(grains.boundary,'linewidth',1)
hold off

%% Plot the grain with their corresponding number example
% % only the very big grains
% all_grains = grains(grains.grainSize>3);
% % plot them
% plot(all_grains,all_grains.meanOrientation,'micronbar','off')
% % plot on top their ids
% text(all_grains,int2str(all_grains.id))

%% --- Plot a property (like Average GND Density) per grain --- %
setaxis = [13.5 17]; %change axis range here
prop=rho_total; % any property with same size as ebsd
% calc grains and ebsd.grain ids:
[grains, ebsd.grainId,ebsd.mis2mean]=calcGrains(ebsd);
%find ebsd.grainID and index of corresponding property
[gid,~,eindex] = unique(ebsd.grainId);
% calc the mean for each grain, ignoring nans
gProp = accumarray(eindex,prop,[],@nanmean);
% or max: or change to whatever you want
%gProp = accumarray(eindex, prop,[],@max);
%add text values to big-ish grains
selected_grains = grains(grains.grainSize>5);
gProp2 = gProp(selected_grains.id);
%plot
figure;plotx2west; plotzOutOfPlane;
plot(grains,log10(gProp)); mtexColorMap('jet'); 
set(gca,'CLim',[setaxis]);
colorbar; 
text(selected_grains,num2str(gProp2, '%10.3e'), 'Color','white','FontSize',12, 'FontWeight', 'bold') %optional add values
%axis image % fit axes to image coordinates

%% --- Plot as pixel image (similar to Asher's plotdensities) --- %
%input
plotme = rho_total;
setaxis = [13.5 15];
%plot
plotme = reshape(log10(plotme), size(ebsd.id,1), size(ebsd.id,2));
hold on
plot(grains.boundary,'linewidth',1)
imagesc(plotme);
colormap(jet);
colorbar;
caxis(setaxis);
set(gca,'YDir','reverse');
axis image % fit axes to image coordinates
%% Plot two things side by side (eg total GND and per grain average)
%inputs
setaxis = [13.5 15];
val1 = rho_total;
val2 = gProp; %run 'Plot property' first
title1 = 'GND Density per pixel';
title2 = 'GND Density average per grain';
%plots
plotx2west; plotzOutOfPlane;
plot(ebsd,log10(val1)); 
hold on; plot(grains.boundary,'linewidth',1); hold off
title(title1); caxis(setaxis);nextAxis ; 
plot(grains,log10(val2)); %change grains to ebsd for some data
text(selected_grains,num2str(gProp2, '%10.3e'), 'Color','white','FontSize',10, 'FontWeight', 'bold') %optional add values
title(title2);caxis(setaxis);
drawNow(gcm,'figSize','large');
mtexColorMap('jet'); mtexColorbar;
axis image % fit axes to image coordinates

%% 1x3 plot (e.g. compare multiple rho or gnd energy values)
%input
val1 = rho_total;
val2 = rho_screw;
val3 = rho_edge;
title1 = 'Total GND energy';
title2 = 'GND energy for screw dislocations';
title3 = 'GND energy for edge dislocations';
setaxis = [13 15];
%plots
plotx2west; plotzOutOfPlane;
plot(ebsd,log10(val1)); %can also plot not on log scale
title(title1);caxis(setaxis);
nextAxis 
plot(ebsd,log10(val2));
title(title2); caxis(setaxis);
nextAxis 
plot(ebsd,log10(val3));
title(title3); caxis(setaxis);
drawNow(gcm,'figSize','large');
mtexColorMap('jet');mtexColorbar;

%% Generic 3x3 plot (example: for with Rank 2 tensors like kappa or alpha)
makeplot = kappa; %change this to whatever you want to plot
%plots
plotx2west; plotzOutOfPlane;
plotmtexFig = newMtexFigure('nrows',3,'ncols',3,'figSize','huge');
plot(ebsd,makeplot{1,1},'micronBar','off');
hold on; plot(grains.boundary,'linewidth',1); hold off;
nextAxis
plot(ebsd,makeplot{1,2},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',1); hold off;
nextAxis
plot(ebsd,makeplot{1,3},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',1); hold off;
nextAxis
plot(ebsd,makeplot{2,1},'micronBar','off');
hold on; plot(grains.boundary,'linewidth',1); hold off
nextAxis
plot(ebsd,makeplot{2,2},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',1); hold off
nextAxis
plot(ebsd,makeplot{2,3},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',1); hold off
nextAxis
plot(ebsd,makeplot{3,1},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',1); hold off
nextAxis
plot(ebsd,makeplot{3,2},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',1); hold off
nextAxis
plot(ebsd,makeplot{3,3},'micronBar','off')
hold on; plot(grains.boundary,'linewidth',1); hold off
drawNow(gcm,'figSize','large')
setColorRange([-0.05,0.05],'zero2white')

%% Pole figures
odf = unimodalODF(orientation.byEuler(0,0,0,cs));
ori = calcOrientations(odf,100);
ori_rotated = calcOrientations(rotate(odf,rotation.byEuler(60*degree,60*degree,0*degree)),100);

h = [Miller(0,0,0,1,cs),Miller(1,0,-1,0,cs)];
plotPDF(ori,h,'antipodal','MarkerSize',4)
hold all % keep plot
plotPDF(ori_rotated,h,'MarkerSize',4);
hold off % next plot command deletes all plots

%% Sample vs Uniform ODF
%input
mdf = calcMDF(odf);
odf2 = uniformODF(orientation.byEuler(0,0,0,cs));
mdf2 = calcMDF(odf2);
%plots
close all
plotAngleDistribution(mdf);
%hold all
%plotAngleDistribution(mdf2);
hold all
plotAngleDistribution(cs,cs);
hold off
legend('model ODF','Uniform2 ODF', 'Uniform ODF');

%% Histograms - GND densities
%inputs
pd = fitdist(log10(rho_total),'Normal');
mu = pd.mu;
sigma = pd.sigma;
f = (exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi)));
%plots
h1 = histogram(log10(rho_total));
h1.Normalization = 'pdf';
hold on
h2 = histogram(log10(rho_edge),197);
h2.Normalization = 'pdf';
hold on
h3 = histogram(log10(rho_screw),197);
h3.Normalization = 'pdf';
hold on
%plot(y,f,'LineWidth',1.5); %if it was a normal distribution
xlim([7, 15]);
xlabel('log10(\rho)')
legend('\rho_{totalGND}', '\rho_{edge}', '\rho_{screw}', 'Location','northwest');
hold off

%% Histograms - Misorientation angles
% materials (eg 'Copper' or 'indexed'
%remove small grain (optional)
ebsd(grains(grains.grainSize<5)) = [];
% repeat grain reconstruction
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'threshold',5*degree);
% smooth the grain boundaries a bit
grains = smooth(grains,5);
% get the misorientations to mean
mori = ebsd('Copper').mis2mean;
mori_mean = mean(mori);
% plot a histogram of the misorientation angles
plotAngleDistribution(mori);
xlabel('Misorientation angles in degree');

