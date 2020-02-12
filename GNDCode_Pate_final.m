% --- INPUTS ---
%optional to autosave workspace results
%save_location = 'D:\Alienware\ResearchGroup\ORNL Project\Outputs\3-25-19\';
%save_filename = 'Dans-area1-1500x1500-5nmstep'
poissons = 0.355;

%% --- Import data with wizard the first time: import_wizard('ebsd') --- %%
%Copy and paste results for future below
%PASTE WIZARD DATA BELOW HERE - Copper included for reference
%% Specify Crystal and Specimen Symmetries
% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('432', [3.615 3.615 3.615], 'mineral', 'Copper', 'color', 'light blue')};
% plotting convention
setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','intoPlane');
%% Specify File Names
% path to files
pname = 'D:\Alienware\ResearchGroup\ORNL Project\Inputs\';
% which files to be imported
fname = [pname 'Cu-K1-03 pnt2ss_3250.ang'];
%% Import the Data
% create an EBSD variable containing the data
ebsd = loadEBSD(fname,CS,'interface','ang',...
  'convertEuler2SpatialReferenceFrame');
%PASTE WIZARD DATA ABOVE HERE
fprintf('Data loaded\n');
% %% If you need to crop the data
% region = [xmin ymin xmax-xmin ymax-ymin]*scale_units %(eg mm: 10E-3);
% plot(ebsd)
% rectangle('position',region,'edgecolor','r','linewidth',2) %verify with box it's correct
% condition = inpolygon(ebsd,region); %set outside rect to FALSE
% ebsd = ebsd(condition) %crop

%% BEGIN GND CODE %%
tic
[grains,ebsd.grainId] = ebsd.calcGrains;
%Filter/denoise to prevent over-estimating
fprintf('Denoising... \n');
f=halfQuadraticFilter; f.alpha = [0.01 0.01]; f.threshold = 5*degree; f.eps = 1e-2;
%ebsd=smooth(ebsd,f); 
ebsd = smooth(ebsd('indexed'),f,'fill',grains);

% --OPTIONAL: Remove all grains with less then 5 pixels and smooth the grain boundaries-- %%
fprintf('Optional: Removing small grains...')
ebsd(grains(grains.grainSize<=3)) = [];
[grains,ebsd.grainId] = calcGrains(ebsd);
% -- end OPTIONAL -- %

ebsd = ebsd.gridify; 
cs = ebsd('indexed').CS; 
a = norm(ebsd.CS.aAxis); % size of the unit cell (lattice parameter)

%% Dislocation systems- Preferred method: Use ebsd file data
dS = dislocationSystem.fcc(cs); % == CHANGE IF NOT FCC ===%
% --- OPTIONAL: Manually set screw and edge ratios (leave dS above)-- %
%Note: result in all cases depends on .u (dislocation energy) 
%which which is off the guess: e_edge/e_screw = 1/(1-poissonsratio) 
%or any realistic value
dS(dS.isEdge).u = 1; %E of edge dislocation (uncomment)
dS(dS.isScrew).u = 1 - poissons; %E of screw dislocations (uncomment)

%% Compute Curvature Tensor
kappa = ebsd.curvature;
%the curvature tensor is directly related to the dislocation density tensor: 
%alpha = kappa.dislocationDensity; %optional-to show relationship btwn a&k
% Rotate Dislocation System into Specimen Coordinates 
dSRot = ebsd.orientations * dS; 
%% Soving the fitting problem: fit disl. system to curvature tensor %%
fprintf('Begin fitting dislocation system to curvature tensor... \n');
[rho,factor] = fitDislocationSystems(kappa,dSRot); 
alpha = sum(dSRot.tensor .* rho,2);
alpha.opt.unit = '1/um'; %manually set units since not stored in rho
%we may also restore the complete curvature tensor with kappa = alpha.curvature;

%% DISLOCATOIN ENERGY calcs (units: 1/m^2)
% gnd total density per pixel is computed by
dislE_tot = factor*sum(abs(rho .* dSRot.u),2); 
%average gnd for area:
dislE_ave = nanmean(dislE_tot);
% gnd screw (13-18) 
dislE_screw = factor*sum(abs(rho .* dSRot.isScrew) .* (dSRot.u .* dSRot.isScrew),2); 
% gnd edge (1-12) 
dislE_edge = factor*sum(abs(rho .* dSRot.isEdge) .* (dSRot.u .* dSRot.isEdge),2); 

%% GND Density Method 1: from Rho directly (units: 1/m^2)
fprintf('Calculating GND Density with Method 1...\n');
rho_screw = factor*sum(abs(rho .* dSRot.isScrew),2);
rho_edge = factor*sum(abs(rho .* dSRot.isEdge),2);
rho_total = abs(rho_edge + rho_screw);
rho_ave = nanmean(rho_total);

%% MAKE PLOTS OF ALL THREE
setaxis = [13.5 15];
plot(ebsd,log10(rho_total)) %can also plot not on log scale
caxis(setaxis)
nextAxis 
plot(ebsd,log10(rho_screw)) 
caxis(setaxis)
nextAxis 
plot(ebsd,log10(rho_edge)) 
caxis(setaxis)
drawNow(gcm,'figSize','large')
mtexColorMap('jet')
mtexColorbar

%% Results based on Asher's Code / Ruggles eq 2.2 (units: 1/m^2) 
fprintf('Calculating L1-norm of alpha for Method 2...\n');
% Calculating L1 norm of alpha
alpha_norm = zeros(max(length(alpha)), 1);
for index = 1:max(length(alpha));
    val = matrix(alpha(index));
    alpha_norm(index) = norm(val,1);
end
fprintf('Calculating GND density with Method 2...\n');
% --- Will need  ADJUSTMETNS for non-FCC materials
% burgers vector magnitude as a function of LP - for BCC = sqrt(3)*LP/2, for FCC = LP/sqrt(2)
% a = lattice parameter (from ebsd data in code above)
VM = a/sqrt(2); 
gnd_method2 = (1/VM)*alpha_norm*factor;
gnd_method2_ave = nanmean(gnd_method2);
fprintf('Results:\n');
fprintf("Average GND for scan based on Asher's method: %d \n", nanmean(gnd_method2));
fprintf('Average GND for scan from gnd_tot (above): %d \n', nanmean(rho_total));

toc
%save(save_location);

%% PLOTS
%Alpha norm (method 2, like Asher's code) vs Rho (method 1)
%optional subtitles
axisrange = [10 15]; %change axis range here
sub1 = sprintf('Average GND density: %d', gnd_method2_ave);
sub2 = sprintf('Average GND density: %d', rho_ave);
%plots 
plotx2east; plotzIntoPlane;
plot(ebsd,log10(rho_total));
title({'Ave GND Density from |\alpha|_1'; ['\it\fontsize{12}' sub1]},'FontWeight','Normal');
set(gca,'CLim',[axisrange]);
nextAxis 
plot(ebsd,log10(gnd_method2));
title({'Ave GND Density from fitting \rho'; ['\it\fontsize{12}' sub2]},'FontWeight','Normal');
set(gca,'CLim',[axisrange]);
drawNow(gcm,'figSize','large');
mtexColorMap('jet');
mtexColorbar