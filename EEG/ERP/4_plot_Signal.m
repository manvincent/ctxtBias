%% Plotting TF results

clear all
clc

home = '/Users/vincentman/Desktop'
addpath '/Users/eeg/Documents/MATLAB/fieldtrip-master';
addpath(genpath('/Users/eeg/Documents/MATLAB/othercolor'));
ft_defaults

% Defining directories
parentfolder = [home filesep 'TF_analysis'];
groupfolder = [parentfolder filesep 'Group_Analysis_Signal'];

% User-defined options
event_lock = 'stim' % options are stim or resp
test_chan = 'all' % 'P','F','all';

% Contrast setup
cell1 = 'OF_G';
cell2 = 'FO_L'; % Indicate the two cells to average for simple effects tests contrastt1

cell3 = 'OF_L';
cell4 = 'FO_G';% Indicate the two cells to average for simple effects tests contrast2

alpha = 0.05
interp_amnt = 2000; % How smooth you want the plots to be, higher is smoother (but takes longer)
plot_content = 'PMaskedDiff' % Options are: Pvals | MeanDiff | PMaskedDiff
plot_type = 'topTFR' % Options are topTFR (topographical) | singleTFR (single or mean channel)
stat_type = 'TFCE_Stats' % Options are: TFCE_Stats | Stats

stat_contrast = 'ME_stim';% Can be ME_stim | ME_ctxt | or INTERAC. Check contrast setup

test_chan = 'all'; % Only 'all' option

% Plot paramaters
plot_time_range = [-0.1:0.1:0.5];
cust_colormap = [165,0,38;215,48,39;244,109,67;253,174,97;254,224,144;255,255,191;224,243,248;171,217,233;116,173,209;69,117,180;49,54,149];
cust_colormap = flipud(cust_colormap/255); % Convert to matlab RGB 


plot_color = 'gray';
plot_title = 'Localize - Main Effect of Cue';

% Specify layout file for top plot
layout_file = [parentfolder filesep 'biosemi64.lay'];

% End user-defined options
    
statsfolder = [groupfolder filesep stat_type ];
cd(statsfolder)

% Save/define output folder 
if ~exist([statsfolder filesep 'Plots'],'dir');
    mkdir([statsfolder filesep 'Plots']);
end
plotout = [statsfolder filesep 'Plots'];
    
    
% Load mean maps
load([groupfolder filesep 'Group_Signal_Cond_' char(event_lock) '_Cell.mat']);  

data_array = [];
data_array = full_cond_array;

% Cell 1: OF_G
cell1_idx = [];
cell1_idx = find(strcmp([data_array{:,2}], [cell1 '_' event_lock '_ind']));
cell1 = data_array{cell1_idx,1};
% Cell 2: FO_L
cell2_idx = [];
cell2_idx = find(strcmp([data_array{:,2}], [cell2 '_' event_lock '_ind']));
cell2 = data_array{cell2_idx,1};

% Cell 3: OF_L 
cell3_idx = [];
cell3_idx = find(strcmp([data_array{:,2}], [cell3 '_' event_lock '_ind']));
cell3 = data_array{cell3_idx,1};
% Cell 4: FO_G
cell4_idx = [];
cell4_idx = find(strcmp([data_array{:,2}], [cell4 '_' event_lock '_ind']));
cell4 = data_array{cell4_idx,1};



group_contrast = [];
if strcmp(stat_contrast,'ME_ctxt') == 1
    % Compute mean between the two cells for contrast 1
    cell13_ave = [];
    cell13_cat = [];
    cell13_groupave = [];
    cell13_cat = cat(5,cell1,cell3);
    cell13_ave = mean(cell13_cat,5);
    cell13_groupave = mean(cell13_ave,4);
    
    % Compute mean between the two cells for contrast 2
    cell24_ave = [];
    cell24_cat = [];
    cell24_groupave = [];
    cell24_cat = cat(5,cell2,cell4);
    cell24_ave = mean(cell24_cat,5);
    cell24_groupave = mean(cell24_ave,4);
    % Get contrast (cell1 + cell3) - (cell2 + cell4)
    % This (context-stim) congruent vs incongruent
    group_contrast = cell13_groupave - cell24_groupave; 
elseif strcmp(stat_contrast,'ME_stim') == 1
    % Compute mean between the two cells for contrast 1
    cell14_ave = [];
    cell14_cat = [];
    cell14_groupave = [];
    cell14_cat = cat(5,cell1,cell4);
    cell14_ave = mean(cell14_cat,5);
    cell14_groupave = mean(cell14_ave,4);
    
    % Compute mean between the two cells for contrast 2
    cell23_ave = [];
    cell23_cat = [];
    cell23_groupave = [];
    cell23_cat = cat(5,cell2,cell3);
    cell23_ave = mean(cell23_cat,5);
    cell23_groupave = mean(cell23_ave,4);
    % Get contrast (cell1 + cell4) - (cell2 + cell3)
    % This (context-stim) congruent vs incongruent
    group_contrast = cell14_groupave - cell23_groupave;
elseif strcmp(stat_contrast,'INTERAC') == 1
    % Compute mean between the two cells for contrast 1
    cell12_ave = [];
    cell12_cat = [];
    cell12_groupave = [];
    cell12_cat = cat(5,cell1,cell2);
    cell12_ave = mean(cell12_cat,5);
    cell12_groupave = mean(cell12_ave,4);
    
    % Compute mean between the two cells for contrast 2
    cell34_ave = [];
    cell34_cat = [];
    cell34_groupave = [];
    cell34_cat = cat(5,cell3,cell4);
    cell34_ave = mean(cell34_cat,5);
    cell34_groupave = mean(cell34_ave,4);
    % Get contrast (cell1 + cell2) - (cell3 + cell4)
    % This (context-stim) congruent vs incongruent
    group_contrast = cell12_groupave - cell34_groupave;        
end

% Load in p-vals matrix
res_data = [];
field = [];
res_data = load([statsfolder filesep 'TFCE_Pvals_allChan_' event_lock '_locked_' stat_contrast '.mat']);
field = fieldnames(res_data);
pvals_map = [];
pvals_bin = [];
pvals_map = double(getfield(res_data,field{1}));
pvals_map(pvals_map > alpha) = 0;
pvals_map(pvals_map > 0) = 1;

% Mask cell-mean matrix by pvals matrix
cell_groupave_masked = zeros(size(group_contrast));
cell_groupave_masked(pvals_map>0) = group_contrast(pvals_map>0);

% Load in default TFR structure as a scaffold to plug in result data 
load([parentfolder filesep 'default_TFR_struct.mat'])

TFRres = [];
TFRres = TFRdef;
TFRres.powspctrm = permute(repmat(double(cell_groupave_masked),1,1,size(TFRres.freq,2)),[1 3 2]);
TFRres.time = linspace(-.2,.5,size(TFRres.powspctrm,3));

% Top plot
cfg = [];
cfg.xlim = plot_time_range;
cfg.zlim = [-.6,.6];
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.colorbar = 'southoutside';
cfg.colormap = cust_colormap
cfg.style = 'straight';
cfg.interactive = plot_color;
cfg.layout = layout_file;
f = figure;
ft_topoplotTFR(cfg,TFRres);
f.Name = plot_title;

