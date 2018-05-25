%% Set-up parameters
clear all
clc

% Suppress figures
set(0,'DefaultFigureVisible', 'off')

normMethod = 'mean'; % options are 'mean' or 'dB' or 'median'

home = '/Volumes/MacintoshHD4/RFlex2';
addpath '/Users/eeg/Documents/MATLAB/fieldtrip-master';
ft_defaults

% Defining directories
datafolder = [home filesep 'ParentRFlex'];
parentfolder = [home filesep 'FFT_analysis'];
cd(parentfolder)

groupfolder = [parentfolder filesep 'Group'];

% Specify layout file for top plot
layout_file = [parentfolder filesep 'biosemi64.lay'];
% Define frequency resolution
% load EEG location file
load([parentfolder filesep 'EEGchanloc.mat']);
% Load in default TFR structure as a scaffold to plug in result data 
load([parentfolder filesep 'default_TFR_struct.mat'])
TFRres = [];
TFRres = TFRdef;

% Setting up subject loop
subject_list = textread([parentfolder filesep 'subjectlist_FFT.txt'], '%s');
 
freq_names = {'delta', 'theta', 'alpha', 'beta','lowGamma','medGamma'};

ctxt_labels = {'OO','FF','OF','FO'};

%% Show only significant results

% Clean up TFRres structure so it's appropriate for FFT analysis
TFRres.label = TFRres.label';
TFRres.freq = 1:length(freq_names);
TFRres.time = 1;


% Get mean across subjects: 64 (channel) by 5 (Frequency band) matrices
concat_condMeans = [];
for s=1:length(subject_list); 
    % Defining the 'subject' variable
    subject = [];
    subjectfolder = [];
    subject = subject_list{s}; 
    subjectfolder = [parentfolder filesep subject];
    % Input subject data
    if strcmp(normMethod,'mean') == 1
        load([subjectfolder filesep subject '_FFT_FreqBand_MeanCenter.mat']);
    elseif strcmp(normMethod,'median') == 1
        load([subjectfolder filesep subject '_FFT_FreqBand_MedianCenter.mat']);
    elseif strcmp(normMethod,'dB') == 1
        load([subjectfolder filesep subject '_FFT_FreqBand_DecibNorm.mat']);
    end
    
    % Get mean for each condition within subject
    
    sub_meanConditions = cell([1,4]);
    for ctxt = 1:length(ctxt_labels)
        currCtxt = ctxt_labels(ctxt);
        
        ctxt_concat = [];
        for b = 1:3
            blockData = [];
            blockData = outTableFreqBand{b,ctxt}{1};        
            ctxt_concat = cat(3,ctxt_concat, blockData);
        end
        % Get mean across instances within a condition
        ctxt_mean = [];
        ctxt_mean = mean(ctxt_concat,3);
        % Store these means in a new cell
        sub_meanConditions{1,ctxt} = num2cell(ctxt_mean);        
    end
    % Assign condition mean matrices across subjects
    concat_condMeans = vertcat(concat_condMeans,sub_meanConditions);    
end


% Get group means for each of the 4 conditions
groupData = cell([1,4]);
for ctxt = 1:length(ctxt_labels)
    groupMatrix_ctxt = [];
    for s = 1:length(subject_list)
        subData = [];
        subData = cell2mat(concat_condMeans{s,ctxt});
        groupMatrix_ctxt = cat(3,groupMatrix_ctxt,subData);
    end
    groupMean_ctxt = [] ;
    groupMean_ctxt = mean(groupMatrix_ctxt,3);
    groupData{1,ctxt} = groupMean_ctxt;
end

% Create differnce matrices based on specified contrasts: 
GainVal_meanDiff = [];
LossVal_meanDiff = [];
Inform_meanDiff = [];
GainVal_meanDiff = (groupData{1,1} + groupData{1,3}) - (groupData{1,2} + groupData{1,4});
LossVal_meanDiff = (groupData{1,1} + groupData{1,4}) - (groupData{1,2} + groupData{1,3});
Inform_meanDiff = (groupData{1,1} + groupData{1,2}) - (groupData{1,3} + groupData{1,4});
% Load in statistical maps
% Read in stat results
%load([groupfolder filesep 'TFCE_Stats' filesep 'TFCE_Pvals_context_stats_ME_LossVal.mat'])
%load([groupfolder filesep 'TFCE_Stats' filesep 'TFCE_Pvals_context_stats_ME_GainVal.mat'])
%load([groupfolder filesep 'TFCE_Stats' filesep 'TFCE_Pvals_context_stats_ME_Inform.mat'])

        

% Read in stat results
load([groupfolder filesep 'TFCE_Stats' filesep 'TFCE_LME_AcrossFreqBand_Results.mat'])

% Create statistics topographical plots per contrast
for band = 1:length(freq_names)
    currBandName = freq_names(band);
    
    % Pull out current pValues
    currPvals = [];
    currPvals = freqStats.(char(currBandName)).tfce_pvalue';
    
    % Pull out the corresponding frequency band from statstical map
    pVal_GainVal = [];
    pVal_LossVal = [];    
    pVal_Inform = [];
    pVal_GainVal = currPvals(:,2);
    pVal_LossVal =  currPvals(:,3);
    pVal_Inform =  currPvals(:,4);
    %pVal_LossVal = context_stats_TFCE_Pvals_ME_LossVal(:,band);
    %pVal_GainVal = context_stats_TFCE_Pvals_ME_GainVal(:,band);
    %pVal_Inform = context_stats_TFCE_Pvals_ME_Inform(:,band);
    
    % Zero-out p-values greater than 0.05 (returns boolean)
    pVal_LossVal = pVal_LossVal < 0.05;
    pVal_GainVal = pVal_GainVal < 0.05;
    pVal_Inform = pVal_Inform < 0.05;
    
    % Multiply masked p-values by corresponding difference map (for this
    % freq band)
    currLossVal_meanDiff = [];
    currGainVal_meanDiff = [];
    currInform_meanDiff = [];
    
    currLossVal_meanDiff = GainVal_meanDiff(:,band);
    currGainVal_meanDiff = LossVal_meanDiff(:,band);
    currInform_meanDiff = Inform_meanDiff(:,band);
    
    GainVal_meanDiff(:,band) = currLossVal_meanDiff.*pVal_GainVal(:,1);
    LossVal_meanDiff(:,band) = currGainVal_meanDiff.*pVal_LossVal(:,1);
    Inform_meanDiff(:,band) = currInform_meanDiff.*pVal_Inform(:,1);

    TFRres.powspctrm = [];
    TFRres.powspctrm = GainVal_meanDiff;
    plot_name = ['Sig. power differences for ' char(currBandName) ' frequencies: Gain Value']
    % Top plot
    cfg = [];
    cfg.ylim = [band-0.001 band+0.001];
    cfg.zlim = [-50 50];
    cfg.comment = 'no';
    cfg.colorbar = 'no';
    cfg.colormap = jet(100);
    cfg.style = 'both';
    cfg.layout = layout_file;
    f1 = figure;
    ft_topoplotER(cfg,TFRres);
    f1.Name = plot_name;
    saveas(f1,[parentfolder filesep 'Group' filesep 'TFCE_Stats' filesep 'Plot_GainValue_tfceLME_f' char(currBandName) '.png'],'png');
    
    TFRres.powspctrm = LossVal_meanDiff;
    plot_name = ['Sig. power differences for ' char(currBandName) ' frequencies: Loss Value']
    % Top plot
    cfg = [];
    cfg.ylim = [band-0.001 band+0.001];
    cfg.zlim = [-50 50];
    cfg.comment = 'no';
    cfg.colorbar = 'no';
    cfg.colormap = jet(100);
    cfg.style = 'both';
    %cfg.colorbar = 'SouthOutside';
    cfg.layout = layout_file;
    f2 = figure;
    ft_topoplotER(cfg,TFRres);
    f2.Name = plot_name;
    saveas(f2,[parentfolder filesep 'Group' filesep 'TFCE_Stats' filesep 'Plot_LossValue_tfceLME_f' char(currBandName) '.png'],'png');
    
	TFRres.powspctrm = Inform_meanDiff;
    plot_name = ['Sig. power differences for ' char(currBandName) ' frequencies: Symmetry']
    % Top plot
    cfg = [];
    cfg.ylim = [band-0.001 band+0.001];
    cfg.zlim = [-50 50];
    cfg.comment = 'no';
    cfg.colorbar = 'no';
    cfg.colormap = jet(100);
    cfg.style = 'both';
    cfg.layout = layout_file;
    f3 = figure;
    ft_topoplotER(cfg,TFRres);
    f3.Name = plot_name;
    saveas(f3,[parentfolder filesep 'Group' filesep 'TFCE_Stats' filesep 'Plot_Sym_tfceLME_f' char(currBandName) '.png'],'png');
end



    