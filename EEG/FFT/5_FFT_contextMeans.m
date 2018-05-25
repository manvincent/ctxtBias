%% Set-up parameters
clear all
clc

home = '/Volumes/MacintoshHD4/RFlex2';
addpath '/Users/vym/Downloads/MathWorks/fieldtrip-master';
ft_defaults

% Defining directories
datafolder = [home filesep 'ParentRFlex'];
parentfolder = [home filesep 'FFT_analysis'];
cd(parentfolder)

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


freq_names = {'delta', 'theta', 'alpha', 'beta', 'lowGamma','medGamma'};
contrast_names = {'GainVal','LossVal','Inform'};
ctxt_labels = {'OO','FF','OF','FO'};


%% Create per-suject tables of sig chans per contrast at group level 

% Read in stat results
load([parentfolder filesep 'Group' filesep 'TFCE_Stats' filesep 'TFCE_LME_AcrossFreqBand_Results.mat'])

for band = 1:length(freq_names)
    currBandName = [];
    currBandName = freq_names(band);
     % Load in contrast p-value   
    fullStats = [];
    fullStats = freqStats.(char(currBandName)).tfce_pvalue';
    
   
    for contrast = 1:length(contrast_names)
        % Specify current contrast 
        currContrast = [];
        currContrast = contrast_names(contrast);

        % Load in concatenated data
        statTable = [];   
        load([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_concatenatedData_mc.mat']);
    
        % Extract p-values for this contrast, across all chanels
        resultsP = [];
        resultsP = fullStats(:,contrast+1);
        % Remove non-significant channels
        currConNonSigChanIdx = [];
        currConNonSigChanIdx = find(resultsP(:,1) > 0.05);
        statTable(:,2+currConNonSigChanIdx') = [];
        
        % Get mean of each channel for each context within subject 
        % Pull out column header
        inVar = [];
        groupVar = [];
        meanTable = [];
        
        inVar = statTable.Properties.VariableNames(1,3:end);
        groupVar = {'subj_idx','context'};
        meanTable = varfun(@mean,statTable,'InputVariables',inVar,'GroupingVariables',groupVar);
        
        writetable(meanTable,[parentfolder filesep 'Group' filesep 'ContextMeans_Contrast' char(currContrast) '_f' char(currBandName) '.csv'])
    end
end


    