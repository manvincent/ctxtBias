%% Statistics on time-frequency analyses

clear all
clc

% Repeated measures (2x2) ANOVA 
    % Gain stim vs. No Loss stim: ME STIM
    % Gain context vs. No Loss context: ME CTXT
    % Context-congruent vs. Incongruent: INTERAC

home = '/Volumes/MacintoshHD4/RFlex2';
addpath(genpath('/Users/eeg/Documents/MATLAB/eeglab13_4_4b'));
addpath '/Users/eeg/Documents/MATLAB/fieldtrip-master';

addpath(genpath('/Volumes/MacintoshHD4/RFlex2/FFT_analysis/ept_TFCE-matlab-master_lme'));
addpath(genpath('/Volumes/MacintoshHD4/RFlex2/FFT_analysis/ept_TFCE-matlab-master'));
ft_defaults

% Defining directories
parentfolder = [home filesep 'FFT_analysis'];
cd(parentfolder)
groupfolder = [parentfolder filesep 'Group'];

% Save/define output folder 
if ~exist([groupfolder filesep 'TFCE_Stats'],'dir');
    mkdir([groupfolder filesep 'TFCE_Stats']);
end
statsout = [groupfolder filesep 'TFCE_Stats'];

num_elec = 64; % Speicfy # of electodes in dataset
freq_names = {'delta', 'theta', 'alpha', 'beta', 'lowGamma','medGamma'};
load([parentfolder filesep 'EEGchanloc.mat']); %load in channel file (EEGLAB)

%% Run TFCE
freqStats = struct();
freqNull = struct();
for band = (1:4) % 1:length(freq_names)
    currBandName = [];
    currBandName = freq_names(band);
    % Load in pValues from lme
    statTable = [];
    load([groupfolder  filesep 'Freq_' char(currBandName) '_concatenatedData_mc.mat']);
    % Organize DVs and Design into separate tables
    statTable.Properties.VariableNames{'subj_idx'} = 'participant_id';
    dvTable =transpose(table2array(statTable(:,3:end-3)));
           
    factors = struct(...
        'name', {'Con1', 'Con2','Con3'}, ... % predictor name to permute
        'flag_within', {1,1,1}, ... ; % which predictors in that list are within-subject
        'num_perm', 1000);
        channel_neighbours = ept_ChN2(loc);
    
    model_description = [' ~ ', ...
                         'Con1 + Con2 + Con3 ', ...
                         '+ (1 | participant_id)']


    % run the model for all channels
    observed = [];
    perm_values = [];
    [observed, perm_values] = ept_slow_lme_permutation(...
        dvTable, ...
        statTable, ...
        model_description, ...
        factors, ...
        channel_neighbours);
    
    freqStats.(char(currBandName)) = observed
    freqNull.(char(currBandName)) = perm_values  
end

% Save result structures
parsave([statsout filesep 'TFCE_LME_AcrossFreqBand_Results.mat'],1,freqStats,'-v7.3')
parsave([statsout filesep 'TFCE_LME_AcrossFreqBand_NullDistrb.mat'],1,freqNull,'-v7.3')




