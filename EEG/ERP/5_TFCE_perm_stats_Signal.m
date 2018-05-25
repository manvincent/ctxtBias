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
addpath '/Volumes/MacintoshHD4/RFlex2/TF_analysis/resampling_statistical_toolkit';
addpath(genpath('/Volumes/MacintoshHD4/RFlex2/TF_analysis/ept_TFCE-matlab-master'));
ft_defaults

% Defining directories
parentfolder = [home filesep 'TF_analysis'];
cd(parentfolder)
groupfolder = [parentfolder filesep 'Group_Analysis_Signal'];

% Save/define output folder 
if ~exist([groupfolder filesep 'TFCE_Stats'],'dir');
    mkdir([groupfolder filesep 'TFCE_Stats']);
end
statsout = [groupfolder filesep 'TFCE_Stats'];

num_elec = 64; % Speicfy # of electodes in dataset
freq_range = 2:1:60; % Specify frequency range in TF decomp on current data
load([parentfolder filesep 'EEGchanloc.mat']); %load in channel file (EEGLAB)

% ANOVA parameters
nPerm      = 10000; % default number of permutations
%rSample    = 250; % default assumed sampling rate (not used in calculation but only for result viewer)
%saveName   = [statsout filesep 'ept_Results_', date, '.mat'];
type       = 'r'; % 'r' = repeated measures; 'm' = mixed ANOVA
nGroups 	= 1;

event_lock_list = {'stim'};
test_chan_list = {'all'};
%% Run Stats
curr_event = [];
curr_event = char(event_lock_list(1));


% Define condition-specifc compile data
FF_G = [];
FF_C = [];
FF_L = [];
FO_G = [];
FO_C = [];
FO_L = [];
OO_G = [];
OO_C = [];
OO_L = [];
OF_G = [];
OF_C = [];
OF_L = [];


load([groupfolder filesep 'Group_Signal_Cond_' curr_event '_Cell.mat']);    
FF_G = full_cond_array{1,1};
FF_C = full_cond_array{2,1};
FF_L = full_cond_array{3,1};
FO_G = full_cond_array{4,1};
FO_C = full_cond_array{5,1};
FO_L = full_cond_array{6,1};
OO_G = full_cond_array{7,1};
OO_C = full_cond_array{8,1};
OO_L = full_cond_array{9,1};
OF_G = full_cond_array{10,1};
OF_C = full_cond_array{11,1};
OF_L = full_cond_array{12,1};

% Set up 2x2 ANOVA matrix (User can modify this depending on their statistical question)
rep_ANOVA_cell = [];
rep_ANOVA_cell = cell(2,2);
% Set it up to be:
%               Gain Stim | No Loss Stim
% Gain Context      {1,1}       {1,2}
% No Loss Context   {2,1}       {2,2}

rep_ANOVA_cell{1,1} = squeeze(permute(OF_G(:,:,:,:),[4 1 2 3 ]));
rep_ANOVA_cell{1,2} = squeeze(permute(OF_L(:,:,:,:),[4 1 2 3 ]));
rep_ANOVA_cell{2,1} = squeeze(permute(FO_G(:,:,:,:),[4 1 2 3 ]));
rep_ANOVA_cell{2,2} = squeeze(permute(FO_L(:,:,:,:),[4 1 2 3 ]));

% Define the stats cell array to use for your test
stats_data = [];
stats_data = rep_ANOVA_cell;

% Run element-wise repeated measures ANOVA  test
[Info, Data, Results] = ept_TFCE_ANOVA(stats_data,loc)   ;

condition_stats_Fvals_ME_ctxt = Results.Obs.A;
condition_stats_Fvals_ME_stim = Results.Obs.B;
condition_stats_Fvals_INTERAC = Results.Obs.AB;

condition_stats_TFCE_Fvals_ME_ctxt = Results.TFCE_Obs.A;
condition_stats_TFCE_Fvals_ME_stim = Results.TFCE_Obs.B;
condition_stats_TFCE_Fvals_INTERAC = Results.TFCE_Obs.AB;

condition_stats_TFCE_Pvals_ME_ctxt = Results.P_Values.A;
condition_stats_TFCE_Pvals_ME_stim = Results.P_Values.B;
condition_stats_TFCE_Pvals_INTERAC = Results.P_Values.AB;

parsave([statsout filesep 'Fvals_allChan_' curr_event '_locked_ME_ctxt.mat'],1,condition_stats_Fvals_ME_ctxt,'-v7.3');
parsave([statsout filesep 'Fvals_allChan_' curr_event '_locked_ME_stim.mat'],1,condition_stats_Fvals_ME_stim,'-v7.3');
parsave([statsout filesep 'Fvals_allChan_' curr_event '_locked_INTERAC.mat'],1,condition_stats_Fvals_INTERAC,'-v7.3');             

parsave([statsout filesep 'TFCE_Fvals_allChan_' curr_event '_locked_ME_ctxt.mat'],1,condition_stats_TFCE_Fvals_ME_ctxt,'-v7.3');
parsave([statsout filesep 'TFCE_Fvals_allChan_' curr_event '_locked_ME_stim.mat'],1,condition_stats_TFCE_Fvals_ME_stim,'-v7.3');
parsave([statsout filesep 'TFCE_Fvals_allChan_' curr_event '_locked_INTERAC.mat'],1,condition_stats_TFCE_Fvals_INTERAC,'-v7.3');

parsave([statsout filesep 'TFCE_Pvals_allChan_' curr_event '_locked_ME_ctxt.mat'],1,condition_stats_TFCE_Pvals_ME_ctxt,'-v7.3');
parsave([statsout filesep 'TFCE_Pvals_allChan_' curr_event '_locked_ME_stim.mat'],1,condition_stats_TFCE_Pvals_ME_stim,'-v7.3');
parsave([statsout filesep 'TFCE_Pvals_allChan_' curr_event '_locked_INTERAC.mat'],1,condition_stats_TFCE_Pvals_INTERAC,'-v7.3');

