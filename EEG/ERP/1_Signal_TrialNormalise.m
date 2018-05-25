
%% Set-up parameters
clear all
clc




home = '/Volumes/MacintoshHD4/RFlex2';
addpath '/Users/eeg/Documents/MATLAB/fieldtrip-master';
addpath(genpath('/Users/eeg/Documents/MATLAB/eeglab13_4_4b'));
ft_defaults

% Defining directories
datafolder = [home filesep 'ParentRFlex'];
parentfolder = [home filesep 'TF_analysis'];
cd(parentfolder)

% Setting up subject loop
subject_list = textread([parentfolder filesep 'subjectlist_TF.txt'], '%s');


% Set up output temporal resolution
temp_res = 0.01; % Temporal resolution in seconds

% Defining your baseline normalization window, relative to locked-event
baseline_min_bound = [];
baseline_max_bound = [];
baseline_min_bound = -1.40;
baseline_max_bound = -0.20;
% Defining trial window, relative to locked-event
trial_min_bound = [];
trial_max_bound = [];
trial_min_bound = -0.18;
trial_max_bound = 0.50;

% Specify trigger codes for baseline normalisation and analyses

% Create acceptable EEG trigger arrays
OO_stim = [];
FF_stim = [];
OF_stim = [];
FO_stim = [];
FB = [];
OO_resp = [];
FF_resp = [];
OF_resp = [];
FO_resp = [];
stim_code_mat = [];
resp_code_mat = [];
    
OO_stim = [11,12,13];
FF_stim = [21,22,23];
OF_stim = [31,32,33];
FO_stim = [41,42,43];
FB = [19,29,39];
OO_resp = [54,55,56,7754,7755,7756];
FF_resp = [64,65,66,7764,7765,7766];
OF_resp = [74,75,76,7774,7775,7776];
FO_resp = [84,85,86,7784,7785,7786];

stim_code_mat = [OO_stim; FF_stim; OF_stim; FO_stim];
resp_code_mat = [OO_resp; FF_resp; OF_resp; FO_resp];

recycle('off')
        


%% Baseline normalisation
parfor s=1:length(subject_list); 


    % Defining the 'subject' variable
    subject = [];
    subjectfolder = [];

    subject = subject_list{s}; 
    subjectfolder = [parentfolder filesep subject];
    if ~exist(subjectfolder,'dir');
        mkdir(subjectfolder)
    end

    if exist([subjectfolder filesep subject '_Signal_TrialBN.mat'],'file');
        fprintf('\nData already baseline normalized for %s\n', subject);
    else

        try
        fprintf('\n\n\n*** Baseline Normalization for subject %d (%s) ***\n\n\n', s, subject);
        
        % Load in preprocessed data from EEGLAB 
        EEG = [];
        EEG = pop_loadset('filename',['preprocessed.' subject '.set'],'filepath', [datafolder filesep subject],'loadmode','all');
        EEG = eeg_checkset( EEG ); 

        % Convert EEGLAB structure to FieldTrip sturcture
        PreprocessedData = [];
        PreprocessedData = eeglab2fieldtrip(EEG,'preprocessing','none');
        
        % Round data to new temporal respolution
        PreprocessedData.timeround = [];
        PreprocessedData.timeround = round(PreprocessedData.time{1,1},2);

        % Load in subject's onset data
        Eventlist = [];
        Eventlist = load([home filesep 'ParentRFlex/Eventlist_Structures' filesep subject '_eventinfo.mat']);
        Eventlist = Eventlist.Eventlist;
        
        for e = 1:length(PreprocessedData.label) % For every electrode                        
            fprintf('\nComputing for subject %d (%s) electrode #%d %s \n',s, subject, e, PreprocessedData.label{e});

            % Create vector of all the locked-trigger (e.g. stim-lock) timings                                              
            locked_event_vec = [];
            for i = 1:numel(Eventlist)
                if sum(any(Eventlist(i).code == stim_code_mat)) > 0
                    eventtime = [];
                    event_time = round(Eventlist(i).time,3);
                    locked_event_vec = [locked_event_vec, event_time];                    
                end
            end
            
            % Extract average power of pre-lock baseline period for each trial
            for j = 1:length(locked_event_vec)
                min_basetime = [];
                max_basetime = [];
                min_basetime = double(locked_event_vec(j) + baseline_min_bound) ;
                max_basetime = double(locked_event_vec(j) + baseline_max_bound);

                min_trialtime = [];
                max_trialtime = [];
                min_trialtime = double(locked_event_vec(j) + trial_min_bound);
                max_trialtime = double(locked_event_vec(j) + trial_max_bound);

                min_basetime_round = [];
                max_basetime_round = [];
                min_trialtime_round = [];
                max_trialtime_round = [];


                % Find index in TFR time field of the baseline period for this trial
                trial_baseline_min_idx = [];
                trial_baseline_max_idx = [];
                trial_baseline_min_idx = find(ismember(PreprocessedData.timeround,min_basetime_round));
                trial_baseline_max_idx = find(ismember(PreprocessedData.timeround,max_basetime_round));


                % Extract mean absolute power for the baseline period for this trial
                ave_baseline_trial = [];
                ave_baseline_trial = mean(PreprocessedData.trial{1,1}(e,trial_baseline_min_idx:trial_baseline_max_idx));

                % Find index in TFR time field of the trial window for this trial
                trial_min_idx = [];
                trial_max_idx = [];
                trial_min_idx = find(ismember(PreprocessedData.timeround,min_trialtime_round));
                trial_max_idx = find(ismember(PreprocessedData.timeround,max_trialtime_round));

                % Compute baseline normalisation using the decibel transform
                PreprocessedData.trial{1,1}(e,trial_min_idx:trial_max_idx) = 10*log10( PreprocessedData.trial{1,1}(e,trial_min_idx:trial_max_idx) ./ ave_baseline_trial); 

            end                                                                                                                       

        end  % end electrode
        
        % Save baseline_normalised time frequency output to mat
        parsave([subjectfolder filesep subject '_Signal_TrialBN.mat'],1,PreprocessedData);
        fprintf('\n Saved Signal BN results for subject %d (%s)\n', s, subject);
        catch ME
        end
    end
end

%% Turn recycling back on
recycle('on'); 
      
             
      



 
