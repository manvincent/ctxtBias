
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

% TF decomp options
wave_width = [];
freq_range = [];
pad_length = [];
temp_res = [];

wave_width = 5; % Width of wavelets (in # of cycles)
                % Lower width means greater temp res
                % Higher width means greater freq res
freq_range = 2:1:90; % min_freq : freq_step : max_freq
temp_res = 0.02; % Temporal resolution in seconds

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
        
%% Time frequency analysis

parfor s=1:length(subject_list); 


    % Defining the 'subject' variable
    subject = [];
    subjectfolder = [];

    subject = subject_list{s}; 
    subjectfolder = [parentfolder filesep subject];
    if ~exist(subjectfolder,'dir');
        mkdir(subjectfolder)
    end

    if exist([subjectfolder filesep subject '_TFR_Raw_Pow.mat'],'file');
        fprintf('\nData already analyzed for %s\n', subject);
    else

        try
        fprintf('\n\n\n*** TF analysis on subject %d (%s) ***\n\n\n', s, subject);
        % Load in preprocessed data from EEGLAB 
        EEG = [];
        EEG = pop_loadset('filename',['preprocessed.' subject '.set'],'filepath', [datafolder filesep subject],'loadmode','all');
        EEG = eeg_checkset( EEG ); 

        % Convert EEGLAB structure to FieldTrip sturcture
        PreprocessedData = [];
        PreprocessedData = eeglab2fieldtrip(EEG,'preprocessing','none')



        % Defining time frequency configuration file
        cfgTF              = [];        

        
        cfgTF.method       = 'wavelet';
        cfgTF.output       = 'pow';
        cfgTF.channel      = 'all';           
        cfgTF.width        = wave_width; 
        cfgTF.foi          = freq_range;	                
        cfgTF.toi          = min(PreprocessedData.time{1,1}):temp_res:max(PreprocessedData.time{1,1});

        % Run time frequency analysis
        TFRwave = ft_freqanalysis(cfgTF, PreprocessedData);

        % Save time frequency output to mat
        parsave([subjectfolder filesep subject '_TFR_Raw_Pow.mat'],1,TFRwave,'-v7.3');
        fprintf('\n Saved TF results for subject %d (%s)\n', s, subject);

        catch ME

        end
    end

end


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

    if exist([subjectfolder filesep subject '_TFR_TrialBN.mat'],'file');
        fprintf('\nData already baseline normalized for %s\n', subject);
    else

        try
        fprintf('\n\n\n*** Baseline Normalization for subject %d (%s) ***\n\n\n', s, subject);
        
        % Load in subject's TFR data 
        TFRwave = [];;              
        TFRwave = load([subjectfolder filesep subject '_TFR_Raw_Pow.mat']);
        TFRwave = TFRwave.TFRwave;

        % Load in subject's onset data
        Eventlist = [];
        Eventlist = load([home filesep 'ParentRFlex/Eventlist_Structures' filesep subject '_eventinfo.mat']);
        Eventlist = Eventlist.Eventlist;
        
        for e = 1:length(TFRwave.label) % For every electrode                        
            for f = 1:length(TFRwave.freq) % For every frequency
                fprintf('\nComputing for subject %d (%s) electrode #%d %s, freq = %d \n',s, subject, e, TFRwave.label{e},TFRwave.freq(f));
                
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
                    TFRwave.timeround = [];
                    % Rounding equation
                    if temp_res == 1
                        min_basetime_round = round(min_basetime,0);
                        max_basetime_round = round(max_basetime,0);
                        min_trialtime_round = round(min_trialtime,0);
                        max_trialtime_round = round(max_trialtime,0);
                        TFRwave.timeround = round(TFRwave.time,0);
                    elseif temp_res == 0.05
                        min_basetime_round = round(min_basetime*20)/20;
                        max_basetime_round = round(max_basetime*20)/20;
                        min_trialtime_round = round(min_trialtime*20)/20;
                        max_trialtime_round = round(max_trialtime*20)/20;
                        TFRwave.timeround = round(TFRwave.time*20)/20;
                    elseif temp_res == 0.02
                        min_basetime_round = round(min_basetime*50)/50;
                        max_basetime_round = round(max_basetime*50)/50;
                        min_trialtime_round = round(min_trialtime*50)/50;
                        max_trialtime_round = round(max_trialtime*50)/50;
                        TFRwave.timeround = round(TFRwave.time*50)/50;
                    else
                        fprintf('\nERROR! Check your temp_res and rounding equation\n');
                        break
                    end
                       
                    % Find index in TFR time field of the baseline period for this trial
                    trial_baseline_min_idx = [];
                    trial_baseline_max_idx = [];
                    trial_baseline_min_idx = find(ismember(TFRwave.timeround,min_basetime_round));
                    trial_baseline_max_idx = find(ismember(TFRwave.timeround,max_basetime_round));
                    
                    
                    % Extract mean absolute power for the baseline period for this trial
                    ave_baseline_trial = [];
                    ave_baseline_trial = mean(TFRwave.powspctrm(e,f,trial_baseline_min_idx:trial_baseline_max_idx));
                    
                    % Find index in TFR time field of the trial window for this trial
                    trial_min_idx = [];
                    trial_max_idx = [];
                    trial_min_idx = find(ismember(TFRwave.timeround,min_trialtime_round));
                    trial_max_idx = find(ismember(TFRwave.timeround,max_trialtime_round));

                    % Compute baseline normalisation using the decibel transform
                    TFRwave.powspctrm(e,f,trial_min_idx:trial_max_idx) = 10*log10( TFRwave.powspctrm(e,f,trial_min_idx:trial_max_idx) ./ ave_baseline_trial); 

                end % end trial                                                                                                                         
            end  % end frequency         
        end  % end electrode
        
        % Save baseline_normalised time frequency output to mat
        parsave([subjectfolder filesep subject '_TFR_TrialBN_Gamma.mat'],1,TFRwave);
        fprintf('\n Saved TF BN results for subject %d (%s)\n', s, subject);
        
        delete([subjectfolder filesep subject '_TFR_Raw_Pow.mat']);
        
        catch ME
        end
    end
end

%% Turn recycling back on
recycle('on'); 
      
             
      



 
