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


event_lock = 'stim'; % Options are 'stim' or 'resp'. Sets the 0 point for the trial-wise TF windows
temp_res = 0.02;
prestim_length = -0.2; % If event_lock = 'stim', time in seconds for the TF window before stimulus.
poststim_length = 0.5; % If event_lock = 'stim', time in seconds for the TF window after stimulus.
preresp_length = -0.5; % If event_lock = 'resp', time in seconds for the TF window before response.
postresp_length = 0.5; % If event_lock = 'resp', time in seconds for the TF window after response.

if strncmp(event_lock,'stim',4) == 1 
    condition_list_label = {'FF_G_stim_ind', 'FF_C_stim_ind', 'FF_L_stim_ind', 'FO_G_stim_ind', 'FO_C_stim_ind', 'FO_L_stim_ind', 'OO_G_stim_ind', 'OO_C_stim_ind', 'OO_L_stim_ind', 'OF_G_stim_ind', 'OF_C_stim_ind', 'OF_L_stim_ind'};
elseif strncmp(event_lock,'resp',4) == 1 
    condition_list_label = {'FF_G_resp_ind', 'FF_C_resp_ind', 'FF_L_resp_ind', 'FO_G_resp_ind', 'FO_C_resp_ind', 'FO_L_resp_ind', 'OO_G_resp_ind', 'OO_C_resp_ind', 'OO_L_resp_ind', 'OF_G_resp_ind', 'OF_C_resp_ind', 'OF_L_resp_ind'};
end
        
        
%% Extract mean TF matrices per condition for each subject
parfor s=1:length(subject_list); 


    % Defining the 'subject' variable
    subject = [];
    subjectfolder = [];

    subject = subject_list{s}; 
    subjectfolder = [parentfolder filesep subject];
    if ~exist(subjectfolder,'dir');
        mkdir(subjectfolder)
    end
    
    if ~exist([subjectfolder filesep 'Cond_Ave_TrialBN'],'dir');
        mkdir([subjectfolder filesep 'Cond_Ave_TrialBN'])
    end

    fprintf('\n\n\n*** TF trial sorting on subject %d (%s) ***\n\n\n', s, subject);


    % Load in subject's TFR data 
    TFRwave = [];;              
    TFRwave = load([subjectfolder filesep subject '_TFR_TrialBN_Gamma.mat']);
    TFRwave = TFRwave.TFRwave;

    % Load in subject's timing/task file
    complete_sub_cell_noNA = [];
    complete_sub_cell_noNA = load([subjectfolder filesep subject '_taskEEG_time_cell_noNA.mat']);
    complete_sub_cell_noNA = complete_sub_cell_noNA.complete_sub_cell_noNA;

    EEG_Stim_Time = [];
    EEG_Resp_Time = [];
    EEG_Stim_Time = double(cell2mat(complete_sub_cell_noNA(:,17)));
    EEG_Resp_Time = double(cell2mat(complete_sub_cell_noNA(:,18)));

    % Find the time indices for each condition 

    if strncmp(event_lock,'stim',4) == 1  
    OF_G_stim_ind = [];
    OF_C_stim_ind = [];
    OF_L_stim_ind = [];
    FF_G_stim_ind = [];
    FF_C_stim_ind = [];
    FF_L_stim_ind = [];
    FO_G_stim_ind = [];
    FO_C_stim_ind = [];
    FO_L_stim_ind = [];
    OO_G_stim_ind = [];
    OO_C_stim_ind = [];
    OO_L_stim_ind = [];
        for l = 1:length(EEG_Stim_Time)
            locations_stim = [];
            locations_resp = [];

            stim_time = [];
            if cell2mat(complete_sub_cell_noNA(l,16)) == 2 % Mask trials to include only even trials (out-of-sampe)
                if ~any(isnan(TFRwave.powspctrm(1,1,knnsearch(TFRwave.timeround',EEG_Stim_Time(l))))) == 1
                    if strncmp(complete_sub_cell_noNA(l,5),'OF',2) && strncmp(complete_sub_cell_noNA(l,13),'Gain  ',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        OF_G_stim_ind = [OF_G_stim_ind, stim_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OF',2) && strncmp(complete_sub_cell_noNA(l,13),'Cntrl ',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        OF_C_stim_ind = [OF_C_stim_ind, stim_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OF',2) && strncmp(complete_sub_cell_noNA(l,13),'NoLoss',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        OF_L_stim_ind = [OF_L_stim_ind, stim_time];   
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FF',2) && strncmp(complete_sub_cell_noNA(l,13),'Gain  ',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        FF_G_stim_ind = [FF_G_stim_ind, stim_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FF',2) && strncmp(complete_sub_cell_noNA(l,13),'Cntrl ',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        FF_C_stim_ind = [FF_C_stim_ind, stim_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FF',2) && strncmp(complete_sub_cell_noNA(l,13),'NoLoss',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        FF_L_stim_ind = [FF_L_stim_ind, stim_time];     
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FO',2) && strncmp(complete_sub_cell_noNA(l,13),'Gain  ',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        FO_G_stim_ind = [FO_G_stim_ind, stim_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FO',2) && strncmp(complete_sub_cell_noNA(l,13),'Cntrl ',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        FO_C_stim_ind = [FO_C_stim_ind, stim_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FO',2) && strncmp(complete_sub_cell_noNA(l,13),'NoLoss',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        FO_L_stim_ind = [FO_L_stim_ind, stim_time];     
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OO',2) && strncmp(complete_sub_cell_noNA(l,13),'Gain  ',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        OO_G_stim_ind = [OO_G_stim_ind, stim_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OO',2) && strncmp(complete_sub_cell_noNA(l,13),'Cntrl ',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        OO_C_stim_ind = [OO_C_stim_ind, stim_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OO',2) && strncmp(complete_sub_cell_noNA(l,13),'NoLoss',6) == 1 % OF Context - Gains
                        stim_time = EEG_Stim_Time(l);
                        OO_L_stim_ind = [OO_L_stim_ind, stim_time];                     
                    end
                end
            %end
        end

        % Define condition lists:
        condition_list = [];
        condition_list = {FF_G_stim_ind, FF_C_stim_ind, FF_L_stim_ind, FO_G_stim_ind, FO_C_stim_ind, FO_L_stim_ind, OO_G_stim_ind, OO_C_stim_ind, OO_L_stim_ind, OF_G_stim_ind, OF_C_stim_ind, OF_L_stim_ind};

    elseif strncmp(event_lock,'resp',4) == 1 
    OF_G_resp_ind = [];
    OF_C_resp_ind = [];
    OF_L_resp_ind = [];
    FF_G_resp_ind = [];
    FF_C_resp_ind = [];
    FF_L_resp_ind = [];
    FO_G_resp_ind = [];
    FO_C_resp_ind = [];
    FO_L_resp_ind = [];
    OO_G_resp_ind = [];
    OO_C_resp_ind = [];
    OO_L_resp_ind = [];
        for l = 1:length(EEG_Resp_Time)
            locations_resp = [];
            locations_resp = [];

            resp_time = [];
            if cell2mat(complete_sub_cell_noNA(l,16)) == 2
                if ~any(isnan(TFRwave.powspctrm(1,1,knnsearch(TFRwave.timeround',EEG_Resp_Time(l))))) == 1
                    if strncmp(complete_sub_cell_noNA(l,5),'OF',2) && strncmp(complete_sub_cell_noNA(l,13),'Gain  ',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        OF_G_resp_ind = [OF_G_resp_ind, resp_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OF',2) && strncmp(complete_sub_cell_noNA(l,13),'Cntrl ',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        OF_C_resp_ind = [OF_C_resp_ind, resp_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OF',2) && strncmp(complete_sub_cell_noNA(l,13),'NoLoss',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        OF_L_resp_ind = [OF_L_resp_ind, resp_time];   
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FF',2) && strncmp(complete_sub_cell_noNA(l,13),'Gain  ',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        FF_G_resp_ind = [FF_G_resp_ind, resp_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FF',2) && strncmp(complete_sub_cell_noNA(l,13),'Cntrl ',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        FF_C_resp_ind = [FF_C_resp_ind, resp_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FF',2) && strncmp(complete_sub_cell_noNA(l,13),'NoLoss',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        FF_L_resp_ind = [FF_L_resp_ind, resp_time];     
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FO',2) && strncmp(complete_sub_cell_noNA(l,13),'Gain  ',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        FO_G_resp_ind = [FO_G_resp_ind, resp_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FO',2) && strncmp(complete_sub_cell_noNA(l,13),'Cntrl ',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        FO_C_resp_ind = [FO_C_resp_ind, resp_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'FO',2) && strncmp(complete_sub_cell_noNA(l,13),'NoLoss',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        FO_L_resp_ind = [FO_L_resp_ind, resp_time];     
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OO',2) && strncmp(complete_sub_cell_noNA(l,13),'Gain  ',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        OO_G_resp_ind = [OO_G_resp_ind, resp_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OO',2) && strncmp(complete_sub_cell_noNA(l,13),'Cntrl ',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        OO_C_resp_ind = [OO_C_resp_ind, resp_time];
                    elseif strncmp(complete_sub_cell_noNA(l,5),'OO',2) && strncmp(complete_sub_cell_noNA(l,13),'NoLoss',6) == 1 % OF Context - Gains
                        resp_time = EEG_Resp_Time(l);
                        OO_L_resp_ind = [OO_L_resp_ind, resp_time];                     
                    end
                end
            %end
        end
        
        % Define condition lists:
        condition_list = [];
        condition_list = {FF_G_resp_ind, FF_C_resp_ind, FF_L_resp_ind, FO_G_resp_ind, FO_C_resp_ind, FO_L_resp_ind, OO_G_resp_ind, OO_C_resp_ind, OO_L_resp_ind, OF_G_resp_ind, OF_C_resp_ind, OF_L_resp_ind};

    end

    pre_trial = [];
    post_trial = [];
    if strncmp(event_lock,'stim',4) == 1  
        pre_trial = prestim_length;
        post_trial = poststim_length;            
    elseif strncmp(event_lock,'resp',4) == 1  
        pre_trial = preresp_length;
        post_trial = postresp_length;
    end

    TFRindex = struct('Start',[],'End',[]);
    TFRtrial = [];
    for c = 1:length(condition_list)
        curr_cond = [];
        curr_cond = condition_list{1,c};
        curr_cond_label = [];
        curr_cond_label = char(condition_list_label(c));

        if ~isempty(curr_cond) == 1
        % Mark the index in the TFR structure where each trial starts and ends, per condition
        curr_cond_start = [];
        curr_cond_start =  knnsearch(TFRwave.timeround',(pre_trial + curr_cond)');                
        TFRindex.Start.(curr_cond_label) = curr_cond_start;

        curr_cond_end = [];                
        curr_cond_end = knnsearch(TFRwave.timeround',(post_trial + curr_cond)'); 
        TFRindex.End.(curr_cond_label) = curr_cond_end;                

        cat_trials = [];
        for t = 1:length(curr_cond); 
            % Extract the TF information for each trial and store according to condition
            trial_name = [];
            trial_name = ['trial_' num2str(t)];
            TFRtrial.(curr_cond_label).(trial_name) = TFRwave.powspctrm(:,:,curr_cond_start(t):curr_cond_end(t));

            % Get the average across all trials within a condition
            trial_mat = [];
            trial_mat = TFRtrial.(curr_cond_label).(trial_name);

            cat_trials = cat(4,cat_trials,trial_mat);                                                                                         
        end
        
        TFRtrial.(curr_cond_label).mean = [];
        TFRtrial.(curr_cond_label).mean = mean(cat_trials,4);
        ave_cond = [];
        ave_cond = TFRtrial.(curr_cond_label).mean;   
        % Construct new fieldtrip-format TFR structure for saving
        TFRcondave  = [];
        TFRcondave.label = TFRwave.label;
        TFRcondave.dimord = TFRwave.dimord;
        TFRcondave.freq = TFRwave.freq;
        TFRcondave.time = prestim_length:temp_res:poststim_length;
        TFRcondave.powspctrm = ave_cond;
        TFRcondave.elec = TFRwave.elec;
        TFRcondave.cfg = TFRwave.cfg;
        TFRcondave.cond = curr_cond_label;      

        parsave([subjectfolder filesep 'Cond_Ave_TrialBN' filesep curr_cond_label '_meanTFR_' event_lock '.mat'],1,TFRcondave,'-v7.3');
        else
            continue
        end

    end
end

fprintf('\n Finished all \n');            
    
    