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
groupfolder = [parentfolder filesep 'Group_Analysis'];
    
    if ~exist([groupfolder filesep 'AllSub_Single_Trial_EEG'],'dir');
        mkdir([groupfolder filesep 'AllSub_Single_Trial_EEG']);
    end
    group_outputfolder = [groupfolder filesep 'AllSub_Single_Trial_EEG'];


% Setting up subject loop
subject_list = textread([parentfolder filesep 'subjectlist_TF.txt'], '%s');


event_lock = 'stim'; % Options are 'stim' or 'resp'. Sets the 0 point for the trial-wise TF windows

% Define the averaged frequency and time windows to extract for each trial
freq_win_min = 46; % In Hz
freq_win_max = 47;

temp_win_min_stim = 0.28; % If event_lock = 'stim', time in seconds for the TF window to extract
temp_win_max_stim = 0.30; 
%temp_win_min_resp = []; % If event_lock = 'resp', time in seconds for the TF window to extract
%temp_win_max_resp = []; 

chan_ave = 'TFCE_Res_Signal'; % Options are 'F' for frontal (CHAN) or 'P' for parietal (CHAN)
                   % Other options may be possible but need to be specified
                   % below, or ask me to set it up!
                   


% Set the labels of the channels to average
chan_labels = [];
if strncmp(chan_ave,'P',1) == 1 
    chan_labels = {'P2','P3','P5'};
elseif strncmp(chan_ave,'F',1) == 1 
    chan_labels = {'AF3','AFz','AF4','F1','Fz','F2','FC1','FCz','FC2'};
elseif strncmp(chan_ave,'FFT',1) == 1
    chan_labels = {'FC3', 'FC1', 'C3', 'CP3', 'P1', 'P3', 'P5', 'Pz'};
elseif strncmp(chan_ave,'TFCE_Res_Signal',1) == 1 
    chan_labels = {'C1','C3','CP3','CP1','P1','P3','Pz','CPz','FC2','Cz','C2','C4','CP4','CP2','P2','P4'};
end        

% Name the current extracted data
single_trial_name = 'flowGamma_TFCE_ResSignal_Stim'

%% Extract mean TF matrices per condition for each subject
all_sub_single_trial = [];
all_sub_single_trial_noNA = [];
for s=1:length(subject_list); 


    % Defining the 'subject' variable
    subject = [];
    subjectfolder = [];

    subject = subject_list{s}; 
    subjectfolder = [parentfolder filesep subject];
    if ~exist(subjectfolder,'dir');
        mkdir(subjectfolder)
    end
    
    if ~exist([subjectfolder filesep 'Single_Trial_EEG'],'dir');
        mkdir([subjectfolder filesep 'Single_Trial_EEG']);
    end
    outputfolder = [subjectfolder filesep 'Single_Trial_EEG'];

    fprintf('\n\n\n*** Single-trial data extraction for subject %d (%s) ***\n\n\n', s, subject);


    % Load in subject's TFR data 
    TFRwave = [];;              
    TFRwave = load([subjectfolder filesep subject '_TFR_BN']);
    TFRwave = TFRwave.TFRwave;

    % Load in subject's timing/task file
    complete_sub_table = [];
    complete_sub_table = load([subjectfolder filesep subject '_taskEEG_time.mat']);
    complete_sub_table = complete_sub_table.complete_sub_table;
      
    % Specify the TFR structure indices for the temporal windows per trial
    EEG_Stim_Time = [];
    EEG_Resp_Time = [];
    EEG_Stim_Time = double(table2array(complete_sub_table(:,17)));
    EEG_Resp_Time = double(table2array(complete_sub_table(:,18)));
             
    pre_wind = [];
    post_wind = [];
    time_vec = [];
    if strncmp(event_lock,'stim',4) == 1  
        time_vec = EEG_Stim_Time;
        pre_wind = temp_win_min_stim;
        post_wind = temp_win_max_stim;            
    elseif strncmp(event_lock,'resp',4) == 1  
        time_vec = EEG_Resp_Time;
        pre_wind = temp_win_min_resp;
        post_wind = temp_win_max_resp;
    end

    TFRindex = struct('Start',[],'End',[]);     
    % Mark the index in the TFR structure where each trial starts and ends
    wind_start_idx = [];
    wind_start_idx =  knnsearch(TFRwave.timeround',(pre_wind + time_vec));                
    TFRindex.Start.(event_lock) = wind_start_idx;

    wind_end_idx = [];                
    wind_end_idx = knnsearch(TFRwave.timeround',(post_wind + time_vec)); 
    TFRindex.End.(event_lock) = wind_end_idx;                

    % Extract TF matrices per trial 
    TFRtrial = struct(event_lock, []);
    TFRtrial.(event_lock) = struct('time',[],'freq',[], 'chan',[],'means',[]); 

    % Extract TF matrices per trial in specified temporal window, for all channels and frequencies
    cat_trials = [];
    for t = 1:length(time_vec);         
            trial_TF = [];
            if isnan(time_vec(t)) == 0
                if time_vec(t) <= max(TFRwave.timeround)       
                % Extract the TF information for each trial                 
                trial_TF = TFRwave.powspctrm(:,:,wind_start_idx(t):wind_end_idx(t));
                end
            elseif isnan(time_vec(t)) == 1
                % If missing stim time, create a trial TF matrix of the
                % same dimensions as the previous trials, filled with NaN
                trial_TF = NaN(size(TFRwave.powspctrm(:,:,wind_start_idx(t-1):wind_end_idx(t-1))),'like',TFRwave.powspctrm(:,:,wind_start_idx(t-1):wind_end_idx(t-1)));
            end                            
            cat_trials = cat(4,cat_trials,trial_TF);                                                             
    end
    
%     if size(cat_trials,4) < length(time_vec)
%         cat_trials = padarray(cat_trials,[0 0 0 length(time_vec)-size(cat_trials,4)],NaN,'post');
%     end
         
    TFRtrial.(event_lock).time = cat_trials;

    % Extract TF matrices per trial in specified frequency window, for all channels
    freq_range = TFRwave.freq ;  
    [diffvalMin,freq_min_idx] = min(abs(freq_range - freq_win_min));
    [diffvalMax,freq_max_idx] = min(abs(freq_range - freq_win_max));
    TFRtrial.(event_lock).freq =  TFRtrial.(event_lock).time(:,freq_min_idx:freq_max_idx,:,:);
   
    
    % Extract TF matrices per trial in specified channel range
    % Define chan_ID of channels to average
    chan_no = [];
    for c = 1:length(chan_labels)
        if any(ismember(TFRwave.label,chan_labels(c)))
            idx = [];
            idx = find(strcmp(TFRwave.label, chan_labels(c)));
            chan_no = [chan_no, idx];
        else 
            fprintf('\n*** Check channel labels! ***\n')
            break
        end
    end
    
    TFRtrial.(event_lock).chan = TFRtrial.(event_lock).freq(chan_no,:,:,:);
       
    % Average across channels, time, and frequency, for each trial
    trial_mats = [];
    single_trial_list = [];
    trial_mats = TFRtrial.(event_lock).chan;
    TFRtrial.(event_lock).means = squeeze(mean(mean(mean(trial_mats,3),2),1));
    single_trial_list = TFRtrial.(event_lock).means;
    
    
    % Concatenate task and EEG values for all trials     
    full_table = [] ;
    full_table = horzcat(complete_sub_table,array2table(single_trial_list));
    full_table.Properties.VariableNames{19} = char(single_trial_name);
    
    % Modify the data table according to data input for stats
    full_table.Context_Trial = [];
    full_table.Onset_task = [];
    full_table.CorrectResponse = [];            
    full_table.Properties.VariableNames{6} = char('rt');
    full_table.Properties.VariableNames{7} = char('response');   
    
    % Modify 'response' column
    temp_response = [];
    temp_response = full_table.response;
    for t = 1:length(temp_response)
        if strcmp(full_table.Response_Type(t,:),'Gain  ') == 1 
            temp_response(t,:) = 1;
        elseif strcmp(full_table.Response_Type(t,:),'NoLoss') == 1 
            temp_response(t,:) = -1;
        elseif strcmp(full_table.stim(t,:),'Cntrl ') == 1
            temp_response(t,:) = 0;
        end
    end    
    full_table.response = temp_response;
                   
    % Modify 'Response_Type' column
    temp_Response_Type = [];
    temp_Response_Type = full_table.Response_Type;                
    for t = 1:length(temp_Response_Type)
        if strcmp(temp_Response_Type(t,:),'NaN   ') == 1
            temp_Response_Type(t,:) = ' ';
        end
    end    
    full_table.Response_Type = temp_Response_Type;
    
    % Modify 'stim' column
    temp_stim = [];
    for t = 1:length(full_table.stim)
        if strcmp(full_table.stim(t,:),'Gain  ') == 1 
            temp_stim(t,:) = 1;
        elseif strcmp(full_table.stim(t,:),'NoLoss') == 1 
            temp_stim(t,:) = -1;
        elseif strcmp(full_table.stim(t,:),'Cntrl ') == 1
            temp_stim(t,:) = 0;
        end
    end
    full_table.stim = temp_stim;
    
    
    % Remove all odd trials (for the out-of-sample approach) 
    %temp_Data_Split = [];
    %temp_Data_Split = full_table.Data_Split;                
    %temp_Data_Split(temp_Data_Split == 1) = nan;  
    %full_table.Data_Split = temp_Data_Split;
    
    % Remove all NaN missing values        
    TF = ismissing(full_table);
    complete_sub_table_noNA = [];
    complete_sub_table_noNA = full_table(~any(TF,2),:);
    
    % Mean-center only noNA trials, per subject
    single_trial_noNA = [];
    single_trial_noNA_cent = [];
    single_trial_noNA = table2array(complete_sub_table_noNA(:,16));
    single_trial_noNA_cent = single_trial_noNA - mean(single_trial_noNA);
    
    complete_sub_table_noNA(:,16) = [];
    complete_sub_table_noNA(:,16) = array2table(single_trial_noNA_cent);
    complete_sub_table_noNA.Properties.VariableNames{16} = char('Ave_TF');
            
    % Output single-trial, mean-centered values for non-NA trials                   
    writetable(complete_sub_table_noNA,[outputfolder filesep subject '_' single_trial_name '_noNA.csv']);
    
    % Concatenate single-trial extracts across all subjects    
    all_sub_single_trial_noNA = vertcat(all_sub_single_trial_noNA,complete_sub_table_noNA);
       
end

% Output all-subject concatenated 
writetable(all_sub_single_trial_noNA,[group_outputfolder filesep 'Group_' single_trial_name '_noNA.csv']);


fprintf('\n Finished all \n');            
    
    