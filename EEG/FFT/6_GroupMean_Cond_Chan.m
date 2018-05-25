%% Set-up parameters
clear all
clc

home = '/Volumes/MacintoshHD4/RFlex2';
addpath '/Users/eeg/Documents/MATLAB/fieldtrip-master';
ft_defaults

% Defining directories
parentfolder = [home filesep 'TF_analysis'];
cd(parentfolder)

% Setting up subject loop
subject_list = textread([parentfolder filesep 'subjectlist_TF.txt'], '%s');


event_lock = 'stim'; % Options are 'stim' or 'resp'. Sets the 0 point for the trial-wise TF windows
chan_ave = 'TFCE_Res_Signal'; 
                   % Other options may be possible but need to be specified
                   % below, or ask me to set it up!
                   

% Set the labels of the channels to average
chan_labels = [];

if strncmp(chan_ave,'TFCE_Res_Signal',1) == 1 
    chan_labels = {'C1','C3','CP3','CP1','P1','P3','Pz','CPz','FC2','Cz','C2','C4','CP4','CP2','P2','P4'};
elseif strncmp(chan_ave,'FFT_Loc',1) == 1 
    chan_labels = {'FC3', 'FC1', 'C3', 'CP3', 'P1', 'P3', 'P5', 'Pz'};
end                   

if strncmp(event_lock,'stim',4) == 1 
    condition_list_label = {'FF_G_stim_ind', 'FF_C_stim_ind', 'FF_L_stim_ind', 'FO_G_stim_ind', 'FO_C_stim_ind', 'FO_L_stim_ind', 'OO_G_stim_ind', 'OO_C_stim_ind', 'OO_L_stim_ind', 'OF_G_stim_ind', 'OF_C_stim_ind', 'OF_L_stim_ind'};
elseif strncmp(event_lock,'resp',4) == 1 
    condition_list_label = {'FF_G_resp_ind', 'FF_C_resp_ind', 'FF_L_resp_ind', 'FO_G_resp_ind', 'FO_C_resp_ind', 'FO_L_resp_ind', 'OO_G_resp_ind', 'OO_C_resp_ind', 'OO_L_resp_ind', 'OF_G_resp_ind', 'OF_C_resp_ind', 'OF_L_resp_ind'};
end
        
%% Compute subject-specific channel-averaged condition averages 
full_cond_array = cell(length(condition_list_label),3);
full_cond_meanchan_array = cell(length(condition_list_label),3);
for cond = 1:length(condition_list_label);
    curr_cond = [];
    curr_cond = condition_list_label(cond);
    
    full_subs_array = [];
    full_subs_meanchan_array = [];
    for s=1:length(subject_list); 
        
       
        % Defining the 'subject' variable
        subject = [];
        subjectfolder = [];

        subject = subject_list{s}; 
        subjectfolder = [parentfolder filesep subject];
        if ~exist(subjectfolder,'dir');
            mkdir(subjectfolder)
        end

        chanavefolder = [];
        if ~exist([subjectfolder filesep 'Chan_Ave_TrialBN'],'dir');
            mkdir([subjectfolder filesep 'Chan_Ave_TrialBN']);
        end
        chanavefolder = [subjectfolder filesep 'Chan_Ave_TrialBN'];
                                       
        % Compile subject data into one cell array containing all conditions
        
        if exist([subjectfolder filesep 'Cond_Ave_TrialBN' filesep char(curr_cond) '_meanTFR_' event_lock '.mat'],'file');            

        fprintf('\n\n\n*** Condition compiling for condition %d, subject %d (%s) ***\n\n\n', cond, s, subject);

        % Load in subject's TFR data 
        TFRcondave = [];;              
        TFRcondave = load([subjectfolder filesep 'Cond_Ave_TrialBN' filesep char(curr_cond) '_meanTFR_' event_lock '.mat']);
        TFRcondave = TFRcondave.TFRcondave;
        
        % Define number of channels to average
        chan_no = [];
        for c = 1:length(chan_labels)
            if any(ismember(TFRcondave.label,chan_labels(c)))
                idx = [];
                idx = find(strcmp(TFRcondave.label, chan_labels(c)));
                chan_no = [chan_no, idx];
            else 
                fprintf('\n*** Check channel labels! ***\n')
                break
            end
        end
        
        % Average the 2D TF matrices across indicated channels
        select_chan_mat = [];
        for ch = 1:length(chan_no)
            curr_chan_TF = [];
            curr_chan_TF = squeeze(TFRcondave.powspctrm(chan_no(ch),:,:));
            fprintf('\n Compiling channel %d \n',chan_no(ch));
            select_chan_mat = cat(3,select_chan_mat, curr_chan_TF);
            
        end
        select_chan_mat = permute(select_chan_mat,[3 1 2]);
        select_chan_mean = [];
        select_chan_mean = mean(select_chan_mat,1);
        % Save these condition- and subject-specific chan_ave datasets
        parsave([chanavefolder filesep char(curr_cond) '_meanTFR_meanChan_' chan_ave '_' event_lock '.mat'],1,select_chan_mean,'-v7.3');
        
        % Compile channel-averaged condition matrices across subjects
        sub_data_meanchan_mat = [];
        sub_data_meanchan_mat = select_chan_mean;
        full_subs_meanchan_array = cat(4,full_subs_meanchan_array,sub_data_meanchan_mat);              
        
        % Compile all-channel condition matrices across subjects
        sub_data_mat = [];
        sub_data_mat = TFRcondave.powspctrm;
        full_subs_array = cat(4,full_subs_array,sub_data_mat);        
        
        else
            fprintf('\n*** Missing condition for %s : %s !!***\n',subject,char(curr_cond));
            full_cond_array{cond,3} = 'Warning! Missing subs with no trials in this condition!!'
            full_cond_meanchan_array{cond,3} = 'Warning! Missing subs with no trials in this condition!!'
            continue
        end
    end
         
    full_cond_array{cond,1} = full_subs_array;
    full_cond_array{cond,2} = curr_cond;
    
    full_cond_meanchan_array{cond,1} = full_subs_meanchan_array;
    full_cond_meanchan_array{cond,2} = curr_cond;
end

% Save output
if ~exist([parentfolder filesep 'Group_Analysis_TrialBN'],'dir');
    mkdir([parentfolder filesep 'Group_Analysis_TrialBN']);
end
groupfolder = [parentfolder filesep 'Group_Analysis_TrialBN'];

%parsave([groupfolder filesep 'Group_meanChan_all_Cond_' event_lock '_Cell.mat'],1,full_cond_array,'-v7.3');
parsave([groupfolder filesep 'Group_meanChan_' chan_ave '_Cond_' event_lock '_Cell.mat'],1,full_cond_meanchan_array,'-v7.3');

fprintf('\n Finished\n');   
