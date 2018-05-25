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


event_lock = 'stim'; % Options are 'stim' or 'resp' or 'fix'. Sets the 0 point for the trial-wise TF windows


% Set the labels of the channels to average
if strncmp(event_lock,'stim',4) == 1 
    condition_list_label = {'FF_G_stim_ind', 'FF_C_stim_ind', 'FF_L_stim_ind', 'FO_G_stim_ind', 'FO_C_stim_ind', 'FO_L_stim_ind', 'OO_G_stim_ind', 'OO_C_stim_ind', 'OO_L_stim_ind', 'OF_G_stim_ind', 'OF_C_stim_ind', 'OF_L_stim_ind'};
elseif strncmp(event_lock,'resp',4) == 1 
    condition_list_label = {'FF_G_resp_ind', 'FF_C_resp_ind', 'FF_L_resp_ind', 'FO_G_resp_ind', 'FO_C_resp_ind', 'FO_L_resp_ind', 'OO_G_resp_ind', 'OO_C_resp_ind', 'OO_L_resp_ind', 'OF_G_resp_ind', 'OF_C_resp_ind', 'OF_L_resp_ind'};
elseif strncmp(event_lock,'fix',3) == 1 
    condition_list_label = {'FF_G_fix_ind', 'FF_C_fix_ind', 'FF_L_fix_ind', 'FO_G_fix_ind', 'FO_C_fix_ind', 'FO_L_fix_ind', 'OO_G_fix_ind', 'OO_C_fix_ind', 'OO_L_fix_ind', 'OF_G_fix_ind', 'OF_C_fix_ind', 'OF_L_fix_ind'};
end
        
%% Average conditions across subjects
full_cond_array = cell(length(condition_list_label),3);
for cond = 1:length(condition_list_label);
    curr_cond = [];
    curr_cond = condition_list_label(cond);
    
    full_subs_array = [];
    for s=1:length(subject_list); 
        
                       
        % Defining the 'subject' variable
        subject = [];
        subjectfolder = [];

        subject = subject_list{s};         
        subjectfolder = [parentfolder filesep subject];
        if ~exist(subjectfolder,'dir');
            mkdir(subjectfolder)
        end

        % Compile subject data into one cell array containing all conditions

        if exist([subjectfolder filesep 'Cond_Ave_TrialBN' filesep char(curr_cond) '_meanTFR_' event_lock '.mat'],'file');            

            fprintf('\n\n\n*** Condition compiling for condition %d, subject %d (%s) ***\n\n\n', cond, s, subject);

            % Load in subject's TFR data 
            TFRcondave = [];;              
            TFRcondave = load([subjectfolder filesep 'Cond_Ave_TrialBN' filesep char(curr_cond) '_meanTFR_' event_lock '.mat']);
            TFRcondave = TFRcondave.TFRcondave;


            % Compile all-channel condition matrices across subjects
            sub_data_mat = [];
            sub_data_mat = TFRcondave.powspctrm;
            full_subs_array = cat(4,full_subs_array,sub_data_mat);        

        else
            fprintf('\n*** Missing condition for %s : %s !!***\n',subject,char(curr_cond));
            full_cond_array{cond,3} = subject;            
            continue
        end
    end
         
    full_cond_array{cond,1} = full_subs_array;
    full_cond_array{cond,2} = curr_cond;
    
end

% Save output
if ~exist([parentfolder filesep 'Group_Analysis_TrialBN'],'dir');
    mkdir([parentfolder filesep 'Group_Analysis_TrialBN']);
end
groupfolder = [parentfolder filesep 'Group_Analysis_TrialBN'];

parsave([groupfolder filesep 'Group_Cond_' event_lock '_Cell.mat'],1,full_cond_array,'-v7.3');

fprintf('\n Finished\n');   
