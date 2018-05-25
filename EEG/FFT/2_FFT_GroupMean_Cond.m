%% Set-up parameters
clear all
clc

home = '/Volumes/MacintoshHD4/RFlex2';
addpath '/Users/eeg/Documents/MATLAB/fieldtrip-master';
ft_defaults

% Defining directories
parentfolder = [home filesep 'FFT_analysis'];
cd(parentfolder)

% Setting up subject loop
subject_list = textread([parentfolder filesep 'subjectlist_FFT.txt'], '%s');


event_lock = 'Block'; % Options are 'stim' or 'resp' or 'fix'. Sets the 0 point for the trial-wise TF windows


% Set the labels of the channels to average
condition_list_label = {'OO', 'FF', 'OF', 'FO'};


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

        if exist([subjectfolder filesep char(subject) '_FFT_FreqBand_MeanCenter_freq.mat'],'file');            

            fprintf('\n\n\n*** Condition compiling for condition %d, subject %d (%s) ***\n\n\n', cond, s, subject);

            % Load in subject's FFT data
            FFTcondave = [];;              
            FFTcondave = load([subjectfolder filesep char(subject) '_FFT_FreqBand_MeanCenter_freq.mat']);
            FFTcondave = FFTcondave.outTableFreqBand;
            
            % Concatenate all instances of a condition 
            sub_data_mat = [];
            %sub_data_mat = cat(4, cell2mat(FFTcondave{1,curr_cond}),cell2mat(FFTcondave{2,curr_cond}),cell2mat(FFTcondave{3,curr_cond}));
            sub_data_mat = mean(cat(3,FFTcondave{:,curr_cond}{:}),3);

            % Compile all-channel condition matrices across subjects        
            full_subs_array = cat(3,full_subs_array,sub_data_mat);        

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
if ~exist([parentfolder filesep 'Group'],'dir');
    mkdir([parentfolder filesep 'Group']);
end
groupfolder = [parentfolder filesep 'Group'];

parsave([groupfolder filesep 'Group_Cond_' event_lock '_freq_Cell.mat'],1,full_cond_array,'-v7.3');

fprintf('\n Finished\n');   
