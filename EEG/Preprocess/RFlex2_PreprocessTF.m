
clear all
clc

% Change working directory to RFlex EEG
cd('/Volumes/MacintoshHD3/RFlex2_Final/')

% Set pathways
addpath('/Users/eeg/Documents/MATLAB/eeglab11_0_5_4b');
addpath('/Volumes/MacintoshHD3/RFlex2_Final/ParentRFlex');
addpath('/Users/eeg/Documents/MATLAB/eeglab13_4_4b/plugins/bdfimport1.00');

% Defining directories
rawfolder = '/Volumes/MacintoshHD3/RFlex2_Final/RawRFlex/';
parentfolder = '/Volumes/MacintoshHD3/RFlex2_Final/ParentRFlex/';
if ~exist([parentfolder 'Eventlist_Structures/'],'dir');
    mkdir([parentfolder 'Eventlist_Structures/']);
end
eventlistfolder = [parentfolder 'Eventlist_Structures/'];


% Start EEGlab GUI
eeglab    

ALLERP = buildERPstruct([]);
CURRENTERP = 0; 

% Defining Filters
highpassfilterfq = 0.5
notchfilterfq = 60
% Define Sampling Rate
downsample_rate = 512
% Defining pre-ICA noise rejection voltage threshold
lowerlimitthres = -200
upperlimitthres = 200


% Setting up subject loop
subject_list_bdf = textread([parentfolder 'subjectlist_rflex2_bdf.txt'], '%s');
subject_list_set = textread([parentfolder 'subjectlist_rflex2_set.txt'], '%s');
subject_list_all = vertcat(subject_list_bdf, subject_list_set);


%% Preprocessing raw (bdf) files
for s=1:length(subject_list_all); 

ALLEEG = [];
ALLERP = [];
ERP = [];
EEG = [];

% Defining the 'subject' variable
subject = [];
subjectfolder = [];

subject = subject_list_all{s}; 
subjectfolder = [parentfolder subject '/'];

if exist([subjectfolder 'preprocessed.' subject '.set'],'file');
    fprintf('\nData already preprocessed for %s\n', subject);
else

try
    
fprintf('\n\n\n*** Processing subject %d (%s) ***\n\n\n', s, subject);
mkdir(subjectfolder)

% Loading in a dataset:

% Determine whether to load in .bdf or .set file. 1 if bdf, 0 if set
matches = [];
is_bdf = [];

matches = strfind(subject_list_bdf, subject);
is_bdf = any(vertcat(matches{:}));

% Load in the subject file
if is_bdf == 1
    EEG.setname = [];
    EEG = pop_readbdf( [rawfolder subject '.bdf'],[],80,[65 66],'on' );    
    EEG.setname=[subject '.BDF file'];
    EEG = eeg_checkset( EEG );
elseif is_bdf == 0
    EEG.setname = [];
    EEG = pop_loadset('filename',[subject '.set'],'filepath', rawfolder);
    EEG.setname=[subject '.SET file'];
    EEG = eeg_checkset( EEG ); 
end

% Channel Locations
EEG = pop_chanedit(EEG, 'lookup','/Users/eeg/Documents/MATLAB/eeglab13_4_4b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
EEG = eeg_checkset( EEG );

% Re-reference to Mastoids
EEG = pop_reref( EEG, [65 66] ); % This references the data to the mastoids and removes the two mastoid channels
EEG = eeg_checkset( EEG );   

% Delete extra channels
EEG = pop_select( EEG,'nochannel',{'EXG7' 'EXG8' 'GSR1' 'GSR2' 'Erg1' 'Erg2' 'Resp' 'Plet' 'Temp'}); % This removes the extra channels
EEG = eeg_checkset( EEG );

% Downsample
EEG = pop_resample(EEG, [downsample_rate]);
EEG = eeg_checkset( EEG );

% Apply High-pass filter (IIR Butterworth) 
EEG  = pop_basicfilter( EEG,  1:68 , 'Boundary', -99, 'Cutoff',  highpassfilterfq, 'Design', 'butter', 'Filter', 'highpass','RemoveDC', 'on' );
EEG.setname=['hipass.chan.reref.' subject '.BDF file'];
EEG = eeg_checkset( EEG );

% Apply Notch filter (IIR Butterworth) 
EEG  = pop_basicfilter( EEG,  1:68 , 'Boundary', -99, 'Cutoff',  notchfilterfq, 'Design', 'notch', 'Filter', 'PMnotch', 'RemoveDC', 'on' );
EEG.setname=['lowpass.hipass.chan.reref.' subject '.BDF file'];
EEG = eeg_checkset( EEG );

% Run ICA 
EEG = pop_runica(EEG, 'icatype','fastica', 'approach', 'symm','stabilization','on','maxNumIterations',500 );

% Run ADJUST automatic ICA-based cleaning!
cd([parentfolder  subject filesep]);
EEG = interface_ADJ(EEG,[subject '_ADJ_report.txt']);

% Deleting noise components
% Identifying the artifact components from the *_report.txt output of ADJUST
cd([parentfolder  subject filesep]);
[status, result] = system(['tail -n 2 ' subject '_ADJ_report.txt']);
if isstrprop(result(end), 'cntrl'), result(end) = []; end
result_num = str2num(result);

% Removing these artifact components
EEG = pop_subcomp( EEG, result_num, 0);

% Delete extra ocular channels
EEG = pop_select( EEG,'nochannel',{'LO1' 'LO2' 'IO1' 'IO2'});
EEG.setname=['nonOccular.cont.ICAcleaned.noiseEpochs.lowpass.hipass.chan.reref.' subject '.BDF file'];
EEG = eeg_checkset( EEG );

% Exporting Event files 
EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, 'Warning', 'off' );
EEG = eeg_checkset( EEG );

% Refer to function to modify triggers!!
cd(parentfolder);
EEG = retrigger(EEG, subject);

% Export and save eventlist .mat structures 
Eventlist = [];
Eventlist = EEG.EVENTLIST.eventinfo ;
save([eventlistfolder subject '_eventinfo.mat'],'Eventlist');

% Saving the cleaned, continuous ICA dataset
EEG = pop_saveset( EEG, 'filename',['preprocessed.' subject '.set'],'filepath', subjectfolder);
EEG = eeg_checkset( EEG );

% Return to home directory
cd('/Volumes/MacintoshHD3/RFlex2_Final/ParentRFlex')

catch ME
    
end

end
end

fprintf('\n\n\n**** FINISHED ****\n\n\n');
