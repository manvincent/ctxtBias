% EEGLAB Interpolation Script
% Uses EEGLAB v13.4.4
% Vincent Man
% October 17th 2015

clear all 
eeglab

rawfolder='/Volumes/MacintoshHD3/RFlex2_Final/RawRFlex'
subject='117591t'
interp_channel_no=[15, 16] % Look up channel label - number correspondance 

EEG.setname = [];
EEG = pop_readbdf( [rawfolder filesep subject '.bdf'],[],80,[65 66],'on' ); 
EEG.setname='BDF file';
EEG = eeg_checkset( EEG );

% Channel Locations
EEG=pop_chanedit(EEG, 'lookup','/Users/eeg/Documents/MATLAB/eeglab13_4_4b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
EEG = eeg_checkset( EEG );

EEG = pop_interp(EEG, interp_channel_no, 'spherical');
EEG = eeg_checkset( EEG );

EEG = pop_saveset( EEG, 'filename', subject ,'filepath', rawfolder);
eeglab redraw
