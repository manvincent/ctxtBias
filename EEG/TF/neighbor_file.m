%% Create neighbour file for Fieldtrip cluster permutation analysis


home = '/Volumes/MacintoshHD3/RFlex2_Final';
addpath '/Users/eeg/Documents/MATLAB/fieldtrip-master';
addpath(genpath('/Users/eeg/Documents/MATLAB/eeglab13_4_4b'));
ft_defaults

% Defining directories
datafolder = [home filesep 'ParentRFlex'];
parentfolder = [home filesep 'TF_analysis'];
cd(parentfolder)


% Use EEGLAB to create a STUDY structure
STUDY = [];
ALLEEG = [];
[STUDY ALLEEG] = std_editset( STUDY, [], 'commands', { ...
	{ 'index' 1 'load' [datafolder filesep '13071t' filesep 'preprocessed.13071t.set'] 'subject' '13071t' 'condition' 1 'group' 1 }, ...
	{ 'index' 2 'load' [datafolder filesep '14070t' filesep 'preprocessed.14070t.set'] 'subject' '14070t' 'condition' 2 'group' 2 }, ...
	{ 'index' 3 'load' [datafolder filesep '17091t' filesep 'preprocessed.17091t.set'] 'subject' '17091t' 'condition' 3 'group' 3 }, ...
	{ 'index' 4 'load' [datafolder filesep '19101t' filesep 'preprocessed.19101t.set'] 'subject' '19101t' 'condition' 4 'group' 4 }, ...
	{ 'index' 5 'load' [datafolder filesep '20100t' filesep 'preprocessed.20100t.set'] 'subject' '20100t' 'condition' 5 'group' 5 }, ...
	{ 'index' 6 'load' [datafolder filesep '22110t' filesep 'preprocessed.22110t.set'] 'subject' '22110t' 'condition' 6 'group' 6 }, ...
	{ 'index' 7 'load' [datafolder filesep '23121t' filesep 'preprocessed.23121t.set'] 'subject' '23121t' 'condition' 7 'group' 7 }, ...
	{ 'index' 8 'load' [datafolder filesep '24120t' filesep 'preprocessed.24120t.set'] 'subject' '24120t' 'condition' 8 'group' 8 }, ...
	{ 'index' 9 'load' [datafolder filesep '25130t' filesep 'preprocessed.25130t.set'] 'subject' '25130t' 'condition' 9 'group' 9 }, ...
	{ 'index' 10 'load' [datafolder filesep '26131t' filesep 'preprocessed.26131t.set'] 'subject' '26131t' 'condition' 10 'group' 10 }, ...
	{ 'index' 11 'load' [datafolder filesep '27140t' filesep 'preprocessed.27140t.set'] 'subject' '27140t' 'condition' 11 'group' 11 }, ...
	{ 'index' 12 'load' [datafolder filesep '28141t' filesep 'preprocessed.28141t.set'] 'subject' '28141t' 'condition' 12 'group' 12 }, ...
	{ 'index' 13 'load' [datafolder filesep '29150t' filesep 'preprocessed.29150t.set'] 'subject' '29150t' 'condition' 13 'group' 13 }, ...
	{ 'index' 14 'load' [datafolder filesep '33170t' filesep 'preprocessed.33170t.set'] 'subject' '33170t' 'condition' 14 'group' 14 }, ...
	{ 'index' 15 'load' [datafolder filesep '39201t' filesep 'preprocessed.39201t.set'] 'subject' '39201t' 'condition' 15 'group' 15 }, ...
	{ 'index' 16 'load' [datafolder filesep '44220t' filesep 'preprocessed.44220t.set'] 'subject' '44220t' 'condition' 16 'group' 16 }, ...
	{ 'index' 17 'load' [datafolder filesep '45231t' filesep 'preprocessed.45231t.set'] 'subject' '45231t' 'condition' 17 'group' 17 }, ...
	{ 'index' 18 'load' [datafolder filesep '46230t' filesep 'preprocessed.46230t.set'] 'subject' '46230t' 'condition' 18 'group' 18 }, ...
	{ 'index' 19 'load' [datafolder filesep '48240t' filesep 'preprocessed.48240t.set'] 'subject' '48240t' 'condition' 19 'group' 19 }, ...
	{ 'index' 20 'load' [datafolder filesep '49250t' filesep 'preprocessed.49250t.set'] 'subject' '49250t' 'condition' 20 'group' 20 }, ...
	{ 'index' 21 'load' [datafolder filesep '50251t' filesep 'preprocessed.50251t.set'] 'subject' '50251t' 'condition' 21 'group' 21 }, ...
	{ 'index' 22 'load' [datafolder filesep '52261t' filesep 'preprocessed.52261t.set'] 'subject' '52261t' 'condition' 22 'group' 22 }, ...
	{ 'index' 23 'load' [datafolder filesep '53270t' filesep 'preprocessed.53270t.set'] 'subject' '53270t' 'condition' 23 'group' 23 }, ...
	{ 'index' 24 'load' [datafolder filesep '54271t' filesep 'preprocessed.54271t.set'] 'subject' '54271t' 'condition' 24 'group' 24 }, ...
	{ 'index' 25 'load' [datafolder filesep '56281t' filesep 'preprocessed.56281t.set'] 'subject' '56281t' 'condition' 25 'group' 25 }, ...
	{ 'index' 26 'load' [datafolder filesep '58291t' filesep 'preprocessed.58291t.set'] 'subject' '58291t' 'condition' 26 'group' 26 }, ...
	{ 'index' 27 'load' [datafolder filesep '59300t' filesep 'preprocessed.59300t.set'] 'subject' '59300t' 'condition' 27 'group' 27 }, ...
	{ 'index' 28 'load' [datafolder filesep '60301t' filesep 'preprocessed.60301t.set'] 'subject' '60301t' 'condition' 28 'group' 28 }, ...
	{ 'index' 29 'load' [datafolder filesep '62310t' filesep 'preprocessed.62310t.set'] 'subject' '62310t' 'condition' 29 'group' 29 }, ...
	{ 'index' 30 'load' [datafolder filesep '63321t' filesep 'preprocessed.63321t.set'] 'subject' '63321t' 'condition' 30 'group' 30 }, ...
	{ 'index' 31 'load' [datafolder filesep '65331t' filesep 'preprocessed.65331t.set'] 'subject' '65331t' 'condition' 31 'group' 31 }, ...
	{ 'index' 32 'load' [datafolder filesep '66330t' filesep 'preprocessed.66330t.set'] 'subject' '66330t' 'condition' 32 'group' 32 }, ...
	{ 'index' 33 'load' [datafolder filesep '69351t' filesep 'preprocessed.69351t.set'] 'subject' '69351t' 'condition' 33 'group' 33 }, ...
	{ 'index' 34 'load' [datafolder filesep '70350t' filesep 'preprocessed.70350t.set'] 'subject' '70350t' 'condition' 34 'group' 34 }, ...
	{ 'index' 35 'load' [datafolder filesep '71361t' filesep 'preprocessed.71361t.set'] 'subject' '71361t' 'condition' 35 'group' 35 }, ...
	{ 'index' 36 'load' [datafolder filesep '74371t' filesep 'preprocessed.74371t.set'] 'subject' '74371t' 'condition' 36 'group' 36 }, ...
	{ 'index' 37 'load' [datafolder filesep '77391t' filesep 'preprocessed.77391t.set'] 'subject' '77391t' 'condition' 37 'group' 37 }, ...
	{ 'index' 38 'load' [datafolder filesep '78390t' filesep 'preprocessed.78390t.set'] 'subject' '78390t' 'condition' 38 'group' 38 }, ...
	{ 'index' 39 'load' [datafolder filesep '79400t' filesep 'preprocessed.79400t.set'] 'subject' '79400t' 'condition' 39 'group' 39 }, ...
	{ 'index' 40 'load' [datafolder filesep '80401t' filesep 'preprocessed.80401t.set'] 'subject' '80401t' 'condition' 40 'group' 40 }, ...
	{ 'index' 41 'load' [datafolder filesep '81410t' filesep 'preprocessed.81410t.set'] 'subject' '81410t' 'condition' 41 'group' 41 }, ...
	{ 'index' 42 'load' [datafolder filesep '82411t' filesep 'preprocessed.82411t.set'] 'subject' '82411t' 'condition' 42 'group' 42 }, ...
	{ 'index' 43 'load' [datafolder filesep '85431t' filesep 'preprocessed.85431t.set'] 'subject' '85431t' 'condition' 43 'group' 43 }, ...
	{ 'index' 44 'load' [datafolder filesep '87441t' filesep 'preprocessed.87441t.set'] 'subject' '87441t' 'condition' 44 'group' 44 }, ...
	{ 'index' 45 'load' [datafolder filesep '89451t' filesep 'preprocessed.89451t.set'] 'subject' '89451t' 'condition' 45 'group' 45 }, ...
	{ 'index' 46 'load' [datafolder filesep '90450t' filesep 'preprocessed.90450t.set'] 'subject' '90450t' 'condition' 46 'group' 46 }, ...
	{ 'index' 47 'load' [datafolder filesep '91461t' filesep 'preprocessed.91461t.set'] 'subject' '91461t' 'condition' 47 'group' 47 }, ...
	{ 'index' 48 'load' [datafolder filesep '92460t' filesep 'preprocessed.92460t.set'] 'subject' '92460t' 'condition' 48 'group' 48 }, ...
	{ 'index' 49 'load' [datafolder filesep '93471t' filesep 'preprocessed.93471t.set'] 'subject' '93471t' 'condition' 49 'group' 49 }, ...
	{ 'index' 50 'load' [datafolder filesep '94470t' filesep 'preprocessed.94470t.set'] 'subject' '94470t' 'condition' 50 'group' 50 }, ...
	{ 'index' 51 'load' [datafolder filesep '95481t' filesep 'preprocessed.95481t.set'] 'subject' '95481t' 'condition' 51 'group' 51 }, ...
	{ 'index' 52 'load' [datafolder filesep '96480t' filesep 'preprocessed.96480t.set'] 'subject' '96480t' 'condition' 52 'group' 52 }, ...
	{ 'index' 53 'load' [datafolder filesep '97490t' filesep 'preprocessed.97490t.set'] 'subject' '97490t' 'condition' 53 'group' 53 }, ...
	{ 'index' 54 'load' [datafolder filesep '98491t' filesep 'preprocessed.98491t.set'] 'subject' '98491t' 'condition' 54 'group' 54 }, ...
	{ 'index' 55 'load' [datafolder filesep '99500t' filesep 'preprocessed.99500t.set'] 'subject' '99500t' 'condition' 55 'group' 55 }, ...
	{ 'index' 56 'load' [datafolder filesep '100501t' filesep 'preprocessed.100501t.set'] 'subject' '100501t' 'condition' 56 'group' 56 }, ...
	{ 'index' 57 'load' [datafolder filesep '101510t' filesep 'preprocessed.101510t.set'] 'subject' '101510t' 'condition' 57 'group' 57 }, ...
	{ 'index' 58 'load' [datafolder filesep '103520t' filesep 'preprocessed.103520t.set'] 'subject' '103520t' 'condition' 58 'group' 58 }, ...
	{ 'index' 59 'load' [datafolder filesep '104521t' filesep 'preprocessed.104521t.set'] 'subject' '104521t' 'condition' 59 'group' 59 }, ...
	{ 'index' 60 'load' [datafolder filesep '105530t' filesep 'preprocessed.105530t.set'] 'subject' '105530t' 'condition' 60 'group' 60 }, ...
	{ 'index' 61 'load' [datafolder filesep '106531t' filesep 'preprocessed.106531t.set'] 'subject' '106531t' 'condition' 61 'group' 61 }, ...
	{ 'index' 62 'load' [datafolder filesep '107540t' filesep 'preprocessed.107540t.set'] 'subject' '107540t' 'condition' 62 'group' 62 }, ...
	{ 'index' 63 'load' [datafolder filesep '109551t' filesep 'preprocessed.109551t.set'] 'subject' '109551t' 'condition' 63 'group' 63 }, ...
	{ 'index' 64 'load' [datafolder filesep '114570t' filesep 'preprocessed.114570t.set'] 'subject' '114570t' 'condition' 64 'group' 64 }, ...
	{ 'index' 65 'load' [datafolder filesep '116580t' filesep 'preprocessed.116580t.set'] 'subject' '116580t' 'condition' 65 'group' 65 }, ...
	{ 'index' 66 'load' [datafolder filesep '117591t' filesep 'preprocessed.117591t.set'] 'subject' '117591t' 'condition' 66 'group' 66 }, ...
	{ 'index' 67 'load' [datafolder filesep '127640t' filesep 'preprocessed.127640t.set'] 'subject' '127640t' 'condition' 67 'group' 67 }, ...
	{ 'index' 68 'load' [datafolder filesep '128641t' filesep 'preprocessed.128641t.set'] 'subject' '128641t' 'condition' 68 'group' 68 }, ...
	{ 'index' 69 'load' [datafolder filesep '130651t' filesep 'preprocessed.130651t.set'] 'subject' '130651t' 'condition' 69 'group' 69 }, ...
	{ 'index' 70 'load' [datafolder filesep '131660t' filesep 'preprocessed.131660t.set'] 'subject' '131660t' 'condition' 70 'group' 70 }, ...
	{ 'index' 71 'load' [datafolder filesep '132661t' filesep 'preprocessed.132661t.set'] 'subject' '132661t' 'condition' 71 'group' 71 }, ...
	{ 'index' 72 'load' [datafolder filesep '134670rt' filesep 'preprocessed.134670rt.set'] 'subject' '134670rt' 'condition' 72 'group' 72 }, ...
	{ 'index' 73 'load' [datafolder filesep '137691t' filesep 'preprocessed.137691t.set'] 'subject' '137691t' 'condition' 73 'group' 73 }, ...
	{ 'index' 74 'load' [datafolder filesep '138690t' filesep 'preprocessed.138690t.set'] 'subject' '138690t' 'condition' 74 'group' 74 }, ...
	{ 'index' 75 'load' [datafolder filesep '139701t' filesep 'preprocessed.139701t.set'] 'subject' '139701t' 'condition' 75 'group' 75 }, ...
	{ 'index' 76 'load' [datafolder filesep '150751t' filesep 'preprocessed.150751t.set'] 'subject' '150751t' 'condition' 76 'group' 76 }, ...
	{ 'index' 77 'load' [datafolder filesep '151760t' filesep 'preprocessed.151760t.set'] 'subject' '151760t' 'condition' 77 'group' 77 }, ...
	{ 'index' 78 'load' [datafolder filesep '153770t' filesep 'preprocessed.153770t.set'] 'subject' '153770t' 'condition' 78 'group' 78 }, ...
	{ 'index' 79 'load' [datafolder filesep '160800t' filesep 'preprocessed.160800t.set'] 'subject' '160800t' 'condition' 79 'group' 79 }, ...
	{ 'index' 80 'load' [datafolder filesep '67341t' filesep 'preprocessed.67341t.set'] 'subject' '67341t' 'condition' 80 'group' 80 }, ...
	{ 'index' 81 'load' [datafolder filesep '72360t' filesep 'preprocessed.72360t.set'] 'subject' '72360t' 'condition' 81 'group' 81 }, ...
	{ 'index' 82 'load' [datafolder filesep '83420t' filesep 'preprocessed.83420t.set'] 'subject' '83420t' 'condition' 82 'group' 82 }, ...
	{ 'index' 83 'load' [datafolder filesep '86430t' filesep 'preprocessed.86430t.set'] 'subject' '86430t' 'condition' 83 'group' 83 }, ...
	{ 'index' 84 'load' [datafolder filesep '88440t' filesep 'preprocessed.88440t.set'] 'subject' '88440t' 'condition' 84 'group' 84 }, ...
	{ 'index' 85 'load' [datafolder filesep '75380t' filesep 'preprocessed.75380t.set'] 'subject' '75380t' 'condition' 85 'group' 85 }, ...
	{ 'index' 86 'load' [datafolder filesep '76381t' filesep 'preprocessed.76381t.set'] 'subject' '76381t' 'condition' 86 'group' 86 }, ...
	{ 'index' 87 'load' [datafolder filesep '12061t' filesep 'preprocessed.12061t.set'] 'subject' '12061t' 'condition' 87 'group' 87 }, ...
	{ 'index' 88 'load' [datafolder filesep '15081t' filesep 'preprocessed.15081t.set'] 'subject' '15081t' 'condition' 88 'group' 88 }, ...
	{ 'index' 89 'load' [datafolder filesep '16080t' filesep 'preprocessed.16080t.set'] 'subject' '16080t' 'condition' 89 'group' 89 }, ...
	{ 'index' 90 'load' [datafolder filesep '30151t' filesep 'preprocessed.30151t.set'] 'subject' '30151t' 'condition' 90 'group' 90 }, ...
	{ 'index' 91 'load' [datafolder filesep '36181t' filesep 'preprocessed.36181t.set'] 'subject' '36181t' 'condition' 91 'group' 91 }, ...
	{ 'index' 92 'load' [datafolder filesep '37191t' filesep 'preprocessed.37191t.set'] 'subject' '37191t' 'condition' 92 'group' 92 }, ...
	{ 'index' 93 'load' [datafolder filesep '40200t' filesep 'preprocessed.40200t.set'] 'subject' '40200t' 'condition' 93 'group' 93 }, ...
	{ 'index' 94 'load' [datafolder filesep '47241t' filesep 'preprocessed.47241t.set'] 'subject' '47241t' 'condition' 94 'group' 94 }, ...
	{ 'index' 95 'load' [datafolder filesep '51260t' filesep 'preprocessed.51260t.set'] 'subject' '51260t' 'condition' 95 'group' 95 }, ...		
	{ 'dipselect' 0.15 } });

% Use this STUDY structure for fieldtrip nearest neighbour file
[STUDY neighbors] = std_prepare_neighbors( STUDY, ALLEEG , 'force', 'on');

% FEed this fieldtrip neighbour file into cluster-based permutation tests
