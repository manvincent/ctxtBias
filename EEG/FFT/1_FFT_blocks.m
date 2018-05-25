
%% Set-up parameters
clear all
clc

home = '/Volumes/MacintoshHD4/RFlex2';
addpath '/Users/eeg/Documents/MATLAB/fieldtrip-master';
addpath(genpath('/Users/eeg/Documents/MATLAB/eeglab13_4_4b'));
ft_defaults

% Defining directories
datafolder = [home filesep 'ParentRFlex'];
parentfolder = [home filesep 'FFT_analysis'];
cd(parentfolder)
% load EEG location file
load([parentfolder filesep 'EEGchanloc.mat']);

% Setting up subject loop
subject_list = textread([parentfolder filesep 'subjectlist_FFT.txt'], '%s');

% Sampling rate of preprocessed data
padN = 200; % number of actual data points per block, after zero-padding
% Define frequency bands
freq_list = [];
delta = [1,3];
theta = [4,7];
alpha = [8,12];
beta = [16,30];
lowGamma = [30,59];
medGamma = [61,90];
highGamma = [90,120];
sanityGamma = [32,100]
freq_list = vertcat(delta, theta, alpha, beta, lowGamma, medGamma, highGamma, sanityGamma);
freq_names = {'delta', 'theta', 'alpha', 'beta', 'lowGamma','medGamma'} %R,'highGamma','sanityGamma'}
normMethod = 'mean' % Options are 'mean' or 'dB' or 'median'

% Specify trigger codes for baseline normalisation and analyses
% Create acceptable EEG trigger arrays
OO_code = [];
FF_code = [];
OF_code = [];
FO_code = [];
OO_end = [];
FF_end = [];
OF_end = [];
FO_end = [];
ctxt_start_codes = [];
ctxt_end_codes = [];

OO_code = 91;
FF_code = 92;
OF_code = 93;
FO_code = 94;
OO_end = 88; %791;
FF_end = 88; %792;
OF_end = 88; %793;
FO_end = 88; %794;
ctxt_start_codes = [OO_code; FF_code; OF_code; FO_code];
ctxt_end_codes = [OO_end; FF_end; OF_end; FO_end];

recycle('off')
        
%% FFT analysis
parfor s=1:length(subject_list); s
    % Defining the 'subject' variable
    subject = [];
    subjectfolder = [];

    subject = subject_list{s}; 
    subjectfolder = [parentfolder filesep subject];
    if ~exist(subjectfolder,'dir');
        mkdir(subjectfolder)
    end

    fprintf('\n\n\n*** FFT analysis on subject %d (%s) ***\n\n\n', s, subject);
    % Load in preprocessed data from EEGLAB 
    EEG = [];
    EEG = pop_loadset('filename',['preprocessed.' subject '.set'],'filepath', [datafolder filesep subject],'loadmode','all');
    EEG = eeg_checkset( EEG ); 

    % Convert EEGLAB structure to FieldTrip sturcture
    PreprocessedData = [];
    PreprocessedData = eeglab2fieldtrip(EEG,'preprocessing','none')

    % Load in subject's onset data
    Eventlist = [];
    Eventlist = load([home filesep 'ParentRFlex/Eventlist_Structures' filesep subject '_eventinfo.mat']);
    Eventlist = Eventlist.Eventlist;

    
    % Initalize structure to hold data for subject
    subCell = cell([3,4]); % Structure will hold 3(row - block instance) * 4(column - context condition) matrices
    subCellNorm = cell([3,4]);    
    subCellFreqBand = cell([3,4]);
    subCellFreqBandNorm = cell([3,4]);
    for c=1:length(ctxt_start_codes) % Loop through the 4 contexts
        currContext = [];
        if ctxt_start_codes(c) == 91
            currContext = 'OO';
        elseif ctxt_start_codes(c) == 92
            currContext = 'FF';
        elseif ctxt_start_codes(c) == 93
            currContext = 'OF';
        elseif ctxt_start_codes(c) == 94
            currContext = 'FO';
        end 
        fprintf('\n*** Analyzing block ID %s ***\n', currContext);

        % Create block start vector                                 
        block_start_vec = [];
        for i = 1:numel(Eventlist)
            starttime = [];
            if sum(any(Eventlist(i).code == ctxt_start_codes(c))) > 0                    
                starttime = round(Eventlist(i).time,3);
                block_start_vec = [block_start_vec, starttime];
            end
        end
        
        % Create block end vector -- create 2 minute (120second) 'blocks' 
        block_end_vec = [];
        for i = 1:numel(block_start_vec)
            endtime = [];
            endtime = block_start_vec(i) + 150;
            block_end_vec = [block_end_vec, endtime];
        end                
     
    
        %Loop through each block within a context (3 block instances per context)
        for b = 1:numel(block_start_vec)
            val = [];
            % Mark the index in the EEG data structure where each block starts and ends
            startTimeDiff = [];
            currBlockStart = [];
            wind_start_idx = [];
            currBlockStart = block_start_vec(b);
            startTimeDiff = abs(EEG.times-currBlockStart);
            [val,wind_start_idx] = min(startTimeDiff);

            endTimeDiff = [];
            currBlockEnd = [];
            wind_end_idx = [];
            currBlockEnd = block_end_vec(b);
            endTimeDiff = abs(EEG.times-currBlockEnd);
            [val,wind_end_idx] = min(endTimeDiff);
          
            blockData = []; 
            blockData = EEG.data(:,wind_start_idx:wind_end_idx);          % Dims here are 64 rows * time columns                
            % Run FFT on block
            fftBlock = fft(blockData',padN)';      % Dims here are 64 rows * freq columns
            % Here we are zero-padding the data to give us padN (200) data points in each block, which increases frequency resolution w/o changing information

            % Compute power spectrum density via: psdsin = 2 .* abs(fftsin) .^ 2 ./ (number of samples .^2)
            % fftIn = fourier values
            psdBlock = 2 .* abs(fftBlock) .^ 2 ./ (size(blockData,2) .^2);       % Dims here are 64 rows * freq columns (same as FFT results)                
            % Trim psdBlock to match frequency resolution
            psdBlock = psdBlock(:,1:(padN/2)+1);
            subCell{b,c} = double(psdBlock);
            
            % Compute frequency vector (max frequency = (number of samples / 2) + 1
            freqaxis = linspace(0,(padN/2),(padN/2)+1);
                                          
             % 1/f normalize the power spectrum (multiply power by its corresponding frequency value)
            psdBlockNorm = bsxfun(@times,freqaxis,psdBlock);
            subCellNorm{b,c} = double(psdBlockNorm);
            
            psdBandNorm = [];
            psdBand = [];
            for band = 1:size(freq_list,1)
                val = [];
                % band number corresponds to: (1)delta, (2)theta, (3)alpha, (4)beta, (5)gamma
                minFreq = []; 
                maxFreq = []; 
                minFreqIdx = [];
                maxFreqIdx = [];                    
                minFreq = min(freq_list(band,:)) ;
                maxFreq = max(freq_list(band,:)) ;                     
                [val,minFreqIdx] = min(abs(freqaxis-minFreq));
                [val,maxFreqIdx] = min(abs(freqaxis-maxFreq));   
                % Get subset of data across all frequencies that fall into the band
                currBandNorm = [];
                currBandNorm = psdBlockNorm(:,minFreqIdx:maxFreqIdx);   
                % Get version that is not 1/f normalized for decibel conversion baseline normalisation below
                currBand = [];
                currBand = psdBlock(:,minFreqIdx:maxFreqIdx);
                % Aggregate power across frequency band
                % Compute mean power across current frequency band
                meanFreqBandNorm = [];
                meanFreqBandNorm = mean(currBandNorm,2);
                meanFreqBand = [];
                meanFreqBand = mean(currBand,2);
                % Concatenate to make a 64*numFreqBand matrix
                psdBandNorm = [psdBandNorm, meanFreqBandNorm];
                psdBand = [psdBand, meanFreqBand];
            end                                                    
            % Save separate FFT results
            subCellFreqBandNorm{b,c} = double(psdBandNorm);                            
            subCellFreqBand{b,c} = double(psdBand);            
        end
    end

    % Normalise using the mean across all block instances and contexts for 'baseline' (keep EEG channel and freq band independent)
    if strcmp(normMethod,'mean') == 1
        meanPowerNorm = mean(cat(3,subCellFreqBandNorm{:}),3);   
        sdPowerNorm = std(cat(3,subCellFreqBandNorm{:}),0,3);
        % Mean center the power for this subject
        for r = 1:size(subCellFreqBand,1)
            for c = 1:size(subCellFreqBand,2)
                subCellFreqBandNorm{r,c} = subCellFreqBandNorm{r,c} - meanPowerNorm; % mean-centre
                % Change outlier (+/- 3e s.d.) power to mean            
                for e = 1:size(sdPowerNorm,1)
                    for f = 1:size(sdPowerNorm,2)
                        if subCellFreqBandNorm{r,c}(e,f) > sdPowerNorm(e,f) | subCellFreqBandNorm{r,c}(e,f) < -1*sdPowerNorm(e,f)
                            subCellFreqBandNorm{r,c}(e,f) = meanPowerNorm(e,f);
                        end
                    end
                end
            end
        end
        % Save output to mat 
        outTable = cell2table(subCellNorm, 'VariableNames',{'OO','FF','OF','FO'});
        parsave([subjectfolder filesep subject '_FFT_1fNorm.mat'],1,outTable,'-v7.3');

        outTableFreqBand = cell2table( subCellFreqBandNorm, 'VariableNames',{'OO','FF','OF','FO'});
        parsave([subjectfolder filesep subject '_FFT_FreqBand_MeanCenter.mat'],1,outTableFreqBand,'-v7.3');    
    elseif strcmp(normMethod,'median') == 1 
        sdPowerNorm = std(cat(3,subCellFreqBandNorm{:}),0,3);
        medianPowerNorm = median(cat(3,subCellFreqBandNorm{:}),3);
        % Mean center the power for this subject
        for r = 1:size(subCellFreqBand,1)
            for c = 1:size(subCellFreqBand,2)
                subCellFreqBandNorm{r,c} = subCellFreqBandNorm{r,c} - medianPowerNorm; % mean-centre
                % Change outlier (+/- 3 s.d.) power to mean            
                for e = 1:size(sdPowerNorm,1)
                    for f = 1:size(sdPowerNorm,2)
                        if subCellFreqBandNorm{r,c}(e,f) > sdPowerNorm(e,f) | subCellFreqBandNorm{r,c}(e,f) < -1*sdPowerNorm(e,f)
                            subCellFreqBandNorm{r,c}(e,f) = medianPowerNorm(e,f);
                        end
                    end
                end
            end
        end
        % Save output to mat 
        outTable = cell2table(subCellNorm, 'VariableNames',{'OO','FF','OF','FO'});
        parsave([subjectfolder filesep subject '_FFT_1fNorm.mat'],1,outTable,'-v7.3');

        outTableFreqBand = cell2table( subCellFreqBandNorm, 'VariableNames',{'OO','FF','OF','FO'});
        parsave([subjectfolder filesep subject '_FFT_FreqBand_MedianCenter.mat'],1,outTableFreqBand,'-v7.3');
    elseif strcmp(normMethod,'dB') == 1 
        meanPower = mean(cat(3,subCellFreqBand{:}),3);
        % Decibel Conversion 
        for r = 1:size(subCellFreqBand,1)
            for c = 1:size(subCellFreqBand,2)
                subCellFreqBand{r,c} = 10*log10( subCellFreqBand{r,c} ./ meanPower); 
            end
        end
        % Save output to mat         
        outTableFreqBand = cell2table( subCellFreqBand, 'VariableNames',{'OO','FF','OF','FO'});
        parsave([subjectfolder filesep subject '_FFT_FreqBand_DecibNorm.mat'],1,outTableFreqBand,'-v7.3');     
    end

    fprintf('\n Saved FFT results for subject %d (%s)\n', s, subject);        
end


%% Create MLM-ready datasheets, concatenate across subjects, separate columns for each electrode, one per frequency band 

parfor band = 1:size(freq_list,1)
	currBandName = freq_names(band)

	sub_concatData = [];
	for s=1:length(subject_list); 

	    % Defining the 'subject' variable
	    subject = [];
	    subjectfolder = [];

	    subject = subject_list{s}; 
	    subjectfolder = [parentfolder filesep subject];
	    if ~exist(subjectfolder,'dir');
	        mkdir(subjectfolder);
	    end

		% Input subject data
        if strcmp(normMethod,'mean') == 1
            subData = load([subjectfolder filesep subject '_FFT_FreqBand_MeanCenter.mat']);
        elseif strcmp(normMethod,'median') == 1 
            subData = load([subjectfolder filesep subject '_FFT_FreqBand_MedianCenter.mat']);
        elseif strcmp(normMethod,'dB') == 1
            subData = load([subjectfolder filesep subject '_FFT_FreqBand_DecibNorm.mat']);
        end

        ctxtConcat = [];
		for c=1:length(ctxt_start_codes) % Loop through the 4 contexts
            currContext = [];
            if ctxt_start_codes(c) == 91
                currContext = 'OO';
            elseif ctxt_start_codes(c) == 92
                currContext = 'FF';
            elseif ctxt_start_codes(c) == 93
                currContext = 'OF';
            elseif ctxt_start_codes(c) == 94
                currContext = 'FO';
            end 

            blockConcat = [];
            for b = 1:3
                band_Data = cell([1,66]);
                band_Data{1,1} = subject;
                band_Data{1,2} = currContext;
                for e = 1:64
                    band_Data{1,2+e} = subData.outTableFreqBand{b,c}{1}(e,band);
                end
                % Concatenate block instances
                blockConcat = vertcat(blockConcat,band_Data);
            end
        % Concatenate context types
        ctxtConcat = vertcat(ctxtConcat, blockConcat);        
        end
    % Concatenate subject
    sub_concatData = vertcat(sub_concatData,ctxtConcat);
    end
    
    
    % Convert to long format table 
    % Read in electrode location labels    
    tableHeader = cell([1,66]);
    tableHeader{1,1} = 'subj_idx';
    tableHeader{1,2} = 'context';
    for e = 1:64
        tableHeader{1,2+e}=loc(e).labels;
    end
    
    % Convert to table
    outTableSub_band = cell2table( sub_concatData,'VariableNames',tableHeader);
  
  
    % Add contrasts for MLM
    promC = rand([size(outTableSub_band,1),1]);
    prevC = rand([size(outTableSub_band,1),1]);   
    interC = rand([size(outTableSub_band,1),1]);   
    % C1: High gain val > low gain val
    promC( strcmp(outTableSub_band.context,'OO')==1) = 1;
    promC( strcmp(outTableSub_band.context,'FF')==1) = -1;    
    promC( strcmp(outTableSub_band.context,'OF')==1) = 1;
    promC( strcmp(outTableSub_band.context,'FO')==1) = -1;
    % C2: High loss val > low loss val
    prevC( strcmp(outTableSub_band.context,'OO')==1) = 1;
    prevC( strcmp(outTableSub_band.context,'FF')==1) = -1;
    prevC( strcmp(outTableSub_band.context,'OF')==1) = -1;    
    prevC( strcmp(outTableSub_band.context,'FO')==1) = 1;    
    % C3: Interaction: Asymmetric vs. Symmetric
    interC( strcmp(outTableSub_band.context,'OO')==1) = 1;
    interC( strcmp(outTableSub_band.context,'FF')==1) = 1;
    interC( strcmp(outTableSub_band.context,'OF')==1) = -1;      
    interC( strcmp(outTableSub_band.context,'FO')==1) = -1;
    % Make contrast table
    contrastTable = table(promC, prevC, interC,'VariableNames',{'Con1','Con2','Con3'});
    % Concatenate tables
    statTable = horzcat(outTableSub_band,contrastTable);

    % Save mat
    if strcmp(normMethod,'mean') == 1
        parsave([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_concatenatedData_mc.mat'],1,statTable,'-v7.3');
    elseif strcmp(normMethod,'median') == 1 
        parsave([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_concatenatedData_medc.mat'],1,statTable,'-v7.3');        
    elseif strcmp(normMethod,'dB') == 1
        parsave([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_concatenatedData_db.mat'],1,statTable,'-v7.3');
    end
               
%     % Run MLM
%     resultsP = zeros(64,3);
%     resultsFStat = zeros(64,3);    
%     for elec = 1:64
%         currElec = [];
%         elecDV = [];
%         formula = []; 
%         model = [];
%         formula = [];
%         currElec = loc(elec).labels;
%         elecDV = statTable.(currElec);
%         formula = [currElec '~Con1+Con2+Con3+(1|subj_idx)'];
%         model = anova(fitlme(statTable,formula),'dfmethod','satterthwaite');
%         
%         % Loop across contrast code
%         for code = 1:3
%             conP = [];
%             conTStat = [];
%             % pull out stats for contrast cond
%             conP = model.pValue(code+1);
%             conFStat = model.FStat(code+1);      
%             % store results in array that is 64 (electrodes) by 3 (contrasts)
%             resultsP(elec,code) = conP;
%             resultsFStat(elec,code) = conFStat;
%         end
%     end
%     % Save in group folder
%     if strcmp(normMethod,'mean') == 1
%         parsave([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_pValues_mc.mat'],1,resultsP,'-v7.3');
%         parsave([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_fStat_mc.mat'],1,resultsFStat,'-v7.3');        
%     elseif strcmp(normMethod,'median') == 1         
%         parsave([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_pValues_medc.mat'],1,resultsP,'-v7.3');
%         parsave([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_fStat_medc.mat'],1,resultsFStat,'-v7.3'); 
%     elseif strcmp(normMethod,'dB') == 1
%         parsave([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_pValues_db.mat'],1,resultsP,'-v7.3');
%         parsave([parentfolder filesep 'Group' filesep 'Freq_' char(currBandName) '_fStat_db.mat'],1,resultsFStat,'-v7.3');
%     end
    
        
    
    
end

    fprintf('\n Finished \n');


%% Turn recycling back on
recycle('on'); 
      
             
      



 
