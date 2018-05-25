%% Retrigger function
function EEG = retrigger(EEG, subject)

    if exist('subject') == 1 && isempty(subject) == 0
        fprintf('Subject file properly loaded\n')
    else fprintf('ERROR: No Subject specified\n')
    end
    sub_ID = [];
    sub_ID = str2num(subject(1:end-4));
    

          
    % Extracting the trigger codes out as a numerical array for manipulation
    code = [];
    for i = 1:length(EEG.EVENTLIST.eventinfo) 
        code_element = [];
        code_element = EEG.EVENTLIST.eventinfo(i).code;
        code = vertcat(code, code_element);
    end

% Re-trigger fixations from '1' to '9' if necessary in early subject
    if sub_ID >= 7 && sub_ID <= 11
        code(code==1) = 9;
    end
    
    
% Re-trigger responses from (2,3,4) to (79,80,81 - numpad) if necessary
    if any(code == 3)
        code(code == 2) = 79;
        code(code == 3) = 80;
        code(code == 4) = 81;    
    end

% Inserting the code array back into the EEG strucutre trigger code field
    for k = 1:length(EEG.EVENTLIST.eventinfo) 
        EEG.EVENTLIST.eventinfo(k).code = code(k,1);
    end    
    
    
 % Create a new vector that specifies the condition
    % Index when each new condition starts
    index_99 = [];
    index_91 = [];
    index_92 = [];
    index_93 = []; 
    index_94 = [];
    all_indices = [];
    
    index_99 = find(code == 99) ;    
    index_91 = find(code == 91) ;
    index_92 = find(code == 92) ;
    index_93 = find(code == 93) ;
    index_94 = find(code == 94) ;
    all_indices = vertcat(index_99, index_91, index_92, index_93, index_94);
   
    % Create the condition vector
    condition_vector = zeros(length(code),1);    
    for condition_onset = 1:numel(all_indices)
        curr_onset_index = [];
        curr_onset_index = all_indices(condition_onset);
        
        next_onset_index = [];
        next_onset_index = min(all_indices(find(all_indices>curr_onset_index)));
        
        if ~isempty(next_onset_index)
            condition_vector(curr_onset_index:next_onset_index,1) = code(curr_onset_index);
        elseif isempty(next_onset_index)
            condition_vector(curr_onset_index:end,1) = code(curr_onset_index);
        end       
    end
    

% Re-trigger triggers based on condition and counterbalance
    if mod(sub_ID,2) == 0 % If sub_ID is even
        sub_condition = 'LCG';
        fprintf('This subject has a counterbalance order of %s\n', sub_condition);
        
        % Add code labels to stimuli onset
        for j = 1:length(EEG.EVENTLIST.eventinfo) 
            if code(j,:) == 10
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_prac';
            elseif code(j,:) == 20
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_prac';
            elseif code(j,:) == 30
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_prac';
            elseif code(j,:) == 11
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_+1/-1';
            elseif code(j,:) == 12
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_+1/-1';
            elseif code(j,:) == 13
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_+1/-1';
            elseif code(j,:) == 21
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_+.5/-.5';
            elseif code(j,:) == 22
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_+.5/-.5';
            elseif code(j,:) == 23
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_+.5/-.5';
            elseif code(j,:) == 31
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_+1/-.5';
            elseif code(j,:) == 32
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_+1/-.5';
            elseif code(j,:) == 33
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_+1/-.5';
            elseif code(j,:) == 41
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_+.5/-1';
            elseif code(j,:) == 42
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_+.5/-1';
            elseif code(j,:) == 43
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_+.5/-1';
            elseif code(j,:) == 19
                EEG.EVENTLIST.eventinfo(j).codelabel = 'FB_LossCorr';
            elseif code(j,:) == 29
                EEG.EVENTLIST.eventinfo(j).codelabel = 'FB_NeutCorr';
            elseif code(j,:) == 39
                EEG.EVENTLIST.eventinfo(j).codelabel = 'FB_GainCorr';
            end                
        end
        
        % Add code labels to response
        for j = 1:length(EEG.EVENTLIST.eventinfo) 
            if condition_vector(j,:) == 91
                if code(j,:) == 79
                    if code(j-1,:) == 11
                        EEG.EVENTLIST.eventinfo(j).code = 56;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_resp_+1/-1';
                    else 
                        EEG.EVENTLIST.eventinfo(j).code = 7756;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Loss_resp_+1/-1';
                    end                    
                elseif code(j,:) == 80
                    if code(j-1,:) == 12
                        EEG.EVENTLIST.eventinfo(j).code = 55;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_resp_+1/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7755;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Neut_resp_+1/-1';
                    end                        
                elseif code(j,:) == 81
                     if code(j-1,:) == 13
                        EEG.EVENTLIST.eventinfo(j).code = 54;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_resp_+1/-1';
                     else
                        EEG.EVENTLIST.eventinfo(j).code = 7754;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Gain_resp_+1/-1';
                     end                     
                end
            elseif condition_vector(j,:) == 92
                if code(j,:) == 79
                    if code(j-1,:) == 21
                        EEG.EVENTLIST.eventinfo(j).code = 66;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_resp_+.5/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7766;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Loss_resp_+.5/-.5';
                    end                    
                elseif code(j,:) == 80
                    if code(j-1,:) == 22
                        EEG.EVENTLIST.eventinfo(j).code = 65;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_resp_+.5/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7765;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Neut_resp_+.5/-.5';
                    end
                elseif code(j,:) == 81
                    if code(j-1,:) == 23
                        EEG.EVENTLIST.eventinfo(j).code = 64;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_resp_+.5/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7764;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Gain_resp_+.5/-.5';
                    end                   
                end
            elseif condition_vector(j,:) == 93
                if code(j,:) == 79
                    if code(j-1,:) == 31
                        EEG.EVENTLIST.eventinfo(j).code = 76;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_resp_+1/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7776;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Loss_resp_+1/-.5';
                    end                    
                elseif code(j,:) == 80
                    if code(j-1,:) == 32
                        EEG.EVENTLIST.eventinfo(j).code = 75;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_resp_+1/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7775;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Neut_resp_+1/-.5';
                    end                    
                elseif code(j,:) == 81
                    if code(j-1,:) == 33
                        EEG.EVENTLIST.eventinfo(j).code = 74;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_resp_+1/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7774;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Gain_resp_+1/-.5';
                    end                                            
                end
            elseif condition_vector(j,:) == 94
                if code(j,:) == 79
                    if code(j-1,:) == 41
                        EEG.EVENTLIST.eventinfo(j).code = 86;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_resp_+.5/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7786;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Loss_resp_+.5/-1';                        
                    end                    
                elseif code(j,:) == 80
                    if code(j-1,:) == 42
                        EEG.EVENTLIST.eventinfo(j).code = 85;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_resp_+.5/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7785;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Neut_resp_+.5/-1';
                    end
                elseif code(j,:) == 81
                    if code(j-1,:) == 43
                        EEG.EVENTLIST.eventinfo(j).code = 84;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_resp_+.5/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7784;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Gain_resp_+.5/-1';                        
                    end
                end
            end
        end                                 
    elseif mod(sub_ID,2) == 1 % If sub_ID is odd
        sub_condition = 'GCL';
         fprintf('This subject has a counterbalance order of %s\n', sub_condition);
         
         % Add code labels to stimuli onset
        for j = 1:length(EEG.EVENTLIST.eventinfo) 
            if code(j,:) == 10
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_prac';
            elseif code(j,:) == 20
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_prac';
            elseif code(j,:) == 30
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_prac';
            elseif code(j,:) == 11
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_+1/-1';
            elseif code(j,:) == 12
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_+1/-1';
            elseif code(j,:) == 13
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_+1/-1';
            elseif code(j,:) == 21
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_+.5/-.5';
            elseif code(j,:) == 22
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_+.5/-.5';
            elseif code(j,:) == 23
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_+.5/-.5';
            elseif code(j,:) == 31
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_+1/-.5';
            elseif code(j,:) == 32
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_+1/-.5';
            elseif code(j,:) == 33
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_+1/-.5';
            elseif code(j,:) == 41
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_stim_+.5/-1';
            elseif code(j,:) == 42
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_stim_+.5/-1';
            elseif code(j,:) == 43
                EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_stim_+.5/-1';
            elseif code(j,:) == 19
                EEG.EVENTLIST.eventinfo(j).codelabel = 'FB_GainCorr';
            elseif code(j,:) == 29
                EEG.EVENTLIST.eventinfo(j).codelabel = 'FB_NeutCorr';
            elseif code(j,:) == 39
                EEG.EVENTLIST.eventinfo(j).codelabel = 'FB_LossCorr';
            end                
        end

                % Add code labels to response
        for j = 1:length(EEG.EVENTLIST.eventinfo) 
            if condition_vector(j,:) == 91
                if code(j,:) == 79
                    if code(j-1,:) == 11
                        EEG.EVENTLIST.eventinfo(j).code = 54;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_resp_+1/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7754;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Gain_resp_+1/-1';
                    end
                elseif code(j,:) == 80
                    if code(j-1,:) == 12
                        EEG.EVENTLIST.eventinfo(j).code = 55;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_resp_+1/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7755;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Neut_resp_+1/-1';
                    end
                elseif code(j,:) == 81
                    if code(j-1,:) == 13
                        EEG.EVENTLIST.eventinfo(j).code = 56;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_resp_+1/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7756;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Loss_resp_+1/-1';
                    end                    
                end
            elseif condition_vector(j,:) == 92
                if code(j,:) == 79
                    if code(j-1,:) == 21
                        EEG.EVENTLIST.eventinfo(j).code = 64;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_resp_+.5/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7764;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Gain_resp_+.5/-.5';
                    end
                elseif code(j,:) == 80
                    if code(j-1,:) == 22
                        EEG.EVENTLIST.eventinfo(j).code = 65;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_resp_+.5/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7765;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Neut_resp_+.5/-.5';
                    end
                elseif code(j,:) == 81
                    if code(j-1,:) == 23
                        EEG.EVENTLIST.eventinfo(j).code = 66;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_resp_+.5/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7766;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Loss_resp_+.5/-.5';
                    end                                            
                end
            elseif condition_vector(j,:) == 93
                if code(j,:) == 79
                    if code(j-1,:) == 31
                        EEG.EVENTLIST.eventinfo(j).code = 74;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_resp_+1/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7774;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Gain_resp_+1/-.5';
                    end
                elseif code(j,:) == 80
                    if code(j-1,:) == 32
                        EEG.EVENTLIST.eventinfo(j).code = 75;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_resp_+1/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7775;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Neut_resp_+1/-.5';
                    end
                elseif code(j,:) == 81
                    if code(j-1,:) == 33
                        EEG.EVENTLIST.eventinfo(j).code = 76;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_resp_+1/-.5';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7776;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Loss_resp_+1/-.5';
                    end
                end
            elseif condition_vector(j,:) == 94
                if code(j,:) == 79
                    if code(j-1,:) == 41
                        EEG.EVENTLIST.eventinfo(j).code = 84;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Gain_resp_+.5/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7784;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Gain_resp_+.5/-1';
                    end                    
                elseif code(j,:) == 80
                    if code(j-1,:) == 42
                        EEG.EVENTLIST.eventinfo(j).code = 85;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Neut_resp_+.5/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7785;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Neut_resp_+.5/-1';
                    end                    
                elseif code(j,:) == 81
                    if code(j-1,:) == 43
                        EEG.EVENTLIST.eventinfo(j).code = 86;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'Loss_resp_+.5/-1';
                    else
                        EEG.EVENTLIST.eventinfo(j).code = 7786;
                        EEG.EVENTLIST.eventinfo(j).codelabel = 'ERR_Loss_resp_+.5/-1';
                    end
                end
            end
        end
         
    end
    
    % Add labels to condition starts, fixations, ends, and clean up boundary events
    for j = 1:length(EEG.EVENTLIST.eventinfo) 
        if code(j,:) == 9
            EEG.EVENTLIST.eventinfo(j).codelabel = 'Fixation';
        elseif code(j,:) == 99
            EEG.EVENTLIST.eventinfo(j).codelabel = 'Prac_Block';
        elseif code(j,:) == 91
            EEG.EVENTLIST.eventinfo(j).codelabel = '+1/-1_Block';
        elseif code(j,:) == 92
            EEG.EVENTLIST.eventinfo(j).codelabel = '+.5/-.5_Block';
        elseif code(j,:) == 93
            EEG.EVENTLIST.eventinfo(j).codelabel = '+1/-.5_Block';
        elseif code(j,:) == 94
            EEG.EVENTLIST.eventinfo(j).codelabel = '+.5/-1_Block';   
        elseif code(j,:) == 7
            if condition_vector(j,:) == 99
                EEG.EVENTLIST.eventinfo(j).codelabel = 'End_Instructions';
                EEG.EVENTLIST.eventinfo(j).code = 799;
            elseif condition_vector(j,:) == 91
                EEG.EVENTLIST.eventinfo(j).codelabel = 'End_Block_+1/-1';
                EEG.EVENTLIST.eventinfo(j).code = 791;
            elseif condition_vector(j,:) == 92
                EEG.EVENTLIST.eventinfo(j).codelabel = 'End_Block_+.5/-.5';
                EEG.EVENTLIST.eventinfo(j).code = 792;
            elseif condition_vector(j,:) == 93
                EEG.EVENTLIST.eventinfo(j).codelabel = 'End_Block_+1/-.5';
                EEG.EVENTLIST.eventinfo(j).code = 793;
            elseif condition_vector(j,:) == 94
                EEG.EVENTLIST.eventinfo(j).codelabel = 'End_Block_+.5/-1';
                EEG.EVENTLIST.eventinfo(j).code = 794;
            end                                  
        elseif code(j,:) == 88
            EEG.EVENTLIST.eventinfo(j).codelabel = 'End_Run';
        elseif code(j,:) == -99
            EEG.EVENTLIST.eventinfo(j).code = 0;           
        end
    end
    
    % Writing the new triggers into the EEG structure
    EEG = pop_overwritevent(EEG,'code');
    %EEG = pop_overwritevent(EEG,'codelabel');
end

