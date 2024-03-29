function detectStruct = detection_plm_calculator(channel_cell_data, optional_params, stageStruct)
% detects PLM using selected detection and preprocessing methods
%
%
%
% HR is calculalated using
% channcel_indices(2) (ECG) (see rr_simple).  Timelocked cardiac
% accelerations (cal for cardiovascular acceleration/arousal lock)
% parameters are computed from matched PLM-cardiac cycles using the methods
% introduced by Winkelman in his 1999 paper.
%
%  channel_indices(1) = LAT/RAT channel
%  channel_indices(2) = ECG 
%
% Written by Hyatt Moore IV, 6/15/2012
%
% updated later on

% 1/9/2013 -  adjustable noise floor with rules
%
% 6/12/2013 - ported to single file for standalone use


%this allows direct input of parameters from outside function calls, which
%can be particularly useful in the batch job mode
if(nargin>=2 && ~isempty(optional_params))
    params = optional_params;
else
    pfile = strcat(mfilename('fullpath'),'.plist');
    
    if(exist(pfile,'file'))
        %load it
        params = plist.loadXMLPlist(pfile);
    else

        params.filter_with_anc = 1; %adaptively filter ECG (0 = no; 1 = yes)
        params.min_duration_sec = 0.75;
        params.average_power_window_sec = 20;  %calculate average power over consecutive windows of this duration in seconds
        params.merge_within_sec = 2.0;
        params.use_summer = 1;
        params.median_removal = 1;        
        plist.saveXMLPlist(pfile,params);
    end
end

samplerate = params.samplerate;

if(~iscell(channel_cell_data))
    channel_cell_data = {channel_cell_data};
end

%Heart rate detections
params.HR_seconds_after = 10; %number of seconds to look at HR variability after LM hit
params.HR_seconds_before = 10; %number of seconds to look at before LM hit

% upper and lower thresholds are based on AASM 2007 scoring criteria here
params.threshold_high_uV = 8; %8 uV above resting baseline
params.threshold_low_uV = 2;
params.max_duration_sec = 10.0;
params.summer_order = 2;
params.noisefloor_uV_to_turnoff_detection = 50;


%adaptive noise cancel if selected
if(params.filter_with_anc)
    %1. adaptively cancel noise from ECG using recursive least squares method
    data = filter.anc_rls(channel_cell_data{1},channel_cell_data{2});    
else
    data = channel_cell_data{1};
end

lm_detectStruct = detection.detection_lm_dualthresh_twopass_variable_noisefloor_mediator(data,params);
        

if(isempty(lm_detectStruct.new_events))
    %just return the emptiness
    detectStruct = lm_detectStruct;
else
    %% remove LMs that occur during Stage 7 and before sleep onset
    params.stages2exclude = 7;
    
    %I think this is faster ....
    firstNonWake = 1;
%     while((stageStruct.line(firstNonWake)==7||stageStruct.line(firstNonWake)==0) && firstNonWake<=numel(stageStruct.line))
%         firstNonWake = firstNonWake+1;
%     end
    
    %This is simpler and vectorized, but it actually has three operations that
    %run through the entire MARKING.sev_STAGES.line vector which is not necessary.
    %firstNonWake = find(MARKING.sev_STAGES.line~=0 & MARKING.sev_STAGES.line~=7,1);
    
    %convert these to the epochs the events occur in

    epochs = sample2epoch(lm_detectStruct.new_events,stageStruct.standard_epoch_sec,params.samplerate);
    
    %obtain the stage scores for these epochs
    stages = stageStruct.line(epochs);
    
    if(size(stages,2)==1)
        stages = stages'; %make it a row matrix
    end
    exclude_indices = false(size(epochs,1),2); %have to take into account start and stops being in different epochs
    
    
    for k=1:numel(params.stages2exclude)
        stage2exclude = params.stages2exclude(k);
        exclude_indices = exclude_indices|stages==stage2exclude;
    end
    
    %remove the firstNonWake indices as well here
    exclude_indices = exclude_indices(:,1)|exclude_indices(:,2)|epochs(:,1)<firstNonWake|epochs(:,2)<firstNonWake;
    
    %remove if LM's are not good... based on new criteria...
    
    
    paramStruct = lm_detectStruct.paramStruct;
    if(~isempty(paramStruct))
        fields = fieldnames(paramStruct);
        for k=1:numel(fields)
            paramStruct.(fields{k}) = paramStruct.(fields{k})(~exclude_indices);
        end
    end
    
    lm_events = lm_detectStruct.new_events(~exclude_indices,:);
    
    
    %% Cardiac Analysis
    %         params.filter_order = 10;  %leave this alone - necessary for rr
    %         detector, but will not allow it to be adjusted here
    
    rr_params.filter_order = 10;
    rr_params.samplerate = samplerate;
    hr_detectStruct = detection.detection_rr_simple(channel_cell_data{2},rr_params);
    hr_start_vec = hr_detectStruct.new_events(:,1);
    inst_hr = hr_detectStruct.paramStruct.inst_hr;
    
    num_lm = size(lm_events,1);
    HR_min = zeros(num_lm,1);
    HR_max = zeros(num_lm,1);
    HR_lead = zeros(num_lm,1);
    HR_lag = zeros(num_lm,1);
    
    num_hr_cycles = numel(inst_hr);
    
    %this loop took 1.31 seconds for an RLS case - don't bother optimizing
    for p=1:num_lm
        j = getPrecedingCardiacCycle(lm_events(p,1),hr_start_vec);
        base_hr = inst_hr(j); %baseline heart rate defined as hr immediately preceeding PLM onset/activation
        start_bound = max(j-9,1);
        stop_bound = min(j+10,num_hr_cycles);
        %     normalizedHR = inst_hr(start_bound:stop_bound)-base_hr;
        try
            [HR_min(p), HR_min_ind] = min(inst_hr(start_bound:j)-base_hr); %get the preceding minimum value
            [HR_max(p), HR_max_ind] = max(inst_hr(j+1:stop_bound)-base_hr); %get the ensuing maximum value
            HR_lead(p) = 11-HR_min_ind; % a result of 1 means that the base_hr was lowest
            HR_lag(p) = HR_max_ind;
        catch ME
            disp(ME);
        end
    end
    
    HR_lead = -HR_lead; %negative cycles precede PLM onset
    HR_delta = HR_max-HR_min;
    HR_run = HR_lag-HR_lead;
    HR_slope = HR_delta./HR_run;
    
    %tack on the additional timelocked HR parameters that we found
    paramStruct.HR_min = HR_min;
    paramStruct.HR_max = HR_max;
    paramStruct.HR_lead = HR_lead;
    paramStruct.HR_lag = HR_lag;
    paramStruct.HR_delta = HR_delta;
    paramStruct.HR_run = HR_run;
    paramStruct.HR_slope = HR_slope;    
   
    detectStruct.new_data = lm_detectStruct.new_data;
    detectStruct.new_events = lm_events;
    
    detectStruct.paramStruct = paramStruct;
    
    %apply PLM rules now - handle this post hoc since various rules can be
    %applied to filter out unwanted LM from PLM inclusion depending on the
    %analysis
    
%     detectStruct = calculatePLMstruct(detectStruct,samplerate);
    
end

function p = getPrecedingCardiacCycle(plm_sample, hr_events)
%plm_sample is the index of a digital time signal
%hr_events is a vector of digital time starting events
%p is the index of hr_events whose value immediately precedes plm_sample
p = find(hr_events<plm_sample,1,'last');
