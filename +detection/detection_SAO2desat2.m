function detectStruct = detection_SAO2desat(channel_cell_data,optional_params,stageStruct)
% INPUT
% signal:   raw signal
% Fs:       sampling frequency, Hz
% desat_settings has the following fields:
%   desat.wbase = length of window for calculating base
%   desat.thresholdO2 = SAO2 threshold in percentage (AASM 2007 defines 4%
%           for hypopnea)
%   desat.noisestd = threshold when an epoch is defined as noise (threshold
%           when a values is below mean-(x*std)
%   desat.noisediff = threshold when an epoch is defined as noise (gradient
%           on desat signal)
%   desat.noisediffclose = threshold when an epoch is defined as noise
%           (points close to large gradient/noise on desat signal)



% Implementation by Henriette Koch, 10/28/2013.
if(nargin>=2 && ~isempty(optional_params))
    params = optional_params;
else
    
    pfile = strcat(mfilename('fullpath'),'.plist');
    
    if(exist(pfile,'file'))
        % Load parameters
        params = plist.loadXMLPlist(pfile);
    else
        % Make parameters and save it for the future
        % Window for baseline (minutes)
        params.wbase = 2;
        % Noise thresholds
        params.noisestd = 5; % threshold when an epoch is defined as noise (threshold
%           when a values is below mean-(x*std)
        params.noisediff = 5; % threshold when an epoch is defined as noise (gradient
%           on desat signal)
        params.noisediffclose = 5; % threshold when an epoch is defined as noise
%           (points close to large gradient/noise on desat signal)
        
        plist.saveXMLPlist(pfile,params);
    end
end

if(~iscell(channel_cell_data))
    channel_cell_data = {channel_cell_data};
end
data = channel_cell_data{1};
Fs = params.samplerate;


%% Lowpass filter
% Filter
params.lowpass = 2;         % lowpass cut frequency
params.filter_order = 10;   % filter order
delay_lp = (params.filter_order)/2;

% Design filter, uses Hamming window of length N+1.
lowpassfilter = fir1(params.filter_order,params.lowpass/Fs*2);

% Apply filter and account for the delay.
filtd = filter(lowpassfilter,1,data);
filtdat = [filtd(delay_lp+1:end); zeros(delay_lp,1)]; % data for plotting
filtdata = filtdat;

%% Thesholds     

% Take out epochs above 100% O2 sat.
inmean = mean(filtdata); % mean used to define noise threshold
o100 = find(filtdata > 100);

% Take out epochs below "x std" or large negative gradient.
% Standard deviation threshold
instd = params.noisestd*std(filtdata);                  % OPTIMER!!!!!!!!!!!!!!!!!
if instd < 10;
    instd = 10;
else
end
ustd = find(filtdata < inmean-instd);

% Negative gradient
indiff = abs(diff(filtdata));
thres_diff = params.noisediff*std(indiff);         % OPTIMER!!!!!!!!!!!!!!!!!
udiff = find(indiff > thres_diff);

filtdata(o100) = NaN;
filtdata(ustd) = NaN;
filtdata(udiff) = NaN;

% Moving window to define threshold (use "wbase" minute window)
w = 100; % to run loop faster and due to low Fs when recording SaO2.
base = zeros(params.wbase*Fs,1); ba = zeros(params.wbase*Fs*60,1); m = 1;
for i = 1:w:length(filtdata)-params.wbase*Fs*60; % step "w" e.g. 100 = 1 sek.
    bain = filtdata(i:i+params.wbase*Fs*60);
    ba = bain(~isnan(bain),1);
    base(m) = mean(ba);
    m = 1+m;
end

%% Desaturation

% Hypopnea definition
thresholdO2 = 4;

% Calculate desaturation threshold.
for i = 1:length(base)
    thres(i) = base(i)-thresholdO2; % pure 4% desaturation
%     thres(i) = base(i)-thresholdO2*base(i); % uncomment if 4% of threshold
end
base_s = reshape(repmat(base',w,1),1,[]);
thres_s = reshape(repmat(thres,w,1),1,[]);

% Extract periods of desaturation.
for i = 1:length(thres_s)
    if filtdata(params.wbase*Fs*60+i) <= thres_s(i);
        desat_shift(i) = 1;
    elseif isnan(filtdata(params.wbase*Fs*60+i))
        desat_shift(i) = desat_shift(i-1);
    else desat_shift(i) = 0;
    end
end
% Zeropad first "wbase" of signal.
desat_signal = [zeros(params.wbase*Fs*60,1) ; desat_shift'];

% Make vector containing desaturated periods
[idx_start val_start] = find(diff(desat_signal)==1);
[idx_end val_end] = find(diff(desat_signal)==-1);
ml = min([length(val_start) length(val_end)]); % ensures equal size
desaturation_sample = [idx_start(1:ml) idx_end(1:ml)]; % extract equal numbers of shift

% Extract noise periods (if gradient p/m Fs desat is large)
thres_close = params.noisediffclose*std(indiff);       % OPTIMER!!!!!!!!!!!!!!!!!
if isempty(desaturation_sample);
    desaturation_sample = [];
else
for i = 1:length(desaturation_sample(:,1));
desat(i,:) = any(indiff(desaturation_sample(i,1)-2*Fs:desaturation_sample(i,2)+2*Fs) < ...
    thres_close)*desaturation_sample(i,:);
end
end
desat(~any(desat,2),:) = [];


%% OUTPUT
% Events
detectStruct.new_events = desat;
% Signal
detectStruct.new_data = filtdat;
% Event features
duration = (desat(:,2)-desat(:,1))/Fs; % event duration in sec.
detectStruct.paramStruct.duration = duration;



% 
% %% PLOT
% eh = 45; % epoch number
% dh = 12; % hour divided
% plot(ini_in(wbase*Fs*60+1:length(in)),'b')
% hold on
% for i = 1:length(desat(:,1))
%     plot(desat(i,1)-wbase*Fs*60:desat(i,2)-wbase*Fs*60,thres_s(desat(i,1)-wbase*Fs*60:desat(i,2)-wbase*Fs*60),'.r','markersize',10)%,'linewidth',12)
% end
% hold off
% axis([(eh-1)*(Fs*3600/dh)+1-wbase*Fs*60 (eh-1)*(Fs*3600/dh)+(Fs*3600/dh)-wbase*Fs*60 90 101]) % one minut
% 
% % PLOT
% % epoch hour
% 
% eh = 74; % epoch number
% dh = 12; % hour divided
% % Input signal, baseline and threshold
% plot(ini_in(wbase*Fs*60+1:length(in)),'b')
% axis([(eh-1)*(Fs*3600/dh)+1-wbase*Fs*60 (eh-1)*(Fs*3600/dh)+(Fs*3600/dh)-wbase*Fs*60 90 101]) % one minut
% hold on
% plot(thres_s,'g')
% plot(base_s,'r')
% % plot([1 length(base_s)], [inmean-instd inmean-instd],'--b') % noise std
% 
% % Diff
% % plot(10*indiff(wbase*Fs*60-1:length(in)-1)+99)
% % plot([1 length(base_s)],[10*thres_diff+99 10*thres_diff+99],'k')
% % plot([1 length(base_s)],[10*thres_close+99 10*thres_close+99],'--k')
% 
% % Desat samples
% for i = 1:length(desat(:,1))
%     plot(desat(i,1)-wbase*Fs*60:desat(i,2)-wbase*Fs*60,thres_s(desat(i,1)-wbase*Fs*60:desat(i,2)-wbase*Fs*60),'.r','markersize',10)%,'linewidth',12)
% end
% % plot(desat(:,2)-wbase*Fs*60,thres_s(desat(:,2)-wbase*Fs*60),'*g')
% % plot(desat(:,1)-wbase*Fs*60,thres_s(desat(:,1)-wbase*Fs*60),'*r')
% 
% hold off
% grid on
% legend('SaO2','Threshold','Baseline','Desaturation','Orientation','horizontal') %'5 std (noise threshold)',
% axis([(eh-1)*(Fs*3600/dh)+1-wbase*Fs*60 (eh-1)*(Fs*3600/dh)+(Fs*3600/dh)-wbase*Fs*60 87 101]) % one minut
