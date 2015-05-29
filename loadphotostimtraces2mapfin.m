%%%%%%%%%%%%%%%%%%%%%%%
%%%% AutoTraceCoding Series (1) - LoadPhotoStimTraces2MapFfin
%%% Xu
%%Last update 11/17/2005
%% Coding strategy based upon Gordon Sheperd's analysis routines (see
%% communications)

%%% Need data files:  photostim data (*P binary trace file), control data (*P_C
%%% binary trace file), photostim site data (*P.DAT ) and control site data
%%% (*P_C.DAT).  The program will prompt you to get the photostim data (*P
%%% binary trace file), and follow the GUI instructions.

%%% Load trace files, select only related trace data and put them into data
%%% array.  Do baseline subtraction and minimal smoothing.  Plot photostim
%%% and control traces.  Sort photostim traces into direct, indirect and
%%% null traces.      save (savename, 'para', 'constant', 'sumPHAmp_array', 'sumCTAmp_array');



% clear all;
% close all; 


%%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters for trace plotting  %%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

% Condition input arguments
def.stimcnts = '1:332'; %stimulation count series
def.aligncnts = '323:337'; %alignment count series
def.stimdur  = '10';          % laser onset duration to plot in ms
def.stimdel = '50'; % stimulation delay in ms
def.indirectwindow = '150';
def.polarity = '1'; % polarity: excitatory or inhibitory
def = struct2cell (def);    % converts def from structure to cell

%%Input parameters for trace plotting

prompt = {'Enter stimulation count series:', 'Enter alignment count series:', 'Enter laser onset duration in ms',...
        'Enter stim delay in ms', 'Enter indirect window period in ms', 'Enter polarity: excitatory 1, inhibitory 2'};
dlg_title = 'Enter parameters for trace plotting';
num_lines= 1;
parameters = inputdlg(prompt, dlg_title, num_lines, def);

stimcnts = str2num(char(parameters(1)));      % stimulation count series
aligncnts = str2num(char(parameters(2)));  % alignment count series
stimdur = str2num(char(parameters(3)));       % plot duration in ms
stimdel = str2num(char(parameters(4)));      % stimulation delay in ms
indirectwindow = str2num(char(parameters(5)));   % indirect window in ms
polarity = str2num(char(parameters(6)));       % soma location


%%constants
downsamp = 1;         % factor to downsample (1=none)
srate   = 10000/downsamp; % samples per second, (AD convert: 10 KHz)
framdur = 400;        % duration of frame / block (ms) // the DOS program records up to 400 ms of data for one stimulation

nframsamp = srate*framdur/1000;  % samples in a frame

basesamplepoints_stimdel = stimdel*srate/1000;
laseronpoints = stimdur*srate/1000;%% 10 ms

base_pts = basesamplepoints_stimdel;
baselaser_pts = basesamplepoints_stimdel + laseronpoints;
indirectwindow_pts = indirectwindow*srate/1000; %% normally 140 ms after laser off(10 ms)

alpha   = 20;          % output gain selection (usu 20 for photostimulation trials)
beta    = 1;          % feedback resistor gain selection (1)
maxV    =  10;        % largest voltage value
minV    = -maxV;      % smallest voltage value
maxADC  =  2^11;      % largest quantized value
minADC  = -maxADC;    % smallest quantized value


% Plotting and display constants
xyinc = 84;           % xy increment is 84  (for 50 um)
xyadj = 1.68;         % adjustment for x and y coordinates
xymic  = xyinc / xyadj;  % increment in microns


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Get data files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the photostimulation trial data file 
topdir=pwd;
[filename,pathname] =uigetfile({'*','All Files (*)'},'Open photostimulation data file');  % get data for c-clamp/v-clamp traces

fid=fopen(strcat(pathname,filename),'r');

tracedat = fread(fid, inf, 'short'); %% 'short'  integer,  16 bits.
fclose(fid);
tracedat = tracedat(1:downsamp:end);

myData.seq.load =1;



%%plot all traces at voltage clamp mode here

% Notes:  % Plot voltage pulses during current clamp
% 
%     curVtrace = 10^3 * (curVtrace * maxV) / ...
%       (maxADC * 10 * alpha); % mV
% 
%   % Plot current pulses during voltage clamp
%     curCtrace = curCtrace - median(curCtrace(baselinecnt));
%     curCtrace = 10^9 * (curCtrace * maxV) / ...
%       (maxADC * alpha * beta * 10^9); % nA

PHtrace = 10^3*(tracedat * maxV) / ...
    (maxADC *alpha * beta); % PA   maxV, largest voltage value; maxADC, largest quantized value  


%% Get the control photostimulation data file 
controltracefile   = fullfile(pathname, [filename, '_C']);

fid=fopen(controltracefile,'r');

controltracedat = fread(fid, inf, 'short');
fclose(fid);

controltracedat = controltracedat(1:downsamp:end);

CTtrace = 10^3*(controltracedat * maxV) / ...
    (maxADC *alpha * beta); % PA   maxV, largest voltage value; maxADC, largest quantized value  

%%Load photostim site count file

sitecntfile   = fullfile(pathname, [filename, '.DAT']);

stimsitemat = load(sitecntfile);

%%Load control site count file

ctsitecntfile   = fullfile(pathname, [filename, '_C.DAT']);

ctstimsitemat = load(ctsitecntfile);

controlcountindex = ctstimsitemat(:, 3);


%%%%%The motorized stage produced some slightly different x y coordinates
%%%%%when move back to the location.  The following lines help to get
%%%%%things straight

button = questdlg('Do you want to adjust stimsite coordinates?','Adjust?', 'Yes','No','Yes');

if strcmp(button,'Yes')
    
    
    for i = 1 : length(stimsitemat);
        
        if abs(stimsitemat (i, 1)- round(stimsitemat(i, 1)/84)*84) > 40; %%mod line)
            stimsitemat(i, 1) = round(stimsitemat(i, 1)/84)*84 -84;
        else
            stimsitemat(i, 1) = round(stimsitemat(i, 1)/84)*84;
        end
        
        if abs(stimsitemat (i, 2)- round(stimsitemat(i, 2)/84)*84) > 40; %%mod line)
            stimsitemat(i, 2) = round(stimsitemat(i, 2)/84)*84 -84;
        else
            stimsitemat(i, 2) = round(stimsitemat(i, 2)/84)*84;
        end
        
        
    end
    
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


stimsitemat(:, 1:2) = stimsitemat(:, 1:2)/xyadj;  % 84 for 50 microns    stimstie: x, y, count 


%%only get good or relevant traces for further processing 

goodtraces = find(ismember(stimsitemat(:,3), stimcnts))';

% aligntraces = find(ismember(stimsitemat(:,3), aligncnts))';


%% Put good traces into data arrays

ni = 0;

for i = goodtraces;
    
    ni = ni +1;
    
    x= ((i-1)*nframsamp +1) : i*nframsamp;
    PHarray(:, ni)= PHtrace(x);
    
    %%control count, exactly the control trace (#cni) related to the photo stim trace (#ni)    
    control_id = find(controlcountindex==i);  %%controlcountindex = ctstimsitemat(:, 3);
    control_id = control_id(1);
    
    ct_x= ((control_id-1)*nframsamp +1) : control_id*nframsamp;
    CTarray(:, ni) = CTtrace(ct_x);
    
    stimsiteValid(ni, :) = stimsitemat(i, :);
    
    %     goodtraceind (ni) = i;
    
end

ni = 0;
for i = 1:size(stimsitemat,1);
    
    ni = ni +1;
    
    %%control count, exactly the control trace (#cni) related to the photo stim trace (#ni) 
    if ~isempty(find(controlcountindex==i))
    control_id = find(controlcountindex==i);  %%controlcountindex = ctstimsitemat(:, 3);
    control_id = control_id(1);
    end 
        
    
    ct_x= ((control_id-1)*nframsamp +1) : control_id*nframsamp;
    ct= ((ni-1)*nframsamp +1) : ni*nframsamp;    
    controltracedat_singal(ct) = controltracedat(ct_x);
    
end

PHarray_org = PHarray;
CTarray_org = CTarray;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Process raw traces %%%%%%%%%%%%%%%%%%%%%%%%%%%traces%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rows cols] = size(PHarray);

%%% smoothing (mean filtering) and baseline substraction

%% smoothing -- it's minimal, and gets rid of high freq noises and does not affect EPSC slopes, etc
PHarray = colfilt(PHarray, [11 1], 'sliding', @mean);    
CTarray = colfilt(CTarray, [11 1], 'sliding', @mean);    


%% baselinesubstraction

baselineadjust = median(PHarray(1:basesamplepoints_stimdel, :), 1);

BaseoffsetArray = repmat(baselineadjust, rows, 1);

PHarray = PHarray - BaseoffsetArray;
PHarray4base = PHarray;%% for plotting sorted response and measurement display in Plotdirecttraces
%%in single electro and also calculate in dual e.g. sumPH_baseamp_array = sum (PHarray4base(1 : base_pts, :), 1)';


%%control trace array
ctbaselineadjust = median(CTarray(1:basesamplepoints_stimdel, :), 1);

CT_baseoffsetArray = repmat(ctbaselineadjust, rows, 1);

CTarray = CTarray - CT_baseoffsetArray;


%Mi added
myData.data.analysis.numberOfMaps=1;
myData.sampleRate=srate;  
myData.data.map.map1.tracedat=tracedat;           %Original data from data file
myData.data.map.map1.controltracedat=controltracedat_singal; %Original control data from data file
myData.data.map.map1.bsArray=PHarray;             %baseline substraction traces (PA)
myData.data.map.map1.CTbsArray=CTarray;           %baseline substraction control traces (PA)
myData.data.map.map1.stimsitemat=stimsitemat;     % site info
myData.data.map.map1.goodtraces=goodtraces;       % cnt of traces to be show
myData.data.map.map1.location=stimsiteValid;      % valid location of the traces
myData.pathname=pathname;
myData.filepath=filename;


% generate mappattern for the old data
% stimxdim = length(unique(round(stimsiteValid(:,1))));        %dimension of x coordinates
% stimydim = length(unique(round(stimsiteValid(:,2))));        %dimension of y coordinates
% mappattern=zeros(stimydim,stimxdim,1); 

% end
