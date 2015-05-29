
function varargout = Peak_Detection(varargin)
% This program enables fast and accurate detections of synaptic events.
% Please refer to help.chm in the program folder to for further
% introduction.
% Developed by Yulin Shi in Dr. Xiangmin Xu's lab, Dept of Anatomy and
% Neurobiology, UC Irvine
% Last Modified by GUIDE v2.5 25-Mar-2011 16:15:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Peak_Detection_OpeningFcn, ...
    'gui_OutputFcn',  @Peak_Detection_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


% End initialization code - DO NOT EDIT


% --- Executes just before Peak_Detection is made visible.
function Peak_Detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Peak_Detection (see VARARGIN)
global SampleRate
global DataPath
global DataRange
global OnsetPoint
global PlotDataRange
global DirectResponseAnalysisWindow
global DirectResponseAnalysisWindow_2nd

% Choose default command line output for Peak_Detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

[pathstr, name, ext] = fileparts(which('Peak_Detection.m')) ;
temp=load ([pathstr '\Parameter.mat'],'-mat');
%
try
    SampleRate=temp.SampleRate;
    DataPath=temp.DataPath;
    DataRange=temp.DataRange;
    PlotDataRange=(temp.PlotDataRange);
    OnsetPoint=temp.OnsetPoint;
    DirectResponseAnalysisWindow=temp.DirectResponseAnalysisWindow;
    DirectResponseAnalysisWindow_2nd=temp.DirectResponseAnalysisWindow_2nd;
catch e
    SampleRate=10000;
    DataPath='data.ephys.trace_1';
    DataRange=1000:2600;
    PlotDataRange=1:4000;
    OnsetPoint=1000;
    DirectResponseAnalysisWindow=100;
    DirectResponseAnalysisWindow_2nd=300;
   
end
set(handles.edit_Data_Range,'String',[num2str(DataRange(1)/SampleRate*1000) ':' num2str(DataRange(end)/SampleRate*1000)]);
set(handles.edit_Onset, 'String', [num2str(OnsetPoint*1000/SampleRate)]);
set(handles.edit_Direct_Analysis_Window, 'String', [num2str(DirectResponseAnalysisWindow*1000/SampleRate)]);
set(handles.edit_Plot_Window, 'String', [num2str(double(PlotDataRange(1))*1000/SampleRate) ':' num2str(double(PlotDataRange(end))/SampleRate*1000)]);
set(handles.edit_Direct_Analysis_Window,'String',num2str(DirectResponseAnalysisWindow*1000/SampleRate));
set(handles.edit_2nd_Direct_Rsps_Window,'String',num2str(DirectResponseAnalysisWindow_2nd*1000/SampleRate));


% --- Outputs from this function are returned to the command line.
function varargout = Peak_Detection_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_Load.
function pushbutton_Load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Original_Trace
global Original_Trace_before_filter
global Original_Control_Trace
global Control_Trace
global DataPath

global SampleRate
global MapPattern
global xSpacing
global ySpacing
global soma1Coordinates
global spatialRotation
global xPatternOffset
global yPatternOffset

global DataRange
global OnsetPoint
global b_outward
global myData
global xlocation
global ylocation

[fileName,pathName]=uigetfile('.xsg', 'Load xsg file');

if isequal(fileName,0) || isequal(pathName,0)
    return
else
    set(handles.text31,'string','Loading...');
    drawnow update
    temp=load (fullfile(pathName,fileName),'-mat');
    SampleRate=temp.header.ephys.ephys.sampleRate;
    
    eval(['DataRange=' get(handles.edit_Data_Range,'String') ';'] );
    DataRange=DataRange(1)*SampleRate/1000:1:DataRange(end)*SampleRate/1000;
    OnsetPoint=str2num(get(handles.edit_Onset,'String'))*SampleRate/1000;
    
    save_parameter();
    Original_Trace=eval(['temp','.',DataPath]);
    MapPattern=[];
    MapPattern=temp.header.mapper.mapper.mapPatternArray;
    xSpacing=temp.header.mapper.mapper.xSpacing;
    ySpacing=temp.header.mapper.mapper.ySpacing;
    spatialRotation=temp.header.mapper.mapper.spatialRotation;
    xPatternOffset=temp.header.mapper.mapper.xPatternOffset;
    yPatternOffset=temp.header.mapper.mapper.yPatternOffset;
    try
        soma1Coordinates=somaPositionTransformer_Peak_Detection(temp.header.mapper.mapper.soma1Coordinates(1),temp.header.mapper.mapper.soma1Coordinates(2),...
            temp.header.mapper.mapper.spatialRotation, temp.header.mapper.mapper.xPatternOffset, temp.header.mapper.mapper.yPatternOffset);
    catch e
        disp('No Soma Found!');
    end
    trace_length=temp.header.ephys.ephys.traceLength*temp.header.ephys.ephys.sampleRate;
    
    Original_Trace=reshape(Original_Trace,trace_length,[]);
    Original_Trace_before_filter=Original_Trace;
    if get(handles.checkbox_Outward,'value')
        b_outward=1;
    else
        b_outward=0;
    end
    Original_Trace=NoiseRemoval(handles,Original_Trace_before_filter,SampleRate,OnsetPoint);
end
set(handles.text31,'string','Loading Complete');


% --- Executes on button press in pushbutton_Plot.
function pushbutton_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global myData
global Original_Trace
global MapPattern
global SampleRate
global PlotDataRange
global OnsetPoint

arrayTracesAsInputMap_peak_detection(Original_Trace,MapPattern,SampleRate,1,get(handles.checkbox_No_Color_No_Text,'Value'),PlotDataRange,OnsetPoint);


% --- Executes on button press in pushbutton_Train.
function pushbutton_Train_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Original_Trace
global Control_Trace
global PlotDataRange
global training_order

global NUM
global DEN
global K
global Tc
global To
global traceNumber
global Model
global Convolution
global cellcount
global DirectResponseAnalysisWindow
global diff2;
global leads
global trailings
global OnsetPoint
global TotalTrained
global StatusBoxHandle
global Coefficient
global orig_training_data
global Phase_Diff
StatusBoxHandle=handles.text31;

button = questdlg(['All existing filter will be removed, Continue?'],'Continue?', 'Yes','No','No');
if strcmp(button,'No')
    return
end

TotalTrained=0;
[pathstr, name, ext] = fileparts(which('Peak_Detection.m')) ;


temp=load ([pathstr '\Filter.mat'],'-mat');
NUM=temp.Num;
DEN=temp.Den;
K=[];
Tc=[];
To=[];
traceNumber=[];Coefficient=[];
orig_training_data={};
Model={};
Convolution={};
cellcount=1;
diff2=[];
leads=[];
trailings=[];
Phase_Diff=[];
total_trace=size(Original_Trace,2);
Train_Trace=[];
Train_Trace=Original_Trace;
training_order=str2num(get(handles.edit_Training_Order,'string'));

for i= 1:25:total_trace
    if ((total_trace-i)<25)
        figure
        set(gcf, 'Name',['Trace ', num2str(i), ':' num2str(total_trace)]);
        for j=1:(total_trace+1-i)
            subplot(5,5,j);
            text(0,55,num2str(i+j-1));
            hold on
            plot(round(0.75*OnsetPoint):PlotDataRange(end),Train_Trace(round(0.75*OnsetPoint):PlotDataRange(end),i+j-1)-Train_Trace(round(0.75*OnsetPoint),i+j-1));
            set(gca,'ButtonDownFcn','plotdirecttrace');
            set(gca,'UserData',i+j-1);
            set(gca,'YLim', [-400, 50]);
            set(gca,'XLim',[round(0.75*OnsetPoint),PlotDataRange(end)]);
            line([OnsetPoint+DirectResponseAnalysisWindow,OnsetPoint,OnsetPoint,OnsetPoint+DirectResponseAnalysisWindow,OnsetPoint+DirectResponseAnalysisWindow],[50,50,-400,-400,50],'color','r');
        end
    else
        figure,
        set(gcf, 'Name',['Trace ', num2str(i), ':',num2str(i+24)]);
        for j=1:25
            subplot(5,5,j);
            text(0,55,num2str(i+j-1));
            hold on
            plot(round(0.75*OnsetPoint):PlotDataRange(end),Train_Trace(round(0.75*OnsetPoint):PlotDataRange(end),i+j-1)-Train_Trace(round(0.75*OnsetPoint),i+j-1));
            set(gca,'ButtonDownFcn','plotdirecttrace');
            set(gca,'UserData',i+j-1);
            set(gca,'YLim', [-400, 50]);
            set(gca,'XLim',[round(0.75*OnsetPoint),PlotDataRange(end)]);
            line([OnsetPoint+DirectResponseAnalysisWindow,OnsetPoint,OnsetPoint,OnsetPoint+DirectResponseAnalysisWindow,OnsetPoint+DirectResponseAnalysisWindow],[50,50,-400,-400,50],'color','r');
        end
    end
end

% --- Executes on button press in pushbutton_Save.
function pushbutton_Save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global traceNumber
global Model
global Convolution
global Coefficient
global filter_type
global Original_Trace
global diff2
global leads
global trailings
global orig_training_data
global Phase_Diff

[fileName,pathName]=uiputfile('.mat', 'Create data m-file:');
if isequal(fileName,0) || isequal(pathName,0)
    return
end
savetofilename=fullfile(pathName,fileName);
filter_type='poly';
save(fullfile(pathName,fileName),'Phase_Diff','orig_training_data','filter_type','traceNumber','Model','Original_Trace','Convolution','Coefficient','diff2','leads','trailings');
set(handles.text31,'string','Filter Saved');



% --- Executes on button press in pushbutton_View.
function pushbutton_View_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_View (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Coefficient
global Original_Trace
global Control_Trace
global Model
global cellcount
global DataRange
global SampleRate
global diff2
global TH
global RH
global leads
global trailings
global OnsetPoint
global PlotDataRange
global myData
global orig_training_data
global DirectResponseAnalysisWindow
global DirectResponseAnalysisWindow_2nd
global amp_th
global Spon_Amp_Th
global Phase_Diff


View_Trace=Original_Trace;

LINECOLOR=lines(8);
DataLength=length(DataRange);
DataRange_MS=double(PlotDataRange)/(SampleRate/1000.0);

sizeC=length(Coefficient);
i=str2num(get(handles.edit_number,'string'));
seq=[];
fitting={};
final_seq=[];

h=figure;
hold on

temp=load ('Filter.mat','-mat');
NUM=temp.Num;
DEN=temp.Den;
sf = filtfilt(NUM,DEN,View_Trace(:,i));
plot(DataRange_MS,View_Trace(PlotDataRange,i),'k');
set(gca,'XLim',[DataRange_MS(1)-1,DataRange_MS(end)]);
xlabel('ms')
for j=1: cellcount-1
    
    m_dmean=Model{j}-mean(Model{j});
    Filter = m_dmean/sum(abs(m_dmean));
    temp=conv(sf,Filter);
    temp_R=flipud(conv(flipud(sf),Filter));
    fitting{j}=temp((1+round(length(Model{j})/2)+Phase_Diff(j)):(size(View_Trace,1)+round(length(Model{j})/2)+Phase_Diff(j)));
    fitting{j}=fitting{j}/2+fitting{j}/2;
    plot(DataRange_MS,fitting{j}(PlotDataRange),'g');
    plot([DataRange_MS(1),DataRange_MS(end)],[TH,TH],'color',LINECOLOR(2,:));
    plot([DataRange_MS(1),DataRange_MS(end)],[RH,RH],'color',LINECOLOR(4,:));
    if j==1
        seq=(fitting{j}>TH);
        seq_tmp{1}=seq;
    else
        tmp=(fitting{j}>TH);
        seq=seq+tmp;
        seq_tmp{j}=tmp;
    end
end
final_seq=(seq>0);
final_seq2=(final_seq(2:end)-final_seq(1:end-1));
final_seq2=(final_seq2>0);
final_seq2(end+1)=0;

Spontaneous_TE=[];
if ~isempty(final_seq)
    
    TE=[];
    start_stop_pair=[];

    for TE_length=1:length(seq_tmp)
        [tmp_TE,tmp_start_stop_pair]=parse_centermass(seq_tmp{TE_length},10,[0.5 (min(leads)+min(trailings))*0.5],fitting{TE_length});
        TE=[TE,tmp_TE];
        start_stop_pair=[start_stop_pair;tmp_start_stop_pair];
    end
    [TE,i_TE]=sort(TE);
    start_stop_pair=start_stop_pair(i_TE,:);
    tmpcount=size(start_stop_pair,1);
    m=1;
    for j=1:tmpcount
        if TE(m)<DataRange(1)
            Spontaneous_TE=[Spontaneous_TE,TE(m)];
            TE(m)=[];
            continue;
        end
        if TE(m)>(OnsetPoint+DirectResponseAnalysisWindow_2nd)
            break
        end
        tmpdata=zeros(start_stop_pair(j,2)-start_stop_pair(j,1)+1,1);
        for k=1: cellcount-1
            tmpfitting=fitting{k};
            tmpdata=tmpdata+(tmpfitting(start_stop_pair(j,1):start_stop_pair(j,2))>RH);
        end
        if sum(tmpdata)>0
            TE(m)=[];
        else
            m=m+1;
        end
    end
    
end
AA=Lines(30);

set(gcf, 'Name',['Trace ', num2str(i)]);

vara=Align_Peaks( View_Trace(:,i),DataRange,TE,leads,trailings,SampleRate/1000,RH,i,TH,fitting,diff2,OnsetPoint,DirectResponseAnalysisWindow,DirectResponseAnalysisWindow_2nd,amp_th);

spon_vara=Align_Peaks(View_Trace(:,i),1:(OnsetPoint-2),Spontaneous_TE,leads,trailings,SampleRate/1000,0,i,TH,fitting,diff2,OnsetPoint,DirectResponseAnalysisWindow,DirectResponseAnalysisWindow_2nd,amp_th,1,Spon_Amp_Th);
max_Trailings=round(min((max(trailings)),mean(trailings)+std(trailings))*10);

startIndex=vara{1};
endIndex=vara{2};
outIndex=vara{3};
amplitude=vara{4};
PEAK_AMPLITUDE=vara{5};

%%%%%% plot detection information and extracted spikes
figure(h)
if isempty(outIndex)
    return
end
plot(outIndex/(SampleRate/1000), -20,'+k');
figure
hold on
set(gca,'XLim',[DataRange_MS(1)-1,DataRange_MS(end)]);
xlabel('ms')
p=1;
for t= 1: length(startIndex)
    plot([startIndex(t):endIndex(t)]/10,View_Trace(startIndex(t):endIndex(t),i), 'Color',AA(t,:));
    text(outIndex(t)/(SampleRate/1000),View_Trace(DataRange(1)+outIndex(t),i)-10,num2str(PEAK_AMPLITUDE(t)));
    p=p+1;
end
set(gcf, 'Name',['View ', num2str(i)]);
plot(outIndex/(SampleRate/1000), -20,'+k');
peak_NUMBER=length(vara{1});
PEAK_TIME=vara{3}/(SampleRate/1000);
PEAK_AMPLITUDE=vara{5};
AMP_OVER_TIME=sum(vara{4})/(DataRange_MS(end)-DataRange_MS(1));



function edit_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_number as text
%        str2double(get(hObject,'String')) returns contents of edit_number as a double


% --- Executes during object creation, after setting all properties.
function edit_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SampleRate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SampleRate as text
%        str2double(get(hObject,'String')) returns contents of edit_SampleRate as a double
global SampleRate

SampleRate=str2num(get(handles.edit_SampleRate,'String'));
save_parameter();

% --- Executes during object creation, after setting all properties.
function edit_SampleRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SampleRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_DataPath_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DataPath as text
%        str2double(get(hObject,'String')) returns contents of edit_DataPath as a double
global DataPath
DataPath=get(handles.edit_DataPath,'String');
save_parameter()

% --- Executes during object creation, after setting all properties.
function edit_DataPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DataPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_Load_Filter.
function pushbutton_Load_Filter_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Load_Filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileName=[];
[fileName,pathName]=uigetfile('.mat', 'Load Filter file');
if isequal(fileName,0)
    return
end
global K
global Tc
global To
global traceNumber
global Model
global Convolution
global Original_Trace
global Coefficient
global exp_fit
global cellcount
global meanC
global stdC

global TH
global RH
global diff2
global leads
global trailings
global orig_training_data
global amp_th
global Spon_Amp_Th
global meanAmp
global stdAmp
global Phase_Diff

if isequal(fileName,0) || isequal(pathName,0)
    return
end

temp=load (fullfile(pathName,fileName),'-mat');
global filter_type
filter_type=temp.filter_type;

if strcmp(filter_type,'exp')
    K=temp.K;
    Tc=temp.Tc;
    To=temp.To;
end

traceNumber=temp.traceNumber;
Model=temp.Model;
Convolution=temp.Convolution;
Coefficient=temp.Coefficient;
diff2=temp.diff2;
cellcount=length(traceNumber)+1;
meanC=mean(Coefficient);
stdC=std(Coefficient);
leads=temp.leads;
trailings=temp.trailings;

orig_training_data=temp.orig_training_data;
Phase_Diff=temp.Phase_Diff;
orig_amp=[];
set(handles.text31,'string',sprintf('Total %d filters',length(orig_training_data)));
for i_orig=1:length(orig_training_data)
    [p,i_o]=min(orig_training_data{i_orig});
    [p1,i_1]=max(orig_training_data{i_orig}(1:i_o-1));
    orig_amp(i_orig)=abs(p1-p);
end
amp_th=0;
amp_th_i=1;
meanAmp=mean(orig_amp);
stdAmp=std(orig_amp);
while amp_th<=0
    amp_th=meanAmp-5*stdAmp/amp_th_i;
    amp_th_i=amp_th_i*2;
end
Spon_Amp_Th=amp_th;
set(handles.edit_RH_MeanC, 'String','1');
set(handles.edit17, 'String','-1.2');
set(handles.edit14, 'String',num2str(meanC-1.2*stdC));
TH=meanC-1.2*stdC;
set(handles.edit19, 'String','1');
set(handles.edit16, 'String','4');
set(handles.edit15, 'String',num2str(meanC+4*stdC));
RH=meanC+4*stdC;
set(handles.edit_Amp_Th_Mean, 'String','1');
set(handles.edit_Amp_Th_Std, 'String',num2str(-5/amp_th_i*2));
set(handles.edit_Amp_Th_Value, 'String',num2str(amp_th));
set(handles.edit_Spon_Amp_Th, 'String',num2str(Spon_Amp_Th));

% --- Executes on button press in pushbutton_Load_Conf.
function pushbutton_Load_Conf_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Load_Conf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SampleRate
global DataPath
global DataRange
global OnsetPoint
global PlotDataRange
global DirectResponseAnalysisWindow
global DirectResponseAnalysisWindow_2nd

SampleRate=10000;
DataPath='data.ephys.trace_1';
DataRange=1000:2600;
OnsetPoint=1000;
PlotDataRange=1:4000;
DirectResponseAnalysisWindow=100;
DirectResponseAnalysisWindow_2nd=300;
set(handles.edit_Data_Range,'String',[num2str(DataRange(1)/SampleRate*1000) ':' num2str(DataRange(end)/SampleRate*1000)]);
set(handles.edit_Onset,'String',[num2str(OnsetPoint/SampleRate*1000)]);
set(handles.edit_Plot_Window,'String',[num2str(PlotDataRange(1)/SampleRate*1000) ':' num2str(ceil(PlotDataRange(end)/SampleRate*1000))]);
set(handles.edit_Direct_Analysis_Window,'String',num2str(DirectResponseAnalysisWindow/SampleRate*1000));
set(handles.edit_2nd_Direct_Rsps_Window,'String',num2str(DirectResponseAnalysisWindow_2nd/SampleRate*1000));
save_parameter();

% --- Executes on button press in pushbutton_Save_Fonf.
function pushbutton_Save_Fonf_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Save_Fonf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edit_Data_Range_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Data_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Data_Range as text
%        str2double(get(hObject,'String')) returns contents of edit_Data_Range as a double
global DataRange
global SampleRate
global DataPath
global OnsetPoint
global PlotDataRange
global DirectResponseAnalysisWindow
eval(['DataRange=' get(handles.edit_Data_Range,'String') ';']);
DataRange=DataRange(1)*SampleRate/1000:1:DataRange(end)*SampleRate/1000;
save_parameter()

% --- Executes during object creation, after setting all properties.
function edit_Data_Range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Data_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_Load_Default.
function pushbutton_Load_Default_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Load_Default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global K
global Tc
global To
global traceNumber
global Model
global Convolution
global Original_Trace
global Coefficient
global exp_fit
global cellcount
global meanC
global stdC
global TH
global RH
global diff2
global leads
global trailings
temp=load ('c_Filter.mat','-mat');
global filter_type
filter_type=temp.filter_type;

if strcmp(filter_type,'exp')
    K=temp.K;
    Tc=temp.Tc;
    To=temp.To;
else
end
leads=temp.leads;
trailings=temp.trailings;
traceNumber=temp.traceNumber;
Model=temp.Model;
Convolution=temp.Convolution;
Coefficient=temp.Coefficient;
diff2=temp.diff2;
cellcount=length(traceNumber)+1;
meanC=mean(Coefficient);
stdC=std(Coefficient);
set(handles.edit_RH_MeanC, 'String','1');
set(handles.edit17, 'String','-1');
set(handles.edit14, 'String',num2str(meanC-stdC));

set(handles.edit19, 'String','1');
set(handles.edit16, 'String','3');
set(handles.edit15, 'String',num2str(meanC+3*stdC));
TH=meanC-stdC;
RH=meanC+3*stdC;

% --- Executes on button press in pushbutton_Detect_All.
function pushbutton_Detect_All_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Detect_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global Coefficient
global Original_Trace
global Control_Trace
global Model
global cellcount
global DataRange
global SampleRate
global TH
global RH
global diff2
global leads
global trailings
global OnsetPoint
global soma1Coordinates
global MapPattern
global xSpacing
global ySpacing
global PlotDataRange
global b_outward
global xlocation
global ylocation
global spatialRotation
global xPatternOffset
global yPatternOffset
global DirectResponseAnalysisWindow
global DirectResponseAnalysisWindow_2nd
global amp_th
global Spon_Amp_Th
global Phase_Diff
global excluded_Human_Pick_rsps_Data
global Human_Pick_rsps_Data
[fileName,pathName]=uiputfile('.mat', 'Create data m-file:');
if isequal(fileName,0)
    return
end
eval(['Trace_To_Exclude=[' get(handles.edit_trace_Exclusion,'string') '];']);
DataLength=length(DataRange);
DataRange_MS=DataRange/(SampleRate/1000);
final_seq={};
CT_final_seq={};
% i=95

All_Detection_Results={};
NO_Results=1;
seq=[];
CT_seq=[];
fitting={};

%%%%%Load High pass filter parameters
temp=load ('Filter.mat','-mat');
NUM=temp.Num;
DEN=temp.Den;


%%%%%Initialization
peak_NUMBER=[];
PEAK_TIME={};
PEAK_AMPLITUDE={};
PEAK_AMPLITUDE_FOR_HIST=[];
AMP_OVER_TIME=[];
SUM_AMPLITUDE={};
ORIG_TRACE={};
PLOT_ENDINDEX={};
FITTING_DATA={};
SPON_PEAK_AMPLITUDE_FOR_HIST=[];
Spon_AMP_OVER_TIME={};
Spon_PEAK_AMPLITUDE={};
Spon_AllDetectionResults={};
Spon_peak_NUMBER=[];
Spon_PEAK_TIME={};
low_th_num=str2num(get(handles.edit17,'string'));
human_pick_rsps_trace=[];
total_peak_offset_trace={};
total_bias_trace={};
total_count_corr_trace={};
total_out_corrected_index_trace={};

%%%%%%%
for i =1 : size(Original_Trace,2)
    sf = filtfilt(NUM,DEN,Original_Trace(:,i));
    total_peak_offset=[];
    final_seq=[];
    seq_tmp={};
    %%%%%% Convolution
    for j=1: cellcount-1
        m_dmean=Model{j}-mean(Model{j});
        Filter = m_dmean/sum(abs(m_dmean));
        temp=conv(sf,Filter);
        temp_R=flipud(conv(flipud(sf),Filter));
        fitting{j}=temp((1+round(length(Model{j})/2)+Phase_Diff(j)):(size(Original_Trace,1)+round(length(Model{j})/2)+Phase_Diff(j)));
        fitting_r{j}=temp_R(1+round(length(Model{j})/2):size(Original_Trace,1)+round(length(Model{j})/2));
        
        % % % % % % % %
        
        %%%%%%% Get suprathreholds
        if j==1
            seq=(fitting{j}>TH);
            seq_tmp{j}=seq;
        else
            tmp=(fitting{j}>TH);
            seq=seq+tmp;
            seq_tmp{j}=tmp;
        end
        
    end
    final_seq=(seq>0);
    
    %%%%%%%
    %%%%%Get all all mass center, then sort in acsending order by the
    %%%%%arrival time
    TE=[];
    start_stop_pair=[];
    abcabcabc=0;
    for TE_length=1:length(seq_tmp)
        if isempty(seq_tmp{TE_length})
            continue
        end
        abcabcabc=1;
        [tmp_TE,tmp_start_stop_pair]=parse_centermass(seq_tmp{TE_length},10,[0.5 (min(leads)+min(trailings))*0.5],fitting{TE_length});
        
        TE=[TE,tmp_TE];
        start_stop_pair=[start_stop_pair;tmp_start_stop_pair];
    end
    if ~abcabcabc
        continue
    end
    [TE,i_TE]=sort(TE);
    start_stop_pair=start_stop_pair(i_TE,:);
    
 
    
    tmpcount=size(start_stop_pair,1);
    m=1;
    Spontaneous_TE=[];
    
    %%%%%%%%%%%%%%%%%%%%%%Peak detection%%%%%%%%%%
    
    for j=1:tmpcount
        if TE(m)<DataRange(1)
            Spontaneous_TE=[Spontaneous_TE,TE(m)];
            TE(m)=[];
            continue
        end
  
        if TE(m)>(OnsetPoint+DirectResponseAnalysisWindow_2nd)
            break
        end
        tmpdata=zeros(start_stop_pair(j,2)-start_stop_pair(j,1)+1,1);
        for k=1: cellcount-1
            tmpfitting=fitting{k};
            
            tmpdata=tmpdata+(tmpfitting(start_stop_pair(j,1):start_stop_pair(j,2))>RH);
        end
        if sum(tmpdata)>0
            TE(m)=[];
        else
            m=m+1;
        end
    end
    
    set(handles.text31,'string',['Processing trace ',num2str(i)]);
    drawnow update
    vara=Align_Peaks(Original_Trace(:,i),DataRange,TE,leads,trailings,SampleRate/1000,RH,i,TH,fitting,diff2,OnsetPoint,DirectResponseAnalysisWindow,DirectResponseAnalysisWindow_2nd,amp_th);
    %%%%%%%%%%%%%%%%%%%%%%%% PEAK DETECTION%%%%%%%%%%%%%%%%%%
    if find(i==Trace_To_Exclude)
        
        vara={};
        
    end
    
    spon_vara=Align_Peaks(Original_Trace(:,i),1:(OnsetPoint-2),Spontaneous_TE,leads,trailings,SampleRate/1000,0,i,TH,fitting,diff2,OnsetPoint,DirectResponseAnalysisWindow,DirectResponseAnalysisWindow_2nd,amp_th,1,Spon_Amp_Th);
    
    
    
    FITTING_DATA{NO_Results}=fitting;
    
    vara{6}=i;
    All_Detection_Results{NO_Results}=vara;
    peak_NUMBER=[peak_NUMBER,length(vara{1})];
    PEAK_TIME{NO_Results}=vara{3}/(SampleRate/1000);
  
    PEAK_AMPLITUDE_FOR_HIST=[PEAK_AMPLITUDE_FOR_HIST,vara{5}];
    AMP_OVER_TIME{NO_Results}=sum(vara{4})/(DataRange_MS(end)-DataRange_MS(1));
    SUM_AMPLITUDE{NO_Results}=vara{4};
    
    SPON_PEAK_AMPLITUDE_FOR_HIST=[SPON_PEAK_AMPLITUDE_FOR_HIST,spon_vara{5}];
    spon_vara{6}=i;
    Spon_AllDetectionResults{NO_Results}=spon_vara;
    
    Spon_AMP_OVER_TIME{NO_Results}=sum(spon_vara{4})/(OnsetPoint*1000/SampleRate);
    
    Spon_PEAK_AMPLITUDE{NO_Results}=spon_vara{5};
    Spon_peak_NUMBER=[Spon_peak_NUMBER , length(spon_vara{1})];
    Spon_PEAK_TIME{NO_Results}=spon_vara{3}/(SampleRate/1000);
    NO_Results=NO_Results+1;
end
set(handles.text31,'string','Done!');
Spon_DataRange=1:(OnsetPoint-2);


%%%Save all the detection results and parameters.
save([pathName,fileName],...
    'Spon_peak_NUMBER',...
    'Spon_AMP_OVER_TIME',...
    'Spon_AllDetectionResults',...
    'Spon_DataRange',...
    'DataRange',...
    'PlotDataRange', ...
    'soma1Coordinates',...
    'SampleRate',...
    'OnsetPoint',...
    'peak_NUMBER',...
    'AMP_OVER_TIME',...
    'All_Detection_Results',...
    'Original_Trace',...
    'low_th_num',...
    'b_outward',...
    'MapPattern',...
    'xSpacing',...
    'ySpacing',...
    'spatialRotation',...
    'xPatternOffset',...
    'yPatternOffset',...
    'Model',...%%%%%%Conv
    'Phase_Diff'...%Conv
    );
a=figure;
set(a,'name','PEAK AMPLITUDE');
hist(PEAK_AMPLITUDE_FOR_HIST);
b=figure;
set(b,'name','SPONTANEOUS PEAK AMPLITUDE');
hist(SPON_PEAK_AMPLITUDE_FOR_HIST);




function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double
global TH
global meanC
global stdC
TH = str2num(get(handles.edit_RH_MeanC, 'String'))*meanC+str2num(get(handles.edit17, 'String'))*stdC;
set(handles.edit14, 'String',num2str(TH));

% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
global RH
global meanC
global stdC
RH = str2num(get(handles.edit19, 'String'))*meanC+str2num(get(handles.edit16, 'String'))*stdC;
set(handles.edit15, 'String',num2str(RH));

% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_RH_MeanC_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RH_MeanC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_RH_MeanC as text
%        str2double(get(hObject,'String')) returns contents of edit_RH_MeanC as a double
global TH
global meanC
global stdC
TH = str2num(get(handles.edit_RH_MeanC, 'String'))*meanC+str2num(get(handles.edit17, 'String'))*stdC;
set(handles.edit14, 'String',num2str(TH));


function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double
global TH
TH = str2num(get(handles.edit14, 'String'));



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double
global RH
global meanC
global stdC
RH = str2num(get(handles.edit19, 'String'))*meanC+str2num(get(handles.edit16, 'String'))*stdC;
set(handles.edit15, 'String',num2str(RH));



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double
global RH
RH = str2num(get(handles.edit15, 'String'));



function edit_Onset_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Onset as text
%        str2double(get(hObject,'String')) returns contents of edit_Onset as a double
global DataRange
global SampleRate
global DataPath
global OnsetPoint
global PlotDataRange
global DirectResponseAnalysisWindow
OnsetPoint=str2num(get(handles.edit_Onset,'String'))*SampleRate/1000;
save_parameter()

% --- Executes during object creation, after setting all properties.
function edit_Onset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileName=[];
[fileName,pathName]=uigetfile('.mat', 'Load file');
if isequal(fileName,0)
    return
end
temp=load (fullfile(pathName,fileName),'-mat');
datavector=[];
for i = 1:length(temp.PEAK_AMPLITUDE)
    if isempty(temp.PEAK_AMPLITUDE{i})
        datavector(i)=0;
    else
        datavector(i)=mean(temp.PEAK_AMPLITUDE{i});
    end
end
plotcolormap(datavector,temp.MapPattern,temp.xSpacing, temp.ySpacing,temp.soma1Coordinates,'Peak Amplitude');


% --- Executes on button press in pushbutton_Save_Results.
function pushbutton_Save_Results_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Save_Results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileName=[];
[fileName,pathName]=uigetfile('.mat', 'Load file');
if isequal(fileName,0)
    return
end
temp=load (fullfile(pathName,fileName),'-mat');
datavector=[];
for i = 1:length(temp.SUM_AMPLITUDE)
    if isempty(temp.SUM_AMPLITUDE{i})
        datavector(i)=0;
    else
        datavector(i)=sum(temp.SUM_AMPLITUDE{i});
    end
end
plotcolormap(datavector,temp.MapPattern,temp.xSpacing, temp.ySpacing,temp.soma1Coordinates,'Sum Amplitude');

% --- Executes on key press with focus on pushbutton12 and none of its controls.
function pushbutton12_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_Freq.
function pushbutton_Freq_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileName=[];
[fileName,pathName]=uigetfile('.mat', 'Load file');
if isequal(fileName,0)
    return
end
temp=load (fullfile(pathName,fileName),'-mat');

plotcolormap(temp.peak_NUMBER,temp.MapPattern,temp.xSpacing, temp.ySpacing,temp.soma1Coordinates,'Freq');


% --- Executes on button press in pushbutton_Onset.
function pushbutton_Onset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileName=[];
[fileName,pathName]=uigetfile('.mat', 'Load file');
if isequal(fileName,0)
    return
end
temp=load (fullfile(pathName,fileName),'-mat');
datavector=[];
for i = 1:length(temp.PEAK_TIME)
    if isempty(temp.PEAK_TIME{i})
        datavector(i)=inf;
    else
        TEMP= temp.PEAK_TIME{i};
        datavector(i)=TEMP(1);
    end
end
temp.soma1Coordinates
plotcolormap(datavector-100,temp.MapPattern,temp.xSpacing, temp.ySpacing,temp.soma1Coordinates,'Onset');


% --- Executes during object creation, after setting all properties.
function edit_RH_MeanC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RH_MeanC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Plot_Window_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Plot_Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Plot_Window as text
%        str2double(get(hObject,'String')) returns contents of edit_Plot_Window as a double
global DataRange
global PlotDataRange
global SampleRate
global DataPath
global OnsetPoint
global DirectResponseAnalysisWindow
str= get(handles.edit_Plot_Window,'string');
colonindex=strfind(str, ':');
indexbegin=str2num(str(1:colonindex-1));
indexend=str2num(str(colonindex+1:end));
PlotDataRange=( indexbegin:1000/SampleRate:indexend)*SampleRate/1000;
PlotDataRange=int32(PlotDataRange);
save_parameter()

% --- Executes during object creation, after setting all properties.
function edit_Plot_Window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Plot_Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_No_Color_No_Text.
function checkbox_No_Color_No_Text_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_No_Color_No_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_No_Color_No_Text


% --- Executes on button press in checkbox_Outward.
function checkbox_Outward_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Outward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Outward
global Original_Trace
if ~isempty(Original_Trace)
    Original_Trace=-Original_Trace;
end


% --- Executes on button press in checkbox_load_old_codes.
function checkbox_load_old_codes_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_load_old_codes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_load_old_codes


% --- Executes during object creation, after setting all properties.
function checkbox_load_old_codes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkbox_load_old_codes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox_control_trace.
function checkbox_control_trace_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_control_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_control_trace



function edit_Direct_Analysis_Window_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Direct_Analysis_Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Direct_Analysis_Window as text
%        str2double(get(hObject,'String')) returns contents of edit_Direct_Analysis_Window as a double
global DataRange
global SampleRate
global DataPath
global OnsetPoint
global PlotDataRange
global DirectResponseAnalysisWindow
DirectResponseAnalysisWindow=str2num(get(handles.edit_Direct_Analysis_Window,'String'))*SampleRate/1000;
save_parameter();

% --- Executes during object creation, after setting all properties.
function edit_Direct_Analysis_Window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Direct_Analysis_Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Training_Order_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Training_Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Training_Order as text
%        str2double(get(hObject,'String')) returns contents of edit_Training_Order as a double


% --- Executes during object creation, after setting all properties.
function edit_Training_Order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Training_Order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_2nd_Direct_Rsps_Window_Callback(hObject, eventdata, handles)
% hObject    handle to edit_2nd_Direct_Rsps_Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_2nd_Direct_Rsps_Window as text
%        str2double(get(hObject,'String')) returns contents of edit_2nd_Direct_Rsps_Window as a double
global DataRange
global SampleRate
global DataPath
global OnsetPoint
global PlotDataRange
global DirectResponseAnalysisWindow
global DirectResponseAnalysisWindow_2nd
DirectResponseAnalysisWindow_2nd=str2num(get(handles.edit_2nd_Direct_Rsps_Window,'String'))*SampleRate/1000;
save_parameter();


% --- Executes during object creation, after setting all properties.
function edit_2nd_Direct_Rsps_Window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_2nd_Direct_Rsps_Window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function save_parameter()
global DataRange
global SampleRate
global DataPath
global OnsetPoint
global PlotDataRange
global DirectResponseAnalysisWindow
global DirectResponseAnalysisWindow_2nd
[pathstr, name, ext] = fileparts(which('Peak_Detection.m')) ;
save ([pathstr '\Parameter.mat'],'SampleRate','DataPath','DataRange','OnsetPoint','PlotDataRange','DirectResponseAnalysisWindow','DirectResponseAnalysisWindow_2nd');



function edit_Amp_Th_Mean_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Th_Mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Amp_Th_Mean as text
%        str2double(get(hObject,'String')) returns contents of edit_Amp_Th_Mean as a double

global meanAmp
global stdAmp
global amp_th

amp_th=str2num(get(handles.edit_Amp_Th_Mean,'string'))*meanAmp+str2num(get(handles.edit_Amp_Th_Std,'string'))*stdAmp;
set(handles.edit_Amp_Th_Value,'string',num2str(amp_th));
% --- Executes during object creation, after setting all properties.
function edit_Amp_Th_Mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Th_Mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_Amp_Th_Value_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Th_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Amp_Th_Value as text
%        str2double(get(hObject,'String')) returns contents of edit_Amp_Th_Value as a double
global amp_th
amp_th=str2num(get(handles.edit_Amp_Th_Value,'string'));


% --- Executes during object creation, after setting all properties.
function edit_Amp_Th_Value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Th_Value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Amp_Th_Std_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Th_Std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Amp_Th_Std as text
%        str2double(get(hObject,'String')) returns contents of edit_Amp_Th_Std as a double

global meanAmp
global stdAmp
global amp_th

amp_th=str2num(get(handles.edit_Amp_Th_Mean,'string'))*meanAmp+str2num(get(handles.edit_Amp_Th_Std,'string'))*stdAmp;
set(handles.edit_Amp_Th_Value,'string',num2str(amp_th));

% --- Executes during object creation, after setting all properties.
function edit_Amp_Th_Std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Amp_Th_Std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Spon_Amp_Th_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Spon_Amp_Th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Spon_Amp_Th as text
%        str2double(get(hObject,'String')) returns contents of edit_Spon_Amp_Th as a double
global Spon_Amp_Th
Spon_Amp_Th=str2num(get(handles.edit_Spon_Amp_Th,'string'));

% --- Executes during object creation, after setting all properties.
function edit_Spon_Amp_Th_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Spon_Amp_Th (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_trace_Exclusion_Callback(hObject, eventdata, handles)
% hObject    handle to edit_trace_Exclusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_trace_Exclusion as text
%        str2double(get(hObject,'String')) returns contents of edit_trace_Exclusion as a double


% --- Executes during object creation, after setting all properties.
function edit_trace_Exclusion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_trace_Exclusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Z=NoiseRemoval(handles, Original_Trace_before_filter,SampleRate,OnsetPoint)

Z=Original_Trace_before_filter;

if get(handles.checkbox_Outward,'value')
    Z = filterTraceArray_peak_detection(...
        -Z,SampleRate, ...
        'mean', 11);
else
    Z = filterTraceArray_peak_detection(...
        Z,SampleRate, ...
        'mean', 11);
end
Z = filterTraceArray_peak_detection(...
    Z,SampleRate, ...
    'mean', 11);
Z = baselineSubtractTraceArrayPhys_peak_detection(...
    Z, ...
    1, OnsetPoint-2);
