function handles = arrayTracesAsInputMap_peak_detection(data,mapPattern,SampleRate,mode,b_plot_color_text,PlotDataRange,OnsetPoint)
% arrayTracesAsInputMap
%
% Plots the set of map traces as a map.
%
% originally arrayTracesAsMap
%
% Editing:
% gs april 2005 -- privatize version for new analysis software 
% -----------------------------------------------------------------

flipFlag = 0;
   
% data
%%%% mode 1 plotting raw data;
%%%%% mode 2 plotting corrected data;
switch mode
    case 1 
    array = data;
    namedisplay = 'arrayTracesAsMap  raw data';
    case 2
    array = data.data.CorrectedRawData;
    namedisplay = 'arrayTracesAsMap  corrected data';
    otherwise
        msgbox('Wrong plotting mode');
        return;
end
traceLength = size(array,1 );
num_of_traces=size(array,2);
if num_of_traces==1
    figure
    set(gcf, 'name','CurrentInjection')
    plot(array(:,1))
    return
end

% arrayinftozero = array;
% arrayinftozero(find(array ==inf))=0;

for n=1:numel(mapPattern)
    newArray(:,find(mapPattern==n)) = array(:,n);
end
% 
% try
%     sr =data.data.acq.sampleRate;
% catch
%     sr = 10000; 
%     disp('Using sample rate of 10K by default');
%     beep
% end
sr=SampleRate;
% startTime = 0.05;
% stopTime = 0.30;
% 
% showStart = 2000;
% showStop = 4000;
showStart = PlotDataRange(1);
showStop = PlotDataRange(end);
% showStart = round(data.data.analysis.traceMapShowStart * sr);
% showStop = round(data.data.analysis.traceMapShowStop * sr);

array = newArray(showStart:showStop, :);
[rows,cols] = size(array);
totalTime = (rows-1)/sr; 
xTimeAxis = linspace(0, totalTime, rows)';
traceAxis = ( 1 : cols );
% quarterPoint = round(cols/4);
% sixteenthPoint = round(cols/16);

% [sizeX, sizeY] = size(state.uncaging.analysis.pulsePattern);
[sizeX, sizeY] = size(mapPattern);

% yFactor = 100; % offset, in pA
yFactor =600;
scaleBarTxt = 'pA';
% if ~state.uncaging.map.cellAttachedCheck
%     yFactor = 100; % 100 pA offset
%     scaleBarTxt = 'pA';
% elseif state.uncaging.map.cellAttachedCheck
%     yFactor = 5; % 10 mV offset 
%     scaleBarTxt = 'mV';
% end
if flipFlag == 1
    yFactor = -yFactor;
end
offsetVector = yFactor * ( 0 : cols-1 );
offsetArray = repmat(offsetVector, rows, 1);
if mode == 2
arrayinftozero = array;
arrayinftozero(find(array ==inf))=0;

arrayinftozero = arrayinftozero-offsetArray;
end
array = array-offsetArray;
[rows,cols] = size(array);
% set up the figure -------------------------------------------------------------
x = .14; 
y = .11; 
w = .5; 
h = .8; 

 figure('Units', 'normalized', ...
    'Position', [x y w h], 'Name', namedisplay, ...
    'NumberTitle', 'off', 'Color', 'w');

% try
%     set(gcf, 'ButtonDownFcn', 'arrayTracesAsMap_tweakFig');
% catch
% end
subplotRows = 1; subplotCols = sizeY; plotnum = 0;
if b_plot_color_text==0
    AA=lines(16);
else
    AA=zeros(16,3);
end
% AA(:,1)=[0,0,0];
for N = 1:sizeY
    startInd = (N-1)*sizeX + 1;
    endInd = N*sizeX;
%     startInd = 15;
    plotnum = plotnum+1;

    %     hsub(plotnum) = subplot(subplotRows,subplotCols,plotnum);
    
    pos1 = 0.025 + (N - 1)*(0.96/sizeY);
    pos2 = 0.02;
    pos3 = 0.05;
    pos4 = 0.96;
    hsub(N) = axes('Position', [pos1 pos2 pos3 pos4]);
    step_color=1/16;
%     set(gca, 'Position', );
abcabc=mapPattern(:,N);
    for M=startInd:endInd
      plot(xTimeAxis, array(:,M),'LineWidth',2,'UserData',{newArray(:,M),abcabc(mod(M-1,sizeX)+1)},'ButtonDownFcn',@plot_small_map_trace,'Color', AA(M-startInd+1,:)); %%%plotting  formatted traces which are modified to fit figure space //differentfrom original traces data.
    hold on;
    end
%%%plotting corresponding time sequence index.
    PosX = ones(sizeX,1)*0.2;
    PosY = [pos2+0.05:pos4/sizeX:(pos2+pos4-0.5*pos4/sizeX+0.05)]';
    if b_plot_color_text==0
        strhandle=text(PosX,PosY,num2str(flipdim(mapPattern(:,N),1)),'Units','normalized');
    end 
    %%% set the ones who were adjusted red.
    tmparray = flipdim(mapPattern(:,N),1);
%     if mode ==2
%         for N1= 1:length(data.data.CorrectLabel)
%             set(strhandle(find(tmparray==data.data.CorrectLabel(N1))),'color','r');
%         end
%     end
    
    
%%%plotting corresponding time sequence index.
    
    
    
    
%     minval = min(mean(array(1:100,startInd:endInd)));
%     maxval = max(mean(array(1:100,startInd:endInd)));
%     tweakFactor = abs(maxval - minval)*0.05;
%     yrange = [minval-tweakFactor maxval+tweakFactor];
    
    if mode ==2
        minval = min(mean(arrayinftozero(1:100,startInd:endInd)));
        maxval = max(mean(arrayinftozero(1:100,startInd:endInd)));
    else
        minval = min(mean(array(:,startInd:endInd)));
        maxval = max(mean(array(:,startInd:endInd)));
    end
    tweakFactor = abs(maxval - minval)*0.05;
    yrange = [minval-tweakFactor maxval+tweakFactor];
    set(gca, 'YLim', yrange);
    
%     set(gca, 'XLim', [(showStart-200)/sr (showStop+200)]);
    
    xlabel('Seconds');
end
set(findobj(gcf, 'Type', 'axes'), 'Visible', 'off');

% % title
% k = strfind(data.data.map.map1.directory, '\');
% % titleStr = [handles.data.map.mapActive.directory(k(end-2)+1:k(end-1)-1) ', ' ...
% %     handles.data.map.mapActive.directory(k(end)+1:end)];
% titleStr = [data.data.analysis.experimentName ', ' ...
%     data.data.map.map1.directory(k(end)+1:end)];
% text('String', titleStr, 'Units', 'Normalized', 'Position', [0 1.005], ...
%     'FontSize', 12, 'FontWeight', 'Bold', 'Parent', hsub(1), ...
%     'Tag', 'singleTraceMap', 'Interpreter', 'none');

% scalebar lines
if mode ==2
Y = mean(arrayinftozero(:,end))+yFactor/4;
else
Y = mean(array(:,end))+yFactor/4;
end
% Y = min(get(gca, 'YLim'));
hscalebar = line([.1 .2], [Y Y]);
set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');
hscalebar = line([.1 .1], [Y Y+yFactor/2]);
set(hscalebar, 'Color', 'k', 'Tag', 'scaleBarLines');

% scalebar text
ht(1) = text(.12, Y+yFactor/6, '100 ms'); 
ht(2) = text(.12, Y+yFactor/3, [num2str(yFactor/2) ' ' scaleBarTxt]); 
set(ht, 'Color', 'k', 'FontSize', 8, 'Tag', 'scaleBarText');

function plot_small_map_trace(hObject, eventdata, handles)
global OnsetPoint
global DirectResponseAnalysisWindow
userdata=get(gco,'UserData');
figure
plot(userdata{1})
YLIM=get(gca,'YLim');
set(gcf,'name',['Trace ' num2str(userdata{2})]);
line([OnsetPoint OnsetPoint+DirectResponseAnalysisWindow OnsetPoint+DirectResponseAnalysisWindow OnsetPoint OnsetPoint],[YLIM(1),YLIM(1),YLIM(2),YLIM(2),YLIM(1)],'Color','r');