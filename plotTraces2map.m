function[a,b,c,d,e,f,g,h,result1,out_std,result2,result3,threshold_onset,nonrsp_array,nonrsp_mean,nonrsp_std,noresponse_mean]= plotTraces2map(data)

%%%%%% sort the traces into direct and inderect.
%%%% roughly get the result of each trace's response, stored in result.
%%%% accept parameters through interactive dlg.

a=0;
b=0;
c=0;
d=0;
e=[];
f=[];
g=0;
h=0;
result =[];
if data.seq.load ==0
    msgbox('load data first');
    return;
end
% def.stimcnts = '1:332'; %stimulation count series
% def.aligncnts = '323:337'; %alignment count series
def.stimdur  = '10';          % laser onset duration to plot in ms
def.stimdel = '50'; % stimulation delay in ms
def.indirectwindow = '150';
def.polarity = '1'; % polarity: excitatory or inhibitory
def.stdNum ='3';
def.threshold_onset='15';
def = struct2cell (def);    % converts def from structure to cell

%%Input parameters for trace plotting

prompt = { 'Enter laser laser direct response window',...
    'Enter stim delay in ms', 'Enter indirect window period in ms', 'Enter polarity: excitatory 1, inhibitory 2','How many deviations','Enter threshold'};
dlg_title = 'Enter parameters for trace plotting';
num_lines= 1;
parameters = inputdlg(prompt, dlg_title, num_lines, def);

% stimcnts = str2num(char(parameters(1)));      % stimulation count series
% aligncnts = str2num(char(parameters(2)));  % alignment count series
stimdur = str2num(char(parameters(1)));       % plot duration in ms
stimdel = str2num(char(parameters(2)));      % stimulation delay in ms
indirectwindow = str2num(char(parameters(3)));   % indirect window in ms
polarity = str2num(char(parameters(4)));       % soma location
stdNum = str2num(char(parameters(5)));
out_std = stdNum;
threshold_onset=str2num(char(parameters(6)));
%%constants
downsamp = 1;         % factor to downsample (1=none)
srate   = 10000/downsamp; % samples per second, (AD convert: 10 KHz)
% framdur = 400;        % duration of frame / block (ms) // the DOS program records up to 400 ms of data for one stimulation
%
% nframsamp = srate*framdur/1000;  % samples in a frame

basesamplepoints_stimdel = stimdel*srate/1000;
laseronpoints = stimdur*srate/1000;%% 10 ms

base_pts = basesamplepoints_stimdel;
baselaser_pts = basesamplepoints_stimdel + laseronpoints;
indirectwindow_pts = indirectwindow*srate/1000; %% normally 140 ms after laser off(10 ms)

%%%%%%%%5 paramers into global myData
a = polarity;
b = basesamplepoints_stimdel;
c = laseronpoints;
h = indirectwindow_pts;
%
% alpha   = 20;          % output gain selection (usu 20 for photostimulation trials)
% beta    = 1;          % feedback resistor gain selection (1)
% maxV    =  10;        % largest voltage value
% minV    = -maxV;      % smallest voltage value
% maxADC  =  2^11;      % largest quantized value
% minADC  = -maxADC;    % smallest quantized value


% % Plotting and display constants
% xyinc = 84;           % xy increment is 84  (for 50 um)
% xyadj = 1.68;         % adjustment for x and y coordinates
% xymic  = xyinc / xyadj;  % increment in microns

if ~isempty(data.data)

    if~isempty(data.data.analysis.numberOfMaps)
        for i = 1:data.data.analysis.numberOfMaps
            string =['data.data.map.map',num2str(i),'.bsArray;'];

            eval(['PHarray =' string]);
            PHarray4base = PHarray;
            %%%% for plotting sorted response and measurement display in Plotdirecttraces
            %%in single electro and also calculate in dual e.g. sumPH_baseamp_array = sum (PHarray4base(1 : base_pts, :), 1)';


            tempstr = ['plotting in map  map',num2str(i)];
            figure('name',tempstr);
            maxYph = max(max(PHarray));%max(max(PHarray(:, 1 : length (goodtraces)) ));
            minYph = min(min(PHarray));  %%(:, 1 : length (goodtraces)) ));
            %     minY = minYph;
            %     maxY = maxYph;
            % maxYCT = max(max(CTarray)) %(:, 1 : length (goodtraces)) ));
            % minYCT = min(min(CTarray)) %(:, 1 : length (goodtraces)) ));
            %
            % maxY = max(maxYph, maxYCT);
            % minY = min(minYph, minYCT);


            subplotnumb = 1;
            subplot(4, 1, subplotnumb);
            for i = 1 : size(PHarray, 2)  %length (goodtraces)% good tracenumbers size(PHarray, 2)  %%colum number - trace number
                xvect = 1: length(PHarray); %%row number /point number
                line( xvect, PHarray(:, i) ); %, 'UserData', { num2str(i), num2str(basesamplepoints_stimdel)},'Tag', ['Traces', num2str(i)], 'ButtonDownFcn', 'PhotoStimPopTrace'
                hold on;
            end
            global lengthofdata
            lengthofdata=size(PHarray, 1);
            %% plot laser stimulation marker
            % stimpx(1) = basesamplepoints_stimdel;
            % stimpx(2) = laseronpoints + basesamplepoints_stimdel;
            %
            % line([stimpx(1), stimpx(1), stimpx(2), stimpx(2), stimpx(1)], [minY, maxY, maxY, minY, minY], 'color', 'r');
            %
            % title('photostimulation trials');
            % ylabel('amplitude PA');
            %
            % Ylim([minYph maxYph]);
            %
            % subplotnumb = subplotnumb +1;
            %
            % subplot(2, 1, subplotnumb);

            % for i = 1 : size(CTarray, 2)  %%colum number - trace number
            %     xvect = 1: length(CTarray);
            %     h(i) = line( xvect, CTarray(:, i)); %, 'UserData', {num2str(i), num2str(basesamplepoints_stimdel)}, 'Tag', ['Traces', num2str(i)], 'ButtonDownFcn', 'PhotoStimPopTrace'
            %     hold on;
            % end
            %
            % title('control trials');
            % xlabel('time ms');
            %
            % Ylim([minY maxY]);
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sort photostimulation responses%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if polarity == 1  %% excitatory, inward currents
                syn_negThreshold = mean(PHarray(1:basesamplepoints_stimdel, :)) - stdNum* std(PHarray(1:basesamplepoints_stimdel, :));%% basesamplepoints_stimdel: 10 kHz 20 ms

            end

            if polarity == 2 %% inhibitory, outward currents
                syn_negThreshold = mean(PHarray(1:basesamplepoints_stimdel, :)) + stdNum* std(PHarray(1:basesamplepoints_stimdel, :));% basesamplepoints_stimdel: 10 kHz 20 ms
            end
            % syn_negThreshold = min(PHarray(3001:4000, :)) - 3* std(PHarray(3001:4000, :));
            %%%%%%parameters into global myData
            d = syn_negThreshold;
            %%%

            %basesamplepoints_stimdel
            [R,C] = size(PHarray);
            colnum = C;
            threshold = repmat(syn_negThreshold, R, 1);
            %size(threshold)
            if polarity == 1  %% excitatory, inward currents

                TH_PHarray = PHarray<=threshold;  %% TH_PHarray: 0 or 1syn_negThreshold = median(PHarray(1:200, :)) - 2* std(PHarray(1:200, :));
            end

            if polarity == 2  %% inhibitory, outward currents
                TH_PHarray = PHarray>=threshold;
            end


            %% a little trick to adjust some data
            for i = 1 : colnum
                if length(find(TH_PHarray(1:basesamplepoints_stimdel,i)))  %% find non-zero before onset and modify the trace
                    PHarray(1:basesamplepoints_stimdel, i) = median(PHarray(1:basesamplepoints_stimdel, i));
                    TH_PHarray(1:basesamplepoints_stimdel, i) = 0;

                end
            end

            %%sort responses into direct, indirect and null responses

            base_pts = basesamplepoints_stimdel;
            baselaser_pts = basesamplepoints_stimdel + laseronpoints;
            indirectwindow_pts = indirectwindow*srate/1000; %% normally 140 ms after laser off(10 ms)


            % figure('name','unknow plotting');

            minY = min(min(PHarray));
            maxY = max(max(PHarray));


            % subplotnumb = 2;

            % subplot(4, 1, subplotnumb);
            %
            % m1 =0;
            % for i = 1 : colnum % good tracenumbers size(PHarray, 2)  %%colum number - trace number
            %
            %     if length(find (TH_PHarray(1:base_pts,i)))  %% before onset
            %         xvect = 1: length(PHarray); %%row number /point number
            %         line(xvect, PHarray4base(:, i), 'UserData', { num2str(i), num2str(base_pts), num2str(syn_negThreshold(i)), laseronpoints }, 'Tag', ['Traces', num2str(i)], 'ButtonDownFcn', 'PhotoStimPopTrace' );
            %         hold on;
            %         m1 = m1 +1;
            %     end
            % end

            Ylim([minY maxY]);

            %  subplotnumb = 2;
            subplotnumb =subplotnumb +1;
            subplot(4, 1, subplotnumb);
            m2 =0; nn=0;


            %%%%%find the direct_array, the index is in the first of each
            %%%%%trace.
            Direct_array = [];
            for i = 1 : colnum % good tracenumbers size(PHarray, 2)  %%colum number - trace number

                if ~length(find (TH_PHarray(1:base_pts,i)))& length(find ( TH_PHarray( (base_pts+1):baselaser_pts, i)))  %% within onset of 10 ms; consided as direct response
                    xvect = 1: length(PHarray); %%row number /point number
                    line(xvect, PHarray4base(:, i), 'UserData', { num2str(i), num2str(base_pts), num2str(syn_negThreshold(i)), laseronpoints}, 'Tag', ['Traces', num2str(i)], 'ButtonDownFcn', 'PhotoStimPopTrace' );
                    hold on;
                    m2= m2+1;
                    Direct_array(:, m2) = [i; PHarray4base(:, i)];

                    nn = nn + 1;
                    dirtrace(nn) = i;

                end
            end
            %%%%%%%% Direct_array into global myData.data
            if ~isempty(Direct_array)
                e = Direct_array;
            end
            %% plot laser stimulation marker
            stimpx(1) = base_pts;
            stimpx(2) = baselaser_pts;

            line([stimpx(1), stimpx(1), stimpx(2), stimpx(2), stimpx(1)], [minY, maxY, maxY, minY, minY], 'color', 'r');

            Ylim([minY maxY]); text(0.7*max(xvect), 0.5*minY+maxY, ['# of direct ones = ', num2str(m2)]);


            %subplotnumb=3
            m3 =0;
            subplotnumb = subplotnumb + 1;
            subplot(4, 1, subplotnumb);


            %%%find the indirect_array;
            Indirect_array=[];
            for i = 1 : colnum% good tracenumbers size(PHarray, 2)  %%colum number - trace number

                if ~length(find (TH_PHarray(1:baselaser_pts, i))) & length(find ( TH_PHarray( (baselaser_pts +1) :(baselaser_pts+indirectwindow_pts) , i) ))  %% >= 30 ms to 180 ms; consided as indirect response
                    xvect = 1: length(PHarray); %%row number /point number
                    line(xvect, PHarray4base(:, i), 'UserData', { num2str(i),  num2str(base_pts), num2str(syn_negThreshold(i)), laseronpoints}, 'Tag', ['Traces', num2str(i)], 'ButtonDownFcn', 'PhotoStimPopTrace' );
                    hold on;
                    m3 = m3 +1;

                    Indirect_array(:, m3) = [i; PHarray4base(:, i)];
                end
            end
            %%%%%Indirect_array into global myData.
            if ~isempty(Indirect_array)
                f = Indirect_array;
            end

            %%%%% find the nullresponse Array.

            nonrsp_array=[];
            k=0;
            %???yu
            for i=1:colnum
                if isempty(find(Direct_array(1,:)==i))&& isempty(find(Indirect_array(1,:)==i))
                    k=k+1;
                    nonrsp_array(:,k)=[i;PHarray4base(:,i)];
                end
            end


            %% plot laser stimulation marker
            stimpx(1) = base_pts;
            stimpx(2) = baselaser_pts;

            line([stimpx(1), stimpx(1), stimpx(2), stimpx(2), stimpx(1)], [minY, maxY, maxY, minY, minY], 'color', 'r');

            Ylim([minY maxY]); text(0.7*max(xvect), 0.5*minY+maxY, ['# of indirect ones = ', num2str(m3)]);


            m4= 0; mm=0;
            %subplotnumb=4;
            subplotnumb = subplotnumb + 1;
            subplot(4, 1, subplotnumb);
            noresponse_mean=0;
            for i = 1 : colnum% good tracenumbers size(PHarray, 2)  %%colum number - trace number

                if ~length(find (TH_PHarray(1:(baselaser_pts+indirectwindow_pts),i)))  %% >= 30 ms to 170 ms; consided as indirect response
                    xvect = 1: length(PHarray); %%row number /point number
                    line(xvect, PHarray4base(:, i), 'UserData', { num2str(i),  num2str(base_pts), num2str(syn_negThreshold(i)), laseronpoints}, 'Tag', ['Traces', num2str(i)], 'ButtonDownFcn', 'PhotoStimPopTrace' );
                    hold on;
                    noresponse_mean =noresponse_mean+ mean(PHarray(baselaser_pts+1:baselaser_pts+indirectwindow_pts, i));
                    m4 = m4+1;
                    mm = mm+1;
                    nulltrace(mm) = i;
                end
            end
            noresponse_mean = noresponse_mean/m4;
            %% plot laser stimulation marker
            stimpx(1) = base_pts;
            stimpx(2) = baselaser_pts;

            line([stimpx(1), stimpx(1), stimpx(2), stimpx(2), stimpx(1)], [minY, maxY, maxY, minY, minY], 'color', 'r');

            Ylim([minY maxY]); text(0.7*max(xvect), 0.5*minY+maxY, ['# of null ones = ', num2str(m4)]);



            sumPHAmp_array = sum (PHarray( (baselaser_pts+1) : (baselaser_pts+indirectwindow_pts), :), 1)';

            if polarity == 1  %% if excitatory
                sumPHarray_corrected = sumPHAmp_array;  %% for manual correction
            end

            if polarity == 2

                sumPHarray_corrected2 = sumPHAmp_array;
            end
            result1 = sumPHAmp_array/indirectwindow_pts;
            %


            %%%%%%%%%%%%%%%%%%%%%    Yulin Added       %%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%  Calculate Onset Latency   %%%%%%%%%%%%%%%%%%%%%
            if polarity == 1
                output= arrayEventFinder(PHarray( (baselaser_pts+1) : (baselaser_pts+indirectwindow_pts), :), -threshold_onset, 'negative'); %%%Get the onset points
            else
                output= arrayEventFinder(PHarray( (baselaser_pts+1) : (baselaser_pts+indirectwindow_pts), :), threshold_onset, 'positive'); %%%Get the onset points
            end
            result3=inf(size(PHarray,2),1);
            result3(output(:,1))=output(:,2)/srate*1000;%%%Calculate Onset Latency based on onset points

            %%%%%%%%%%%%%%%%%%%%%%%%  Calculate Min value   %%%%%%%%%%%%%%%
            [result2,tempresult3] = min(PHarray( (baselaser_pts+1) : (baselaser_pts+indirectwindow_pts), :));%%result2,min value.
            result2 =result2';
            %%%%%%%%%%%%%%%%%%%%%%%  End Yulin Added       %%%%%%%%%%%%%%%
            %
            %%%%%%%%%%%%%%%%%%%    Yulin Deleted       %%%%%%%%%%%%%%%%%%

            %        [result2,result3] = min(PHarray( (baselaser_pts+1) : (baselaser_pts+indirectwindow_pts), :))%%result2,min value.
            %         result3 = result3/srate*1000;%%% result3,onset latency;
            %         result2 =result2'
            %         result3 = result3';
            %%%%%%%%%%%%%%%%%%%%%%%%%%   End Yulin Deleted   %%%%%%%%%%%%%%%%%%%%


            %5 extract mean and std from nonrsp_array;

            nonresp_mean=[];
            nonrsp_std=0;
            nonrsp_mean=0;
            for i = 1: size(nonrsp_array,2)
                tmp = nonrsp_array((baselaser_pts+2):(baselaser_pts+indirectwindow_pts+1),i);
                nonresp_mean =[nonresp_mean,mean(tmp(:))];
            end
            nonrsp_mean = mean(nonresp_mean(:));
            nonrsp_std =std(nonresp_mean);


        end
        g=1;

    else
        msgbox('wrong map number');
        return;
    end
else
    msgbox('no data input');
    return;
end

% sumCTAmp_array = sum (CTarray( (baselaser_pts+1) : (baselaser_pts+indirectwindow_pts), :), 1)';
% sumAmp_array = sumPHAmp_array; %%%  - sumCTAmp_array;  %% summed photostimulation traces - summed control traces (30 - -170 ms)


%
% %% Save corrected
% button = questdlg('Do you want to save input paras and summed PHarray and CTarray?','Save?', 'Yes','No','Yes');
%
% if strcmp(button,'Yes')
%
%     para.stimcnts = stimcnts;
%     para.aligncnts = aligncnts;
%     para.stimdur = stimdur;
%     para.stimdel = stimdel;
%     para.indirectwindow = indirectwindow;
%     para.polarity = polarity;
%
%     constant.downsamp = downsamp;
%     constant.srate  = srate;
%     constant.framdur = framdur;
%     constant.nframsamp = nframsamp;
%     constant.basesamplepoints_stimdel = basesamplepoints_stimdel;
%     constant.laseronpoints = laseronpoints;
%     constant.base_pts = base_pts;
%     constant.baselaser_pts = baselaser_pts;
%     constant.indirectwindow_pts = indirectwindow_pts;
%
%     cd(pathname);
%     savename = [filename, '_', num2str(polarity), 'LoadStimtraces.mat']
%     [savefilename, savepathname] = uiputfile(savename, 'Save as');
%     savefile = [savepathname, savefilename]
%     save (savefile, 'para', 'constant', 'sumPHAmp_array', 'sumCTAmp_array');
% else
%     msgbox('No file save');
%
% end
%