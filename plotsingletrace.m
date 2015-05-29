function [] = plotsingletrace


h_obj = get(gcf, 'CurrentObject');

ud_obj = get(h_obj, 'UserData');

ind=findstr(ud_obj, '_');
cnt = ud_obj(ind(end)+1:end);

x = get(h_obj, 'Xdata');
% x = x -min(x);
y = get(h_obj, 'Ydata');

figure; 
plot(x, y);
title(['trace ', cnt]);
xlim([min(x) max(x)]);

%%% need to plot stim onset and synaptic windown
% stimOn = handles.data.analysis.stimOn;
% line([stimOn stimOn], [min(get(gca,'YLim')) max(get(gca,'YLim'))], 'LineStyle', '-'); 
% edEnd = stimOn + handles.data.analysis.directWindow;
% line([edEnd edEnd], [min(get(gca,'YLim')) max(get(gca,'YLim'))], 'LineStyle', '-');
% synEnd = stimOn + handles.data.analysis.synapticWindow;
% line([synEnd synEnd], [min(get(gca,'YLim')) max(get(gca,'YLim'))], 'LineStyle', '-');

% % Get plot handles and compute trace ranges
% ud = get(gco, 'UserData');
% curcol = get(gco, 'Color');
% if curcol == ecol, curcol = mcol; end
% cnt = ud(3);  % count is always third entry
% pdh    = findobj('Tag', ['Current Trace Duration', num2str(gcf)]);
% stimh  = findobj('UserData', [gcf, stimind, cnt]);
% curh   = findobj('UserData', [gcf, curind, cnt]);
% eventh = findobj('UserData', [gcf, eventind, cnt]);
% 
% plotdur = get(pdh, 'UserData');  plotdur = plotdur(1:end);
% plotsamp  = plotdur/1000*srate;  % samples to plot
% curx = get(curh, 'XData');
% cury = get(curh, 'YData');
% %cury = medsm(cury, 5);
% newx = linspace(0, plotdur, plotsamp);
% 
% gdh = findobj('Tag', ['Grid Size', num2str(gcf)]);
% gridszmat = get(gdh, 'UserData');
% axht = 1/(gridszmat(2)+1);
% newy = maxIdisp*(cury - cury(1)) / axht; % pA
% 
% % Plot indicator box
% if hilitmode == 1 % plot single trace box
%     ncx = min(curx); xcx = max(curx);
%     ncy  = median(cury) - axht/2;
%     xcy  = median(cury) + axht/2;
%     boxx = [ncx xcx xcx ncx ncx];
%     boxy = [ncy ncy xcy xcy ncy];
%     
%     boxh = findobj('Tag', ['TraceBox', num2str(gcf)]);
%     if isempty(boxh)  % no box indicating trace
%         line(boxx, boxy, 'Color', lcol, ...
%             'Tag', ['TraceBox', num2str(gcf)]);
%     else
%         set(boxh, 'XData', boxx, 'YData', boxy);
%     end % if isempty(boxh)
%     
% elseif hilitmode == 2  % highlight selected trace
%     set(get(gca, 'Children'), 'Selected', 'off');
%     set(curh, 'Selected', 'on');
%     
% elseif hilitmode == 3  % change color of selected trace
%     allcurh = findobj('Tag', 'Current Traces');
%     set(allcurh, 'Color', mcol, 'LineWidth', lwid);
%     set(curh, 'Color', ncol, 'LineWidth', 4*lwid);
% end % if hilitmode
% 
% % Set up single trace plot figure
% if stickyfig == 1 & ...
%         ~isempty(intersect(get(0, 'Children'), fignum(singleptr)))
%     pfigpos{singleptr} = get(fignum(singleptr), 'Position');
% end  % if ~isempty...
% singlefig = figure(fignum(singleptr)); clf
% set(singlefig, 'NumberTitle', 'off', ...
%     'Name', 'Single Trace', ...
%     'Position', [pfigpos{singleptr}], 'Color', [1 1 1], ...
%     'MenuBar', 'none');
% 
% % Adjust plot position and y axis scale
% set(gca, 'Position', [0.16 0.13 0.75 0.8]);
% marginperc = 0.1;
% margin = (max(newy) - min(newy)) * marginperc;
% ylim = [min(newy)-margin max(newy)+margin];
% if margin > 0, set(gca, 'YLim', ylim); end
% 
% % Plot stimulus marker
% if     plotstim == 1 % plot stimulus marker line
%     stimy = diff(ylim)*0.1 + ylim(1)*[1,1];
%     stimx(1) = stimdel;
%     stimx(2) = stimdel + stimdur;
%     line(stimx, stimy, 'Color', lcol);
% elseif plotstim == 2 % plot stimulus marker box
%     stimy = [ylim(1) ylim(1) ylim(2) ylim(2)];
%     stimx = stimdel + [0 stimdur stimdur 0];
%     patch(stimx, stimy, lscol, 'EdgeColor', lscol, ...
%         'UserData', [randnum, stimind, cnt]);
% end % if plotstim
% 
% % Plot single trace
% lh = line(newx, newy, 'Color', curcol, 'LineWidth', 2*lwid);
% set(gca, 'TickDir', 'out', 'FontName', labfont, ...
%     'FontSize', asiz, ...
%     'XLim', [0 plotdur], 'ButtonDownFcn', 'grid');
% xlabel('Time (ms)'); ylabel('Current (pA)');
% 
% % Write text info
% %ty(1) = max(cury);
% %ty(2) = 7*(max(cury) - min(cury))/8 + min(cury);
% ty(1) = max(get(gca, 'YLim'));
% ty(2) = min(get(gca, 'YLim'));
% ty(3) = 15*(ty(1)-ty(2))/16 + ty(2);
% th(1) = text(5*plotdur/8, ty(1), ...
%     ['Count = ', num2str(cnt)], ...
%     'Color', lcol, 'FontName', labfont, ...
%     'FontSize', asiz);
% th(2) = text(5*plotdur/8, ty(3), ...
%     ['Time = ', num2str((cnt-1)*minidur)], ...
%     'Color', lcol, 'FontName', labfont, ...
%     'FontSize', asiz);
% 
% % Plot MiniAnal events
% if eventmode == 1 % plot synaptic events
%     eventy = get(eventh, 'YData');
%     eventy = maxIdisp*(eventy - cury(1)) / axht; % pA
%     
%     eventx = get(eventh, 'XData');
%     oldcurx = get(curh, 'XData');
%     eventx = plotdur * (eventx - min(oldcurx)) / ...
%         (max(oldcurx) - min(oldcurx));
%     figure(fignum(singleptr));
%     line(eventx, eventy, ...
%         'Color', ecol, 'LineStyle', 'none', ...
%         'Marker', emrk, 'MarkerSize', 2*mrksiz, ...
%         'MarkerEdgeColor', ecol, 'MarkerFaceColor', ecol);
% end % if eventmode

return;