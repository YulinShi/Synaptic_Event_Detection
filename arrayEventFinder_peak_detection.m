function output = arrayEventFinder_peak_detection(array, threshold, polarity)
% arrayEventFinder
%
% array -- set of traces
% threshold -- 
% polarity -- 'positive' or 'negative'
%
% Notes ... cols, rows, indices, values 
%
% Editing:
% gs april/may 2005 -- modified for mapAnalysis
% ----------------------------------------------------
[R,C] = size(array);
[r,c] = size(threshold);

if r*c > 1
    if r > c
        threshold = threshold';
    end
    threshold = repmat(threshold, R, 1);
    array = array-threshold;
    threshold = 0;
end

switch polarity
case 'positive' %| 'up'
    [r,c] = find(array >= threshold);
case 'negative' %| 'down'
    [r,c] = find(array <= threshold);
end

d = [c r];
dflip = flipdim(d,1);

if isempty(d)
    output = [];
    return
end

[B,I,J] = unique(dflip(:,1)); % B is a list of all the traces with events

output(:,1) = B;
temp = dflip(:,2);
output(:,2) = temp(I);
output(:,3) = (output(:,1)-1)*size(array,1)+output(:,2);
output(:,4) = array(output(:,3));
events = zeros(C,1);
events(output(:,1)) = output(:,2);

% now the new way, with duration threshold
% for synaptic input maps, the threshold is set in the uicontrol as synDuration
