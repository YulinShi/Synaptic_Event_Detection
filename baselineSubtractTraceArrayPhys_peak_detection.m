function bsArray = baselineSubtractTraceArrayPhys_peak_detection(dataArray,  baselineStart, baselineEnd)
% baselineSubtractTraceArrayPhys
%
% takes in a set of raw traces (e.g. 7000 x 256 array)
% returns a baseline-subtracted array
% Private version for mapAnalysis2p0
%
% Editing:
% misc changes gs may 2004
% gs april 2005 -- privatized for mapAnalysis2p0
% --------------------------------------------------------------------------
[rows, numTraces, sets] = size(dataArray); % must be 2D data array; col num = trace num
if sets > 1
    dataArray = reshape(dataArray, rows, []); 
end
baselineStart = round(baselineStart)+1;
baselineEnd = round(baselineEnd)+1;
baselineMeans = mean(dataArray(baselineStart:baselineEnd,:));
% baselineSDs = std(dataArray(baselineStartIndex:baselineEndIndex,:));
baselineArray = repmat(baselineMeans, rows, 1);
bsArray = dataArray - baselineArray; 
if sets > 1
    bsArray = reshape(bsArray, rows, numTraces, sets);
end
