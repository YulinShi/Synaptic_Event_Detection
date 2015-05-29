function fArray = filterTraceArray_peak_detection(dataArray, sampleRate, filterType, value)
% filterTraceArray
%
% takes in a set of raw traces (e.g. 7000 x 256 array)
% returns a filtered array
% Private version for mapAnalysis2p0
%
% Editing:
% misc changes: gs may 2004
% gs april 2005 -- privatized for mapAnalysis2p0
% --------------------------------------------------------------------------------------
[rows, numTraces, sets] = size(dataArray); % must be 2D data array; col num = trace num
% value = 11;
% filterType = 'mean';
if sets > 1
    dataArray = reshape(dataArray, rows, []); 
end
offset = mean(dataArray(round(0.045*sampleRate)+1:round(0.055*sampleRate)+1, :, :)); % first 10 msec used
offset = repmat(offset, rows, 1);
dataArray = dataArray - offset;
fArray = colfilt(dataArray, [value 1], 'sliding', filterType);
fArray = fArray + offset; % NB: the dataArray is updated with the filtered array
if sets > 1
    fArray = reshape(fArray, rows, numTraces, sets);
end
