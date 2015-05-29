f1=figure
f2=figure
for i = 1: length(traceNumber)
    figure
    hold off
    plot(data.data.CorrectedRawData(1100:2600,traceNumber(i)));
    hold on
    plot(M{i},);
end