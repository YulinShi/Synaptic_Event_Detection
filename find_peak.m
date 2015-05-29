function f = find_peak(Data,PeakFlag,EmptyFlag, TH,fitting,DataRange,Onset,DirectResponseAnalysisWindow,TE)

%this function finds the global peak in a data of certain length T
%Data - 1 x T time series

if size(Data,1) > 1
    Data = Data';
end

%first-order difference
Deriv = diff(Data(DataRange));
Ind = (Deriv>0);
IndP = find(diff(Ind) == -1);

IndP = IndP + 1;

f=[];
if ~isempty(IndP)
    %sort peaks (largest first)
    peak_scores=score_peaks(IndP,Data(DataRange), TH,fitting,DataRange,TE-DataRange(1)+1);
    m=1;

    if isempty(peak_scores)
        f=[];
        return 
    end
    [Peak IndPS] = sort(peak_scores);
    IndPS = fliplr(IndPS);
    f = IndP(IndPS);
        m=1;
    f=f+DataRange(1)-1;

    if PeakFlag == 1
          i=1;
        a=length(f);
        tmp_f=[];
        while i<=length(f)
            if f(i)>=Onset &&  f(i)<=Onset+DirectResponseAnalysisWindow
                            i=i+1;
               continue
            else
              tmp_f=[tmp_f,f(i)];
              i=i+1;
            end
        end
        f=tmp_f;

    else
        f = f;
    end
else
    if EmptyFlag
        [Peak IndPS] = sort(Data);
        IndPS = fliplr(IndPS);
        f = IndPS;
    else
        f=[];
    end
end

