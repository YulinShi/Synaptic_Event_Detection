function f = find_peak_for_edge(Data,PeakFlag,EmptyFlag, TH,fitting,Index,DataRange,Onset,DirectResponseAnalysisWindow,prev_end,amp_th,b_spon,Spon_Amp_Th)

%this function finds the global peak in a data of certain length T
%Data - 1 x T time series

if exist('b_spon','var')
amp_th=Spon_Amp_Th;
end

if size(Data,1) > 1
    Data = Data';
end

%first-order difference
Deriv = diff(Data);
Ind = (Deriv>0);
IndP = find(diff(Ind) == -1);

IndP = IndP + 1;

% 
% figure
% plot(Deriv);
% hold on
% plot(Deriv2,'g')

f=[DataRange(1),DataRange(end),DataRange(1),DataRange(end)];
if ~isempty(IndP)
    %sort peaks (largest first)
    peak_scores=IndP;
       if isempty(peak_scores)
f=[DataRange(1),DataRange(end),DataRange(1),DataRange(end)];
return
       end
       Direct_Lead_IndP=DataRange(1);
    Lead_IndP=DataRange(1);
    if peak_scores(1)>Index
        f=[DataRange(1),DataRange(end),DataRange(1),DataRange(end)];
    else
        Lead_IndP_tmp=IndP(find((IndP>=DataRange(1)) & (IndP<Index)));
        if ~isempty(Lead_IndP_tmp)
        
            [Lead_IndP,Lead_IndP_i]=max(Data(Lead_IndP_tmp));
            while 1
            Lead_IndP=Lead_IndP_tmp(Lead_IndP_i);
           break
            end
        else
            tmp=find(peak_scores<Index);
            if ~isempty(tmp)
                Direct_Lead_IndP=(peak_scores(tmp(end)));
            end
        end
    end
    end_IndP=DataRange(end);
    plot_end_IndP=DataRange(end);
    if peak_scores(end)>=Index

        end_IndP_tmp=IndP(find((IndP>Index) & (IndP<=DataRange(end))));
        if ~isempty(end_IndP_tmp)
        
            [end_IndP,end_IndP_i]=max(Data(end_IndP_tmp));
            
           
%             if exist('amp_th','var')
             end_IndP=end_IndP_tmp(1);
             end_end_IndP_tmp=[end_IndP_tmp,DataRange(end)];
             for end_i=2:length(end_end_IndP_tmp)
                if Data(end_end_IndP_tmp(end_i))>Data(end_IndP)
                    data_end_end=-Data(end_IndP:end_end_IndP_tmp(end_i))+Data(end_IndP);
                    if sum(data_end_end>amp_th)
                        break
                    else
                        end_IndP=end_end_IndP_tmp(end_i);
                    end
                end
             end
%             else
%               end_IndP=end_IndP_tmp(end_IndP_i);
%             end
        end
        [tmp,plot_end_IndP]=max(Data(Index:end_IndP));
        plot_end_IndP=Index+plot_end_IndP-1;
    end
     f=[Lead_IndP,end_IndP,Direct_Lead_IndP,plot_end_IndP];
%     for i = 1:length(peak_scores)
%         if peak_scores(m)< DataRange(1) | peak_scores(m)>DataRange(end)
%             peak_scores(m)=[];(
%         else
%             m=m+1;
%         end
%     end
 

%     for i = 1:length(f)
%         if f(m)< DataRange(1) | f(m)>DataRange(end)
%             f(m)=[];
%         else
%             m=m+1;
%         end
%     end
 
else
    if EmptyFlag
        [Peak IndPS] = sort(Data);
        IndPS = fliplr(IndPS);
        f = IndPS;
    else
f=[DataRange(1),DataRange(end),DataRange(1),DataRange(end)];
    end
end

