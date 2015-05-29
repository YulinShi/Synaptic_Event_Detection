function peak_scores=score_peaks(Index, Data, TH, fitting,DataRange,TE)
Deriv = diff(Data);
Deriv2= diff(Deriv);
%  Index(end)=[];
%  Index(end)=[];
% Index(1)=[];
%  Index(1)=[];
% mark_Amp=(abs(Data(Index))-min(abs(Data(Index))))/(max(abs(Data(Index)))-min(abs(Data(Index))))*100;
mark_Amp=abs(Data(Index))/max(abs(Data(Index)))*100;

abs_Deriv2=abs(Deriv2);
sum_Deriv2=zeros(1,length(Index));
% size(sum_Deriv2)
% size(abs_Deriv2(Index-1))
% 
% for i = -5:5
%     sum_Deriv2=sum_Deriv2+abs_Deriv2(Index-1-i)
% end
% mark_Deriv2=sum_Deriv2/max(sum_Deriv2)*100;
mark_Deriv2=abs_Deriv2(Index-1)/max(abs_Deriv2(Index-1))*100;


Dist=1./abs(Index-TE);
Dist=Dist/max(Dist)*200;
% sum(abs_Deriv2((Index-2):(Index)))
c_score=zeros(size(Index));
step=100/size(fitting,2);
for i = 1:size(fitting,2)
   tmpfitting =fitting{i}(DataRange)';
%    tmpfitting=tmpfitting(DataRange)';
   c_score=c_score+(tmpfitting(Index)>TH)*step;
end

peak_scores=c_score+mark_Deriv2+Dist+mark_Amp;
% Index
% c_score
% mark_Amp
% mark_Deriv2
%  peak_scores(find(Index<DataRange(1)|Index>DataRange))=0;