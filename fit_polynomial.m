function [Model,c,SNR,range,diff2, Leads, Trailings,orig_data,Phase_Diff] = fit_polynomial(DataRange,Data,Order,Demean,PlotFlag, SFr)

%This function looks at the time series given by Data and
%it fits a polynomial of the form:
%   y = p(1)*t^Order + p(2)*t^(Order-1) + ... + p(Order)*t + p(Order+1)
%Order is the order of the polynomial fit
%Demean - string 'demean' subtract the mean of the filter
%
%PlotFlag = 1 --> generate plot, otherwise don't

%sampling frequency
% SFr = 10000; %[Hz]

figure
plot(DataRange,Data)

% set(gca,'Xlim',[1 length(Data)]);
set(1,'Position',[50 50 800 600])
try
[t1, t2] = ginput(2);
catch
end
    if ~exist('t1','var') | ~exist('t2','var')
        return
    end
t1 = round(t1)-DataRange(1)+1;
t2 = round(t2)-DataRange(1)+1;
orig_data=Data(t1(1):t1(2));
dData = Data(t1(1):t1(2)) - max(Data(t1(1):t1(2)));
[minvalue,index]=min(dData);
Leads=index*1000/SFr;
Trailings=(length(dData)-index)*1000/SFr;
BG = setdiff(Data,Data(t1(1):t1(2)));
snr(1) = max(abs(dData))/mad(BG,1);
diff2_dData=diff(dData,1);
Int_Diff2_dData=sum(diff2_dData);
diff2=Int_Diff2_dData;
t = [0:length(dData)-1]*1000/SFr;    %[milliseconds]
Nt = length(t);

T = [];
for i = 0:Order
    T = [T; t.^i];
end

W = sqrt(abs(dData))';
W = 1 + W/max(W);  %gives more weight to matching larger values of data

T = [W*ones(1,Order+1)].*fliplr(T');
Y = W.*dData';

%P = inv(T'*T)*T'*Y;

%pinv the same as above, but it doesn't complain about singularity
P = pinv(T) * Y;
%P = polyfit(t,dData,Order);

if PlotFlag == 1
    figure
    %t = [0:length(dData)-1]*1000/SFr;    %[milliseconds]
    subplot(211)
    bh = plot(t,dData,'-o','LineWidth',1,'MarkerSize',4,'MarkerFaceColor','b');
    hold on
    Model = zeros(1,Nt);
    for j = 1:length(P)
        Model = Model + P(j)*t.^(Order+1-j);
    end
    rh = plot(t,Model,'r','LineWidth',2);
    xlabel('time [ms]')
    set(gca,'Xlim',[t(1) t(end)])
    legend([bh, rh],'data','model');
end

%normalize the filter


tmp_Model = Model - mean(Model);
Filter = tmp_Model/sum(abs(tmp_Model));

c = conv(Filter,Data);
in = (length(Filter)/2);
ct = c(round(in):round(length(c)-in));
% length(ct)
% length(DataRange)
BG = setdiff(ct,ct(t1(1):t1(2)));
snr(2) = max(abs(ct))/mad(BG,1);
disp(['SNR gain factor: ' num2str(snr(2)/snr(1))])
SNR=snr(2)/snr(1);
if PlotFlag == 1
    subplot(212)
    %plot(ct/max(abs(ct)),'g')
    gh = plot(DataRange,ct,'g','LineWidth',1);
    hold on
    %     plot(Data/max(abs(Data)),'b')
    bh = plot(DataRange,Data,'b-','LineWidth',1);
    rh = plot((t1(1):t1(2))+DataRange(1)-1,Data(t1(1):t1(2)),'r','LineWidth',2);
    plot([t1(1) t1(1)]+DataRange(1)-1,[min(Data(t1(1):t1(2))), max(ct(t1(1):t1(2)))],'k:'); 
    plot([t1(2) t1(2)]+DataRange(1)-1,[min(Data(t1(1):t1(2))), max(ct(t1(1):t1(2)))],'k:');
    [MAX_C,IC]=max(ct(t1(1):t1(2)));
    [MAX_O,IO]=min(Data(t1(1):t1(2)));

    Phase_Diff=IC-IO;

%     set(gca,'Xlim',[DataRange(1) DataRange(end)]);
    legend([gh,bh,rh],'convolution','data','fitted data')
end

[xc, lags] = xcov(Data,ct);

[tmp, lag_ind] = max(abs(xc));
lag_m = lags(lag_ind);
disp(['Phase difference: ' num2str(lag_m) ' samples'])
mc = max(ct);
disp(['Maximum coefficient: ' num2str(mc)])
range=t1(1) :t1(2);
