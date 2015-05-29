
function plotdirecttrace()
global K
global Tc
global To
global traceNumber
global Model
global Convolution
global Original_Trace
global cellcount
global exp_fit
global diff2
global Coefficient

global SampleRate
global type_str
global DataRange
global trailings
global leads
global StatusBoxHandle
global TotalTrained
global orig_training_data
global training_order
global Phase_Diff
l=get(gco,'UserData');

% 
% if exp_fit==1
% 
%     type_str='exp';
%     
%     [m,xopt]=fit_double_exponential(Original_Trace(DataRange,l)','inh',1, SampleRate);
% 
%     % tempC= tempC(ceil(length(Filter)/2):2101+floor(length(Filter)/2));
%     %
%     % subplot(2,1,2)
%     % plot(Original_Trace(900:3000,l)/max(abs(Original_Trace(900:3000,l))));
%     % hold on
%     % plot(tempC/max(abs(tempC)),'k');
%     tmp_type='exponential';
%     % button = questdlg(['Do you want to save exponential Model and Convolution Value?'],'Save?', 'Yes','No','Yes');
% 
% else
    type_str='pol';
    try
    [m,Conv, SNR,range,diff2_tmp,leads_tmp,trailings_tmp,tmp_orig_training_data,tmp_Phase_Diff] = fit_polynomial(DataRange,Original_Trace(DataRange,l)',training_order,'demean',1,SampleRate);
    catch
        return
    end

    tmp_type='polynomial';
    % test_high_pass(Original_Trace(1000:1500,l),m,NUM,DEN);
% end




m_dmean=m-mean(m);
Filter = m_dmean/sum(abs(m_dmean));
tempC=conv(Original_Trace(DataRange,l),Filter);
C_M_Diff =round(length(Filter)/2);
CT = tempC(C_M_Diff:length(tempC)-C_M_Diff);

MAX_CT=max(tempC(range));
button = questdlg(['Do you want to save ',tmp_type,' Model and Convolution Value?'],'Save?', 'Yes','No','Yes');
if strcmp(button,'Yes')
    if exp_fit==1
        K=[K,xopt(1)];
        Tc=[Tc,xopt(2)];
        To=[To,xopt(3)];
    end
    diff2=[diff2,diff2_tmp];
    traceNumber=[traceNumber,l];
    Model{cellcount}=m;
    Convolution{cellcount}=CT;
    orig_training_data{cellcount}=tmp_orig_training_data;
    Coefficient=[Coefficient,MAX_CT];
    cellcount=cellcount+1;
    trailings=[trailings,trailings_tmp];
    leads=[leads,leads_tmp];
    Phase_Diff=[Phase_Diff,tmp_Phase_Diff];
    TotalTrained=TotalTrained+1;
    set(StatusBoxHandle,'string',['Total Trained:' num2str(TotalTrained)]);
end
