% function [spikes,C,outIndex]
function varaout= Align_Peaks(Signal,DataRange, TE, Leads, Trailings, SampleRate,TH_UP,traceno,TH, fitting,diff2,OnsetPoint,DirectResponseAnalysisWindow,DirectResponseAnalysisWindow_2nd,amp_th,b_spon,Spon_Amp_Th)

max_Leads=round((min(max(Leads),mean(Leads)+std(Leads)))*SampleRate);
max_Trailings=round(min((max(Trailings)),mean(Trailings)+std(Trailings))*SampleRate);
min_Leads=round((min(Leads))*SampleRate);
min_Trailings=round((min(Trailings))*SampleRate);
Leads=round(max([(mean(Leads)-std(Leads)),min(Leads)])*SampleRate);
Trailings=round(max([(mean(Trailings)-std(Trailings)),min(Trailings)])*SampleRate);

orig_Leads=Leads;
orig_Trailings=Trailings;
Signal_Length=length(Signal);
Spike_Length=Leads+Trailings+1;
mean_diff2=mean(diff2);
std_diff2=std(diff2);
ind1= find(TE>Leads);
ind2= find(TE<(Signal_Length-Trailings));
[C,IA,IB]=intersect(ind1,ind2);
spikes=inf(Spike_Length, length(TE));
varaout={};
outIndex=[];
startIndex=[];
endIndex=[];
Plot_endIndex=[];
shift_ms=[];
count_valid_spikes=0;
%Align Spikes
m=1;
prev_end=1;
fail_reason=cell(length(TE),1);
for i =1:length(TE)
    
    Leads=orig_Leads;
    Trailings=orig_Trailings;
    count_valid_spikes=count_valid_spikes+1;
    if TE(i)-Leads>DataRange(end)
        break;
    end
    if TE(i)<prev_end
        fail_reason{i}='0';
        continue
    end
    if (TE(i)-Leads<prev_end)
        Leads=TE(i)-prev_end;
    end
    
    if (TE(i)+Trailings)>Signal_Length
        Trailings=Signal_Length-TE(i);
    end
    sig=Signal((TE(i)-Leads):end);
    total_Index=find_peak(-Signal,1,0,TH,fitting,((TE(i)-Leads)) : (TE(i)+Trailings),OnsetPoint,DirectResponseAnalysisWindow,TE(i));
    
    if isempty(total_Index)
        fail_reason{i}='Bad Score';
        continue;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Peak_exist=0;
    for Index_Length=1:length(total_Index)
        
        
        Index=total_Index(Index_Length);
        
        if Index>DataRange(end)
            fail_reason{i}='Out DataRange';
            continue;
        end
        
        Leads=max_Leads;
        Trailings=max_Trailings;
        if Index-Leads<=prev_end
            sig=Signal(prev_end : Index+Trailings);
            tmpstart=prev_end;
            Leads=Index-prev_end;
            tmpend=Index+Trailings;
        elseif (Index+Trailings>Signal_Length)
            tmpstart=Index-Leads;
            sig=Signal(Index-Leads : end);
            tmpend=Signal_Length;
            Trailings=Signal_Length-Index;
        else
            sig=Signal(Index-Leads : Index+Trailings);
            tmpstart=Index-Leads;
            tmpend=Index+Trailings;
        end
        
        
        if exist('b_spon','var')
            Beginning_index=find_peak_for_edge(Signal,1,0,TH,fitting,Index,((Index-Leads)) : (Index+Trailings),OnsetPoint,DirectResponseAnalysisWindow,prev_end,amp_th,b_spon,Spon_Amp_Th);
        else
            Beginning_index=find_peak_for_edge(Signal,1,0,TH,fitting,Index,((Index-Leads)) : (Index+Trailings),OnsetPoint,DirectResponseAnalysisWindow,prev_end,amp_th);
        end
        firstPeak=Beginning_index(1)-(Index-Leads)+1;
        SecondPeak=Beginning_index(2)-Index+1;
        Plot_SecondPeak=Beginning_index(4)-Index+1;
        v_start=Signal(Beginning_index(1));
        v_end=Signal(Beginning_index(2));
        if ~exist('b_spon','var')
            if Beginning_index(3)<=OnsetPoint+DirectResponseAnalysisWindow
                fail_reason{i}='Direct';
                continue
            end
        end
        
        if isempty (firstPeak)
            fail_reason{i}='amp fail 1';
            continue
        end
        if isempty(SecondPeak)
            fail_reason{i}='amp fail 2';
            continue
        end
        
        v_start=v_start-Signal(Index);
        v_end=v_end-Signal(Index);
        if v_start<=0 | v_end<=0
            fail_reason{i}='amp fail 3';
            continue
        end
        
        Plot_sig=sig(firstPeak:Plot_SecondPeak);
        sig=sig(firstPeak:Leads+SecondPeak);
        
        diff2_sig=diff(sig,1);
        sum_diff2_sig=sum((diff2_sig));
  
        if exist('b_spon','var')
            if Spon_Amp_Th >abs(v_start)
            fail_reason{i}='Amp_Th';
            continue
        end
        else
        if amp_th >abs(v_start)
            fail_reason{i}='Amp_Th';
            continue
        end
        end
        over_RH=0;
        for fitting_i=1:length(fitting)
            if (Index>=OnsetPoint)&&(fitting{fitting_i}(Index)>=TH_UP) && (Index<=(OnsetPoint+DirectResponseAnalysisWindow_2nd))
                over_RH=1;
                break
            end
        end
        if over_RH
            fail_reason{i}='5';
            continue
        end
        Peak_exist=1;
        
        break
    end
    if Peak_exist==0
        fail_reason{i}='4';
        continue
    end
    
    shift_ms=[shift_ms,(TE(i)-Index)/(SampleRate)];
    outIndex=[outIndex,Index];
    startIndex=[startIndex,tmpstart+firstPeak-1];
    endIndex=[endIndex,tmpstart+Leads+SecondPeak-1];
    Plot_endIndex=[Plot_endIndex,tmpstart+Leads+Plot_SecondPeak-1];
    prev_end=endIndex(end);
    spikes(1:(endIndex(end)-startIndex(end)+1),i)=sig(1:end);
    m=m+1;
end
p=1;
for spikes_i=1:size(spikes,2)
    if ~isinf(spikes(1,spikes_i))
        spikes(:,p)=spikes(:,spikes_i);
        p=p+1;
    end
end


varaout{1}=startIndex;
varaout{2}=endIndex;
varaout{3}=outIndex;
amplitude=[];
sumamplitude=[];

for i = 1 :length(startIndex)
    sumamplitude=[sumamplitude,2*abs(sum((Signal(startIndex(i):outIndex(i))-Signal(startIndex(i)))))];
    amplitude=[amplitude,-(Signal(outIndex(i))-Signal(startIndex(i)))];
end
varaout{4}=sumamplitude;
varaout{5}=amplitude;
varaout{6}=fail_reason;
