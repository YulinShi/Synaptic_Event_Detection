
function [fcn,start_stop_pair] = parse_centermass(Index,SFr,Wid,fitting)

%This is a special function, it takes the vector Index which has
%the structure [0 0 0 1 1 1 0 ... 0 1 0 ... 0]. This vector was obtained
%by coincidence detection of certain events (lower and upper threshold
%crossing for threshold detection, and the appearance of coefficients at
%different scales for wavelet detection).
%The real challenge here is to merge multiple 1's that belong to the same
%spike into one event and to locate that event



Merge = mean(Wid);      %[ms] merge spikes that are closer than Merge, since
%it is likely they belong to the same spike
Merge=0;
% Merge = round(Merge * SFr);


Index([1 end]) = 0;   %discard spikes located at the first and last samples

ind_ones = find(Index == 1);    %find where the ones are
start_stop_pair=[];
if isempty(ind_ones)
    TE = [];
else
    temp = diff(Index);  %there will be 1 followed by -1 for each spike
    N_sp = sum(temp == 1); %nominal number of spikes
    
    lead_t = find(temp == 1);  %index of the beginning of a spike
    lag_t = find(temp == -1);  %index of the end of the spike
    
    for i = 1:N_sp
        
        tE(i) = center_mass(lead_t(i), lag_t(i),fitting);
    end
    
    i = 1;        %initialize counter
    while 0 < 1
        if i > (length(tE)-1)
            
           
                start_stop_pair(i,:)=[lead_t(i),lag_t(i)];
            
            break;
        else
            Diff = tE(i+1) - tE(i);
            
            %             if Diff < Refract & Diff > Merge
            %                 tE(i+1) = [];
            %                 lead_t(i+1)=[];
            %                 lag_t(i+1)=[];
            %                 %discard spike too close to its predecessor
            %             elseif Diff <= Merge
            if Diff <= Merge
                tE(i) = ceil(mean([tE(i) tE(i+1)]));  %merge
                tE(i+1) = [];
                lead_t(i+1)=[];
                lag_t(i)=[];%discard
            else
                start_stop_pair(i,:)=[lead_t(i),lag_t(i)];
                i = i+1;
                
            end
        end
    end
    TE = tE;
end

fcn = TE;