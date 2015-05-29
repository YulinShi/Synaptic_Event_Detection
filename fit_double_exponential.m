function [Model,xopt] = fit_double_exponential(Data,SynType,PlotFlag, SFr)

%This function looks at the time series given by Data and
%it fits a double exponential of the form y = K*[exp(-t/to)-exp(-t/tc)]
%where to < tc and both time constants are positive
%SynType = 'inh' or 'exc'
%for inhibitory synapse use: K > 0
%for excitatory synapse use: K < 0;
%PlotFlag = 1 --> generate plot, otherwise don't

%sampling frequency
SFr = 10000; %[Hz]
Niter = 10;


t = [0:length(Data)-1]*1000/SFr;    %[milliseconds]
figure
plot(t,Data)
[t1, t2] = ginput(2);
s1 = round(t1)/1000*SFr;


switch num2str(SynType)
    case 'inh'
        lb(1) = eps;
        x0(1) = 1;
    case 'exc'
        lb(1) = -eps;
        x0(1) = -1;
    otherwise
        error('unknown synapse type')
end

lb(2:3) = eps*[1 1];
ub = [];

x0(2:3) = rand(1,1)*[1 2];

dData = Data(s1(1):s1(2)) - max(Data(s1(1):s1(2)));

xopt = x0;
Error = double_exp(x0,dData,SFr);
OPT = optimset('lsqnonlin');
OPT = optimset('Display','off');

Error = sum(Error.^2)/length(dData);
disp(['Initial Error: ' num2str(Error)])
for i = 1:Niter
    [x, rnorm] = lsqnonlin(@(x) double_exp(x,dData,SFr),x0,lb,ub,OPT);
    Err = rnorm()/length(dData);
    if Err < Error
        disp(['better solution found;   Error: ' num2str(Err)])
        Error = Err;
        xopt = x;
    else
        disp('...')
    end
    x0 = x + rand(1,3);
end
disp(['Final Error: ' num2str(Error)])
disp(['K=' num2str(xopt(1)) '; to=' num2str(xopt(2)) ...
    '; tc=' num2str(xopt(3))])

 figure 
 subplot(2,1,1)
    t = [0:length(dData)-1]*1000/SFr;    %[milliseconds]
    plot(t,dData)
    hold on
    Model = xopt(1)*[exp(-t./xopt(2))-exp(-t./xopt(3))];
    plot(t,xopt(1)*[exp(-t./xopt(2))-exp(-t./xopt(3))],'r')
    xlabel('time [ms]')

