outdata=[];
K=[];
Tc=[];
To=[];
traceNumber=[];
traceNumberC=[];
M={};
c={};
cellcount=1;
Ccellcount=1;
if (exist('savetofilename','var'))
   clear savetofilename
end
if (exist('saveCtofilename','var'))
   clear saveCtofilename
end
for i= 1: 25:256
    if ((256-i)<25)
            figure
      for j=1:(257-i)
             subplot(5,5,j);
             text(0,55,num2str(i+j-1));
             hold on
             plot(data.data.CorrectedRawData(1:3000,i+j-1));
             
set(gca,'ButtonDownFcn','plotdirecttrace');
set(gca,'UserData',i+j-1);
set(gca,'YLim', [-400, 50]);
line([1000,1100,1100,1000,1000],[50,50,-400,-400,50],'color','r');
      end
    else
    figure,
    for j=1:25
       subplot(5,5,j);
       text(0,55,num2str(i+j-1));
       hold on
       plot(data.data.CorrectedRawData(1:3000,i+j-1));

set(gca,'ButtonDownFcn','plotdirecttrace');
set(gca,'UserData',i+j-1);
set(gca,'YLim', [-400, 50]);
line([1000,1100,1100,1000,1000],[50,50,-400,-400,50],'color','r');
    end
    end
end
