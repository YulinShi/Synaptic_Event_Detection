l=get(gco,'UserData');

tempC=conv(data.data.CorrectedRawData(1100:2600,l),m);

figure
plot(data.data.CorrectedRawData(1100:2600,l));
hold on
plot(tempC/1000,'k');

button = questdlg(['Do you want to save C Value?'],'Save?', 'Yes','No','Yes');
if strcmp(button,'Yes')
c{Ccellcount}=tempC;
Ccellcount=Ccellcount+1;
traceNumberC=[traceNumberC,l];
if (exist('saveCtofilename','var'))
save(saveCtofilename,'traceNumberC','c','-append');
else
   clear fileName
   [fileName,pathName]=uiputfile(['trace_C.mat'], 'Create data m-file:');
   if isequal(fileName,0) || isequal(pathName,0)
    disp('User selected Cancel')
   else
    disp(['User selected',fullfile(pathName,fileName)])
   end
   saveCtofilename=fullfile(pathName,fileName);
   save(saveCtofilename,'traceNumberC','c');
end
end 