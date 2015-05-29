function newData=ReformData(oldData,Mappatern)
global Row
Row=size(Mappatern,1);
k=1;
if iscell(oldData)
newData={};
for i = 1:size(Mappatern,1)
    for j= 1:size(Mappatern,2)
       newData{k}=oldData{Mappatern(i,j) };
        k=k+1;
    end
end
else
    newData=[];
for i = 1:size(Mappatern,1)
    for j= 1:size(Mappatern,2)
       newData(:,k)=oldData(:,Mappatern(i,j));
        k=k+1;
    end
end
end