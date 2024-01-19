function IP2D

%Collecting the IP2D data for the cross-sections and write data to an XLS
%file
%Andras A. Sipos
%2023
%%
warning('off')
currentFolder=pwd;
matfiles=dir(fullfile(currentFolder, '\SecByHandResults','*.mat'))
no=1:length(matfiles);

str=["angles"];
for i=1:length(no)
    str2=[matfiles(i).name];
    str2=erase(str2,'__.mat')
    str=[str;str2];
    load(strcat(currentFolder,'\SecByHandResults\',matfiles(i).name),'sec2')    
    for j=1:size(sec2,2)
        M(j)=(j-1)*36;
        try
            N(i,j)=sec2(j).IP;
        catch
            N(i,j)=-1;
        end
    end
end

T=table([M;N],'rownames',str);
writetable(T,'IpData.xlsx',"WriteRowNames",true,"WriteVariableNames",false)