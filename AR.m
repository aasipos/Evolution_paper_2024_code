function AR=AR(no)

%Collecting the AspectRatio for the cross-sections and write data to an XLS file
%Input:
%- no:      Number of file in the SecByHandResults folder (-1: all files)
%Output:
%- AR:      AspectRatio for all slices
%Andras A. Sipos
%2023
%%
if nargin<1 no=-1; end
currentFolder=pwd;
matfiles=dir(fullfile(currentFolder, '\SecByHandResults','*.mat'))
if no==-1 no=1:length(matfiles); end

for i=1:length(no)
    matfiles(no(i)).name
    load(strcat(currentFolder,'\SecByHandResults\',matfiles(i).name))
    %sec=data.sec;
    N=length(sec2);
    vertices=[];
    S=0;
    V=0;
    for j=1:floor(N)
        minx=min(sec2(j).cp2(:,1));
        maxx=max(sec2(j).cp2(:,1));
        miny=min(sec2(j).cp2(:,2));
        maxy=max(sec2(j).cp2(:,2));                
        AR(i,j)=abs(maxy-miny)/abs(maxx-minx);
    end
    
    res(i).name=matfiles(no(i)).name;
    res(i).AR=AR(i,:);   
end

writetable(struct2table(res), 'AR.xlsx');
