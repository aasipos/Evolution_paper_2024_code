function Preprocess(no,folder0)

%preprocessing scanned gastropd shell data.
%Inout:
%- 'folder0' contains the *.stl input AND a data.txt file for correction
%- 'no' is the number of the file to be processed. 
%Note that no=-1 processes all the files in the folder
%
%The code uses the geom3d toolbox:
%https://nl.mathworks.com/matlabcentral/fileexchange/24484-geom3d?requestedDomain=
%
%The STL data converted to *.mat file (without simlification) is written to
%the folder path\surface
%
%The simpilified midsurface, both in *.STL and *.mat is written to the
%folder path\midsurface
%
%Andras A. Sipos
%2023
%%
currentFolder=pwd;
warning('off')
mkdir('surface')
mkdir('midsurface')
warning('on')

if nargin<2 folder0='CT'; end
if nargin<1 no=1; end

stlfiles=dir(fullfile(currentFolder, folder0, '*.stl'));
shapeData=importdata(strcat(currentFolder,'\data.txt'));
if no==-1 no=1:length(stlfiles); end

for i=1:length(no)
    %save surface data in *.mat file in the \surface folder
    addpath(stlfiles(no(i)).folder);
    name=stlfiles(no(i)).name(1:end-4)    
    p=shReader(stlfiles(no(i)).name,shapeData.data(no(i),1:7));
    p.pars=shapeData.data(no(i),:)
    filename=strcat(currentFolder,'\surface\',name,'_.mat');
    save(filename, '-struct','p');    

    %save surface data in *.mat file in the \surface folder
    gap=max(max(p.vertices(:,[1,2]))-min(p.vertices(:,[1,2])))/50;
    ptCloud=ShellSimplifier(p.vertices,1,gap);
    v=ptCloud.Location;
    f=MyCrustOpen(v);
    r.vertices=v;
    r.faces=f;
    r.pars=shapeData.data(no(i),:);
    filename=strcat(currentFolder,'\midsurface\',name,'_.mat');
    save(filename, '-struct','r');    
    filename=strcat(currentFolder,'\midsurface\',name,'_.stl');
    stlwrite(filename,f,v)    
end