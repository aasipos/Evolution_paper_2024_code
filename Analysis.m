function [ps,sec,pars,filename]=Analysis(no,steps,mode)
%The main file for the analysis of a preprocessed shell
%The simplifed midsurfaces is needed in the path\midsurface folder
%
%Input:
%- no:      The number of file to be analyzed (assuming alphabetic order)
%- steps:   Number of cross-sections to be computed
%- mode:    Make png file from the figure
%Output:
%- ps:      spatial data of the shell
%- sec:     sections
%- pars:    relevant parameters
%- filename:filename
%
%The code uses the geom3d toolbox:
%https://nl.mathworks.com/matlabcentral/fileexchange/24484-geom3d?requestedDomain=
%
%Andras A. Sipos
%2023
%%
if nargin<3 mode=1; end    %make png file from the figure
if nargin<2 steps=21; end  %number of sections computed
currentFolder=pwd;

%analysis of preprocessed files
warning('off')
matfiles=dir(fullfile(currentFolder, 'midsurface', '*.mat'));
filename=matfiles(no).name
ps=load(strcat(currentFolder,'\midsurface\',matfiles(no).name));
pars0=ps.pars;

if length(pars0)>= 8    %parameter for curve reconstruction
    R=pars0(8);
else
    R=5;
end      
pars(1)=pars0(7)*pi;
pars(2)=sign(pars0(7))*pi/5;
pars(3)=steps;
pars(4)=pars0(5);
pars(5)=pars0(6);
pars(6)=R;

sec=CSTracer(ps,pars);
