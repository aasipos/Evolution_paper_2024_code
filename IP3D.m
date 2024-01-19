function [result,V,S]=IP3D(no,umb,folder0)

%Spatial isoperimetric ratio from rotational slices
%Input:
%- no:      Number of file in the CT folder (-1: all files)
%- umb:     Is there umbilicus or not
%- folder0: Folder containing the original scans
%Output:
%- result:  IP3D
%- V:       Volume
%- S:       Surface Area
%Andras A. Sipos
%2023
%%
close all
warning off
if nargin<1 no=1; end               %Number of file in the CT folder
if nargin<2 umb=0; end              %Is there umbilicus or not
if nargin<3 folder0='CT'; end
N=10;   %number slices

currentFolder=pwd;
stlfiles=dir(fullfile(currentFolder,'\', folder0, '*.stl'));
shapeData=importdata(strcat(currentFolder,'\data.txt'));
if no==-1 no=1:length(stlfiles); end

for i=1:length(no)
    addpath(stlfiles(no(i)).folder);
    name=stlfiles(no(i)).name(1:end-4)    
    p=shReader(stlfiles(no(i)).name,shapeData.data(no(i),1:7));
    phi=linspace(0,pi,N);
    dphi=pi/(N-1);
    S=0;
    V=0;
    for j=1:N-1
        pl=createPlane([0,0,0],[sin(phi(j)),cos(phi(j)),0]);
        VBPl_LI = isBelowPlane(p.vertices, pl);
        FBP_LI = VBPl_LI(p.faces);
        inside = removeMeshFaces(p.vertices, p.faces, ~((sum(FBP_LI, 2) > 0 & sum(FBP_LI, 2) < 3)));
        [points,~,dS]=produceCurve(inside,pl,dphi);    
        S=S+dS/2;   
        pts2=[points(:,1)*cos(phi(j))-points(:,2)*sin(phi(j)),points(:,3)];
        DD=1;
        minX=min(pts2(:,1))-DD;maxX=max(pts2(:,1))+DD;
        minY=min(pts2(:,2))-DD;maxY=max(pts2(:,2))+DD;
        rgon = polyshape([0,maxX,maxX,0],[minY,minY,maxY,maxY]);
        lgon = polyshape([0,minX,minX,0],[maxY,maxY,minY,minY]);
        shp = boundary(pts2(:,1),pts2(:,2),0.30);
        pgon = polyshape(pts2(shp,1),pts2(shp,2));
        pgon1= intersect(pgon,rgon);
        pgon2= intersect(pgon,lgon);
        [xs1,~]=centroid(pgon1);
        [xs2,~]=centroid(pgon2);
        
        V=V+(abs(pgon1.area*xs1)+abs(pgon2.area*xs2))*dphi;
        A2=0;
        if umb==0 
            dA=0;
        else
            [dA,upol]=umbil(points(:,1:2),shp);
        end
        
%         if length(no)==1
%         figure(j)
%         scatter(pts2(:,1),pts2(:,2));
%         hold on
%         plot(pgon)
%         plot(pgon1)
%         plot(pgon2)
%         
%         if dA>0 plot(upol); end
%             hold off               
%         end
%         figure(N+j)
%         hold on
%         for i2=1:2:size(pts2,1)
%             plot(pts2(i2:i2+1,1),pts2(i2:i2+1,2),'k','LineWidth',2);
%         end
%         hold off
    end
    result(i)=6*sqrt(pi)*V/(S^(3/2));
end


function [points,P,S]=produceCurve(p,pl,dphi)
M=size(p.faces,1);
points=[];    %3D data
P=0;
S=0;
for i=1:M
    [~,pInt] = clipConvexPolygon3dHP(p.vertices(p.faces(i,:),:), pl);
    points=[points;pInt];
    if size(pInt,1)>1
        v=pInt(1:end-1,:)-pInt(2:end,:);        
        P=P+sum(vecnorm(v'));
        r1=(pInt(1,1)^2+pInt(1,2)^2)^(1/2);
        r2=(pInt(2,1)^2+pInt(2,2)^2)^(1/2);
        beta=asin((r2*dphi-r1*dphi)/2/vecnorm(v'));
        S=S+(r1*dphi+r2*dphi)/2*vecnorm(v')*cos(beta);
    end
end

function [dA,upol]=umbil(P,shp)
dA=0;
upol=0;
if isPointInPolygon([0,0],P(shp,:))==1 
    P=[0,0;P];
    TR = delaunayTriangulation(P);
    V=[TR.ConnectivityList(TR.ConnectivityList(:,1)==1,:);TR.ConnectivityList(TR.ConnectivityList(:,2)==1,:);TR.ConnectivityList(TR.ConnectivityList(:,3)==1,:)];
    for k=1:size(V,1)        
        dA=dA+1/2*abs((TR.Points(V(k,2),1)-TR.Points(V(k,1),1))*(TR.Points(V(k,3),2)-TR.Points(V(k,1),2))-(TR.Points(V(k,3),1)-TR.Points(V(k,1),1))*(TR.Points(V(k,2),2)-TR.Points(V(k,1),2)));
    end
    PTS=[TR.Points(V(:,1),1),TR.Points(V(:,1),2);TR.Points(V(:,2),1),TR.Points(V(:,2),2);TR.Points(V(:,3),1),TR.Points(V(:,3),2)];
    shp = boundary(PTS(:,1),PTS(:,2),0);
    upol= polyshape(PTS(shp,1),PTS(shp,2));
end