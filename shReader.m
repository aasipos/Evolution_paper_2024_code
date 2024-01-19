function p=shReader(filename,trdata,mode)

%Read and prepare scaned data for analysis
%Input:
%   filename:   stl file to be read
%   trdada:     correction of coordiante system 
%   mode:       plot results (mode=1: yes)
%Output:
%   p:          data structures containing the vertices
%    
%2023
%%
if nargin<3 mode=0; end
if nargin<2 trdata=[0,0,0,0,0,0,0];  end

roty=trdata(1);
rotx=trdata(2);
dY=trdata(3);
dX=trdata(4);
mp=[trdata(5),trdata(6)];
theta=trdata(7)*pi;

[v,f]=stlread(filename,true);
p.vertices=v;
p.faces=f;

if abs(roty)>0 
    RotY=createRotationOy(roty/180*pi); 
    p.vertices=transformPoint3d(p.vertices, RotY);
end
if abs(rotx)>0 
    RotX=createRotationOx(rotx/180*pi); 
    p.vertices=transformPoint3d(p.vertices, RotX);
end

TRANS = createTranslation3d(dX, dY, 0);
p.vertices=transformPoint3d(p.vertices, TRANS);

pl = createPlane([0,0,0], [sin(theta),cos(theta),0], [0,0,10]);
[cp3,cp2]=shSlicer(p,pl,0,0);

if mode==1
    plx= createPlane([0,0,0], [10,0,0], [0,0,10]); 
    v1x = projPointOnPlane(p.vertices, plx);
    ply= createPlane([0,0,0], [0,10,0], [0,0,10]); 
    v1y = projPointOnPlane(p.vertices, ply);
    plz= createPlane([0,0,0], [10,0,0], [0,10,0]); 
    v1z = projPointOnPlane(p.vertices, plz);
    clf
    subplot(2,2,1)
    patch('Faces',p.faces,'Vertices',v1x(:,[1,3]),'FaceColor','none','EdgeColor','m')
    hold on
    plot([0,0],[-20,20],'k')
    scatter(cp3(:,1),cp3(:,3),20,'r','filled')
    axis equal

    subplot(2,2,2)
    patch('Faces',p.faces,'Vertices',v1y(:,[2,3]),'FaceColor','none','EdgeColor','m')
    hold on
    plot([0,0],[-20,20],'k')
    scatter(cp3(:,2),cp3(:,3),20,'r','filled')
    axis equal

    subplot(2,2,3)
    patch('Faces',p.faces,'Vertices',p.vertices,'FaceColor','none','EdgeColor','m')
    hold on
    scatter3(cp3(:,1),cp3(:,2),cp3(:,3),20,'r','filled')
    axis equal

    subplot(2,2,4)
    hold on
    scatter(cp2(:,1),cp2(:,2),2,'k','filled')    
    scatter(mp(1),mp(2),20,'r','filled')
    axis equal
end
