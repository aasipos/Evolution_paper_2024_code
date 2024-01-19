function [C,sp,curX,curY,F,Ip]=fourierFit2(sec)

%Analysis of curve segment
%requires the Chebfun Toolbox
%Andras A. Sipos
%2021

if isfield(sec,'sp2s')==0
    %alpha shape associated with the point set
    shp0=alphaShape(sec.cp2(sec.idxh,1),sec.cp2(sec.idxh,2),200000);         %convex hull
    bf0 =boundaryFacets(shp0);
    pc = criticalAlpha(shp0,'one-region');
    shp=alphaShape(sec.cp2(sec.idxh,1),sec.cp2(sec.idxh,2),1.05*pc);        %alpha hull
    bf = boundaryFacets(shp);

    bf1=unique(reshape(bf,size(bf,1)*size(bf,2),1));
    bf01=unique(reshape(bf0,size(bf0,1)*size(bf0,2),1));

    %polygons from the alpha shape
    pgon0 = polyshape(shp0.Points(bf01,1),shp0.Points(bf01,2));
    pgon1 = polyshape(shp.Points(bf1,1),shp.Points(bf1,2));
    A=area(pgon1);
    [sx,sy]=centroid(pgon1);
    sp=[sx,sy];

    %curv segment with only positive curvature
    final=polyshape;
    polyout = subtract(pgon0,pgon1);
    polyout2= regions(polyout);
    for i=1:polyout.NumRegions
        if area(polyout,i)>0.005*A
        pgon2=polyout2(i);
        final=union(final,pgon2);
        end
    end

    SP=setdiff(shp0.Points(bf01,:),final.Vertices,'rows');
    sp(2,:)=mean(SP);

    %Representation in polar coordiantes
    phi=atan2(SP(:,1)-sp(2,1),SP(:,2)-sp(2,2));
    [phi,ind]=sort(phi);
    SP=SP(ind,:)-sp(2,:);
    rad=vecnorm(SP')';

    %excluding part associated overlap
    [~,i]=max(abs(diff(phi)));
    rad=[rad(i+1:end);rad(1:i-1)];
    phi=[phi(i+1:end);phi(1:i-1)+2*pi];
    SP=[SP(i+1:end,:);SP(1:i-1,:)];

    %removing points along the global Z
    DX=max(SP(:,1))-min(SP(:,1));
    ind=find(abs(SP(:,1)+sp(2,1))<0.05*DX);
    SP(ind,:)=[];
    rad(ind,:)=[];
    phi(ind,:)=[];    
else
    pgon1=polyshape(sec.sp2s);
    A=area(pgon1);
    [sx,sy]=centroid(pgon1);
    sp=[sx,sy];
    
    %n = numel(sec.sp2sc);
    %perim = 0.0;
    %for i = 1:n-1
    %    perim= perim+ sqrt( (sec.sp2sc(i+1,1)-sec.sp2sc(i,1))^2 + (sec.sp2sc(i+1,2)-sec.sp2sc(i,2))^2 );
    %end
    %Ip=4*pi*A/perim^2;
    
    sp(2,:)=mean(sec.sp2sc);
    phi=atan2(sec.sp2sc(:,1)-sp(2,1),sec.sp2sc(:,2)-sp(2,2));
    if max(abs(diff(phi)))>6
        [~,i]=max(abs(diff(phi)));
        if phi(i+1)-phi(i)>0
            phi(1:i)=phi(1:i)+2*pi;
        else
            phi(i+1:end)=phi(i+1:end)+2*pi;
        end
    end
    [phi,ind]=sort(phi);
    sp2sc=sec.sp2sc(ind,:)-sp(2,:);
    rad=vecnorm(sp2sc')';
end
    
%chebfun representation
phie=linspace(min(phi),max(phi),100);
rads=interp1(phi,rad,phie','pchip');
contR=chebfun(rads,[min(phi),max(phi)],'equi');
perim=sum(sqrt(contR.*contR+diff(contR).*diff(contR)));
Ip=4*pi*A/perim^2;
contR=contR/perim;      %scaling for unit length;

%archlength 
archl=cumsum(sqrt(contR.*contR+diff(contR).*diff(contR)));
arce=linspace(0,1,100);
x0=interp1(archl(phie),contR(phie).*sin(phie),arce','pchip');
y0=interp1(archl(phie),contR(phie).*cos(phie),arce','pchip');
xs=chebfun(x0,[0,1],'equi');
ys=chebfun(y0,[0,1],'equi');
kappa=(diff(xs,2).^2+diff(ys,2).^2).^(1/2);

curX=xs*perim+sp(2,1);
curY=ys*perim+sp(2,2);

C=chebcoeffs(kappa,50);
[A,B]=trigcoeffs(kappa,49);
F=[A,[0;B]];

%circle fitting
opts=optimoptions('lsqnonlin','Display','off');
x=lsqnonlin(@(x)circSeeker(x,curX,curY),[sp(2,1),sp(2,2),mean(contR)],[],[],opts);

sp(3,1)=x(1);
sp(3,2)=x(2);
    

function F=circSeeker(x,xs,ys)
d=domain(ys);
pd=linspace(d(1),d(2),1000);
F=[(xs(pd)-x(1)).^2+(ys(pd)-x(2)).^2-x(3)^2]';
