function  [idxh, hole, mp, Ar]=HoleSeeker(cp2,mp,wp,R)

%cp2    point set in 2D
%mp     guess on a point inside the hole
%%
if nargin<4  R=10; end    %resolution of the grid.

idx0=boundary(cp2,0.0);                %Boundary identification
idx1=boundary(cp2,1.0);                %Boundary identification
Ach=polyarea(cp2(idx0,1),cp2(idx0,2)); %Total area covered by the convex hull of the cross-section
At=polyarea(cp2(idx1,1),cp2(idx1,2));  %Total area covered by the cross-section
Ar=At/Ach;
dx=max(cp2(:,1))-min(cp2(:,1));
dy=max(cp2(:,2))-min(cp2(:,2));
xs=linspace(min(cp2(:,1))-dx*0.10,max(cp2(:,1))+dx*0.10,R);
ys=linspace(min(cp2(:,2))-dy*0.10,max(cp2(:,2))+dy*0.10,R);
k=1;
for i=1:R
    for j=1:R
        gp(k,:)=[xs(i),ys(j)];
        k=k+1;
    end
end

ID2=size(cp2,1)+1:size([cp2;gp;mp],1);
DT=delaunayTriangulation([cp2;gp;mp],wp);
tri = DT.ConnectivityList;

hole=[];
idxh=[];
gpts=[ID2(end)];
k=0;
while size(gpts,1)>0 && k<200
    rows=any(ismember(tri,gpts),2);
    pts=tri(rows,:);
    tri(rows,:)=[];
    pts=reshape(pts,size(pts,1)*3,1);
    pts=unique(pts);
    pts=setdiff(pts,gpts);
    ind=find(pts<=size(cp2,1));
    idxh=[idxh;pts(ind)];
    gpts=setdiff(pts,pts(ind));  
    k=k+1;
end

if k>190 display('error in hole recognition'); end

if size(idxh,1)>1
    idxh=unique(idxh,'rows');
    xc=mean(cp2(idxh,1));
    yc=mean(cp2(idxh,2));

    hole=cp2(idxh,:);
    [~,~,~,indout]=points2contour(hole(:,1)-xc,hole(:,2)-yc,1,'cw');
    idxh=idxh(indout(1:end));
    hole=cp2(idxh,:);

    xc=mean(cp2(idxh,1));
    yc=mean(cp2(idxh,2));
    mp=[xc,yc];
end