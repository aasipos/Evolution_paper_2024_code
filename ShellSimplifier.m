function ptCloud=ShellSimplifier(v,mode,gap)

%Simplification of the data, and the approximation of the midsurface 
%Andras A. Sipos
%2023
%%
if nargin<3 gap=1; end
ptCloud = pointCloud(v);

if mode==1
    M=200;
    R=reshape(randi(ptCloud.Count,M),M*M,1);
    R=unique(R,'rows');
    SR=max(size(R));
else
    SR=floor(ptCloud.Count/2);
    R=1:2:ptCloud.Count;
end

w=[];
for i=1:SR
        idx=R(i);
        dist=1.0;
        pt0=ptCloud.Location(idx,:);
        k=0;
        while dist>0.05*gap && k<50
            indices= findNeighborsInRadius(ptCloud,pt0,gap);
            pts=ptCloud.Location(indices,:);
            pt1=mean(pts);
            dist=norm(pt1-pt0,2);
            pt0=pt1;
            k=k+1;
        end
        if k<50
            w=[w;pt0]; 
        else 
            display('failed recognition')
        end
end

v=unique(w,'rows');
ptCloud=pointCloud(v);
ptCloud=pcdenoise(ptCloud);