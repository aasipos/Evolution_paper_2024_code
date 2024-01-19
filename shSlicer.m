function [cp3,cp2,wp]=shSlicer(p,pl,mode,mod)

%Intersection and projection to the plane pl
%Andras. A. Sipos
%2023
%%
if nargin<4 mod=1; end   %slicing the patch p with the plane pl

[cp3,wp] = crossingPoints(p.vertices,p.faces,pl,mod);
cp3= projPointOnPlane(cp3, pl);
cp2= planePosition(cp3, pl); %projection onto 2D

try
if mode==1
    zmin=min(cp2(:,2));
    zmax=max(cp2(:,2));
    zs= linspace(zmin,zmax,100);
    cp2=[cp2;zeros(100,1),zs'];
    cp3=[cp3;zeros(100,1),zeros(100,1),zs'];
end
catch
end

function [cp,wp]=crossingPoints(v,f,plane,mod)
% computing the edge list
e = meshEdges(f);
edges = [ v(e(:,1), :) v(e(:,2), :) ];

% identify which edges cross the mesh
inds = isBelowPlane(v, plane);
edgeCrossInds = find(sum(inds(e), 2) == 1);
cp = intersectEdgePlane(edges(edgeCrossInds, :), plane);

if mod==1
    nFaces = length(f);
    faceEdges = cell(1, nFaces);
    nCrossEdges = length(edgeCrossInds);
    crossEdgeFaces = zeros(nCrossEdges, 2);

    for iEdge = 1:length(edgeCrossInds)
        edge = e(edgeCrossInds(iEdge), :);
        indFaces = find(sum(ismember(f, edge), 2) == 2);
        crossEdgeFaces(iEdge, :) = indFaces;
    
        for iFace = 1:length(indFaces)
            indEdges = faceEdges{indFaces(iFace)};
            indEdges = [indEdges iEdge];
            faceEdges{indFaces(iFace)} = indEdges;
        end
    end
    
    k=1;
    for i=1:length(crossEdgeFaces)
        w0=[faceEdges{crossEdgeFaces(i,1)}];
        wp(k,:)=w0;
        w0=[faceEdges{crossEdgeFaces(i,2)}];
        wp(k+1,:)=w0;
        k=k+2;
    end
    wp=unique(wp,'rows');
else
    wp=-1;
end   