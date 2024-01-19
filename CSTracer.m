function sec = CSTracer(p,pars)

%Tracing the consequtive cross-sections of the shell
%Andras A. Sipos
%2023
%%
if nargin<2
   theta0=3*pi/4;
   dtheta=-pi/25;
   steps=90;
   mp=[-100,50];      
   R=10;
else
   theta0=pars(1);
   dtheta=pars(2);
   steps=pars(3);
   mp=[pars(4),pars(5)];
   R=pars(6);   
end

SliceMode=1;

for i=1:steps
    theta=theta0+(i-1)*dtheta;
    pl = createPlane([0,0,0], [sin(theta),cos(theta),0], [0,0,10]);
    [cp3,cp2, wp]=shSlicer(p,pl,SliceMode);
    
    %Identification of the actual, closed curve
    [idxh, ~, mp, Ar]=HoleSeeker(cp2,mp,wp,R);
    
    sec(i).mp=mp;
    sec(i).cp3=cp3;
    sec(i).cp2=cp2;
    sec(i).pl=pl;
    sec(i).theta=theta;
    sec(i).idxh=idxh;
    sec(i).Ar=Ar;
    
    %Identification of the actual, open curve
    %In some cases automatic recognation fails
    try
    [C,sp,curX,curY,F,IP]=fourierFit(sec(i));
    sec(i).spect=C;
    sec(i).trigs=F;
    sec(i).IP=IP;
    sec(i).sp2=sp;
    sec(i).sp3=[sin(theta)*sp(:,1),cos(theta)*sp(:,1),sp(:,2)];
    sec(i).curX=curX;
    sec(i).curY=curY;
    sec(i).mp=[sp(2,1),sp(2,2)];    
    end
end
