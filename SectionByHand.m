function sec2=SectionByHand(sec,ps,pars,filename0,ik0)

%Rewrite sec data by user identified cross-section
%The graphical interface works as follows: 
%- points can be placed by the left mouse button
%- from point 1 until the key SPACE is pushed the points belong to the
%outer surface
%- points placed after the SPACE is pushed belong to the overalpped surface
%- placing points is finished by pushing the ENTER
%The input can be produced by Analysis.m 
%- ik0  : the number of the cross section to be corrected 
%- ik0=-1: all cross-sections are corrected manually (default)
% 
%Andras A. Sipos
%2023
%%
close all
N=size(sec,2);
sec2=sec;
if nargin<5 ik0=-1; end

for ik=1:N
    if ik0~=-1 ik=ik0; end
    h=figure;
    h.WindowState='maximized';
    hold on
    scatter(sec(ik).cp2(:,1),sec(ik).cp2(:,2),6,'k','filled')    
    scatter(sec(ik).cp2(sec(ik).idxh,1),sec(ik).cp2(sec(ik).idxh,2),3,'r')    
    plot(sec(ik).curX,sec(ik).curY,'g','LineWidth',1)
    axis equal
    title('(e) Cross-section','Interpreter','latex','FontSize',16)
    x=[];
    y=[];
    forever=1;
    csize=0;
    while forever
        [x1, y1] = ginput(1);
        x=[x;x1];
        y=[y;y1];
        scatter(x1,y1,20,'b','filled');
        cr=get(h,'CurrentCharacter');
        KeyPressed = get(h,'CurrentCharacter');
        if KeyPressed==' '
            spcv = cscvn( [x, y].' );
            fnplt( spcv );
            [x_ad] = keszito( spcv );
            csize=max(size(x_ad));   
            break;
        end
    end
    
    while forever
        [x1, y1] = ginput(1);
        x=[x;x1];
        y=[y;y1];
        scatter(x1,y1,20,'b','filled');
        cr=get(h,'CurrentCharacter');
        KeyPressed = get(h,'CurrentCharacter');
        if KeyPressed==char(13)
             break
        end       
    end
        
    x(end+1)=x(1);
    y(end+1)=y(1);
    spcv = cscvn( [x, y].' );
    [ x_ad,y_ad] = keszito( spcv );
    scatter(x_ad,y_ad,'r','filled');
    hold off 
    sec2(ik).sp2s=[x_ad,y_ad];
    sec2(ik).sp2sc=[x_ad(1:csize),y_ad(1:csize)];   
    [C,sp,curX,curY,F,IP]=fourierFit(sec2(ik));
    
    sec2(ik).spect=C;
    sec2(ik).trigs=F;
    sec2(ik).IP=IP;
    sec2(ik).sp2=sp;
    sec2(ik).sp3=[sin(sec2(ik).theta)*sp(:,1),cos(sec2(ik).theta)*sp(:,1),sp(:,2)];
    sec2(ik).curX=curX;
    sec2(ik).curY=curY;
    sec2(ik).mp=[sp(2,1),sp(2,2)];    
    close;
    ik
    if ik==ik0 break; end
end

if nargin>1
    if ik0==-1 name=filename0(1:end-4);
    else name=filename0(1:end-5);
    end
    currentFolder=pwd;
    warning('off')
    mkdir('SecByHandResults')
    warning('on')
    filename=strcat(pwd,'\SecByHandResults\',name,'_.mat');
    save(filename,'ps','sec','sec2','pars','filename0');
end

function [ x_ad,y_ad,kpi ] = keszito(f)
szam=0;
x_ad=zeros(3000,1);
y_ad=zeros(3000,1);
kpi=zeros(f.pieces,1);
for i=1:f.pieces
    kpi(i)=szam+1; 
    for t=0:0.1:abs(f.breaks(1,i)-f.breaks(1,i+1))
        szam=szam+1;
        x_ad(szam,:)=f.coefs(2*i-1,1)*t^3+f.coefs(2*i-1,2)*t^2+f.coefs(2*i-1,3)*t^1+f.coefs(2*i-1,4)*t^0;
        y_ad(szam,:)=f.coefs(2*i,1)*t^3+f.coefs(2*i,2)*t^2+f.coefs(2*i,3)*t^1+f.coefs(2*i,4)*t^0;
    end
    if i==f.pieces && t~=abs(f.breaks(1,i)-f.breaks(1,i+1)) % this part is for the last point (curve length / t is not an integer usually)
        t=abs(f.breaks(1,i)-f.breaks(1,i+1));
        szam=szam+1;
        x_ad(szam,:)=f.coefs(2*i-1,1)*t^3+f.coefs(2*i-1,2)*t^2+f.coefs(2*i-1,3)*t^1+f.coefs(2*i-1,4)*t^0;
        y_ad(szam,:)=f.coefs(2*i,1)*t^3+f.coefs(2*i,2)*t^2+f.coefs(2*i,3)*t^1+f.coefs(2*i,4)*t^0;
    end
end
kpi(i+1)=szam;
x_ad=x_ad(1:szam,1);
y_ad=y_ad(1:szam,1); 