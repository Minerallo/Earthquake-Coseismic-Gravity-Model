function varargout=Sumatra2005(varargin)
Model=varargin{1};
dxgrid=varargin{2};
longT2=varargin{3};
latT2=varargin{4};
coef_tile=varargin{5};
button=varargin{6};
outdata=varargin{7};
coef_tile=80;

load ('GeomSumatra05.mat')

OPEN=0;     % OPEN   : dislocation in tensile component [m]
RHO=2670;   % RHO    : density of the medium [kg/m^3]
RHOP=2670;  % RHOP   : optional density of the cavity-filling matter [kg/m^3]; takes RHO value if not specified.
mu=30e9;    % Shear modulus (Pa) for Leonard2010 model
BETA=0;
NU=0.25;

numMwint=36;
rake1a=rakesubseg;
D=slipsubseg./1e2;
strike1a=strikesubseg;
dip1a=dipsubseg;
depth1a=depthsubseg.*1e3;

for i=1:numel(dip1a)
W(i)=16000;
L(i)=16000;
end

%% Parameters to define dx, x
dx_bestDSREV=importdata('dx_bestDSREV.mat');EDGDSREV=importdata('EDGmatDSREV.mat');
dx_bestDSNORM=importdata('dx_bestDSNORM.mat');EDGDSNORM=importdata('EDGmatDSNORM.mat');
dx_bestSS=importdata('dx_bestSS.mat');EDGSS=importdata('EDGmatSS.mat');
Mwint=0:0.25:9 ;%used later for the interpolation of dx depending on the Magnitude

[X,Y,zone]=ll2utm(Latsubseg,Lonsubseg);

zone=zeros(size(X))+zone;

%% Building the grid for all earthquakes
longT=-180:dxgrid:180;
latT=-90:dxgrid:90;
%Construction of the total grid point matrix
matot=zeros(numel(latT),numel(longT));
if outdata==2
    matot2=zeros(numel(latT),numel(longT));
end
[longT2,latT2]=meshgrid(longT,latT);

%% Okubo calculation
for R=1:numel(dip1a)
      for j=Model
    %disp([num2str(R),'/',num2str(numel(W)),' Geometry calculation for Sumatra2005'])
        %% Define the grid size and resolution for each earthquake
    % Change dataset depending on the Rake
        if (rake1a(R) <=135 && rake1a(R) >=45)% Dip-slip |reverse fault
            EDGmat=EDGDSREV(j,numMwint);disp ('la1');dx_best=dx_bestDSREV(j,numMwint);
        elseif (rake1a(R) <=-45 && rake1a(R) >=-135)% Dip-slip |normal fault
            EDGmat=EDGDSNORM(j,numMwint);disp ('la2');dx_best=dx_bestDSNORM(j,numMwint);
        elseif (rake1a(R) >-45 && rake1a(R)<=45)%Strike slip
            EDGmat=EDGSS(j,numMwint);disp ('la3');dx_best=dx_bestSS(j,numMwint);
        elseif (rake1a(R)>135 && rake1a(R)<=180)%Strike slip
            EDGmat=EDGSS(j,numMwint);disp ('la4');dx_best=dx_bestSS(j,numMwint);
        elseif (rake1a(R) <-135 && rake1a(R)>=-180)%Strike slip
            EDGmat=EDGSS(j,numMwint);disp ('la5');dx_best=dx_bestSS(j,numMwint);
        end

            dxint1(R)=dx_best;
            xmin1(R)=-EDGmat;
            xmax1(R)=EDGmat;
      end
end
for R=1:numel(dip1a)
    disp([num2str(R),'/',num2str(numel(W)),' Okubo calculation for Sumatra2005'])
    %% Define Est North-vector for Okubo
        E=xmin1(R):dxint1(R):xmax1(R);
        N=xmin1(R):dxint1(R):xmax1(R);
%     E=xmin:dxint:xmax;
%     N=xmin:dxint:xmax;
    [E1,N1]=meshgrid(E,N);
    %Georeferencement of the earthquake
    lon1=E+X(R);
    lat1=N+Y(R);
    %Okubo funtion
         if D(R)==0
            dH = zeros (size(E1,1),size(E1,2));
            dG = zeros (size(E1,1),size(E1,2));
        else
            if strcmp(button,'No')
                [dG,dH]=okubo92(E1,N1,depth1a(R),strike1a(R),dip1a(R),L(R),W(R),rake1a(R),D(R),OPEN,RHO,RHOP,BETA,NU);
            else
                [dG,dH]=okubo92(E1,N1,depth1a(R),strike1a(R),dip1a(R),L(R),W(R),rake1a(R),D(R),OPEN,RHO,RHOP);
            end
        end
[lon2,lat2]=meshgrid(lon1,lat1);
[LA2,LO2] = utm2ll(lon2(1:numel(lon2)),lat2(1:numel(lat2)),zone(R));
indlocallon=find(wrapTo360(longT2(1,:))>=wrapTo360(min(LO2)) & wrapTo360(longT2(1,:))<=wrapTo360(max(LO2)));
% indlocallon=find(wrapTo360(longT2(1,:))>=min(LO2) & wrapTo360(longT2(1,:))<=max(LO2));
indlocallat=find(latT2(:,1)>=min(LA2) & latT2(:,1)<=max(LA2));
LO2 = wrapTo180(LO2);
if outdata==1
    dvar1=dH;
elseif outdata==0
    dvar1=dG;
elseif outdata==2
    dvar1=dG;
    dvar2=dH;
end
Gqlocal = gridfit(LO2,LA2,dvar1(1:numel(dvar1)),longT2(1,indlocallon),latT2(indlocallat,1),'tilesize',round(numel(longT)/coef_tile),'overlap',0.25);
Gqlocal(isnan(Gqlocal)==1)=0;
Gq=zeros(size(longT2));
Gq(indlocallat,indlocallon)=Gq(indlocallat,indlocallon)+Gqlocal;
    matot=matot+Gq;
    if outdata==2
        Gqlocal2 = gridfit(LO2,LA2,dvar2(1:numel(dvar2)),longT2(1,indlocallon),latT2(indlocallat,1),'tilesize',round(numel(longT)/coef_tile),'overlap',0.25);
Gqlocal2(isnan(Gqlocal2)==1)=0;
Gq2=zeros(size(longT2));
Gq2(indlocallat,indlocallon)=Gq2(indlocallat,indlocallon)+Gqlocal2;
    matot2=matot2+Gq2;
    end
end
Gq=matot;
varargout{1}=Gq;
if outdata==2
Gq2=matot2;varargout{2}=Gq2;
else
    Gq2=0;varargout{2}=Gq2;
end


% m_proj('Miller Cylindrical','long',[-85 60],'lat',[-45 25]);
% % fig points
% set(figure,'Position',[1 1 2000 1500]) % [coin_x coin_y hauteur(px) largeur(px)]
% set(gcf,'PaperPositionMode','auto')
% m_pcolor(longT,latT,matot.*1e8),shading flat;
% colormap(jet);
% colormap(b2r((min(min(matot.*1e8))),max(max(matot.*1e8))));
% m_coast('color','k','LineWidth',1);
% m_grid('box','fancy','tickdir','out');colorbar;
% title('Chile Beta=0');

% 
% save('matot','matot')
% save('longT','longT')
% save('latT','latT')



