function varargout=Sumatra(varargin)
%Fault parameters from Banerjee, 2007
Model=varargin{1};
dxgrid=varargin{2};
longT2=varargin{3};
latT2=varargin{4};
coef_tile=varargin{5};
button=varargin{6};
outdata=varargin{7};
%% Earthquake catalogue
load('GeomSumatra.mat');
%% Input model parameters
OPEN=0;     % OPEN   : dislocation in tensile component [m]
RHO=2670;   % RHO    : density of the medium [kg/m^3]
RHOP=2670;  % RHOP   : optional density of the cavity-filling matter [kg/m^3]; takes RHO value if not specified.
mu=30e9;    % Shear modulus (Pa) for Leonard2010 model
BETA=0;
NU=0.25;
%The seismic moment of the preferred model is M0 =
%7.62 * 10^22 N m, corresponding to Mw 9.22. (Again, we
%use the formula Mw  (2/3) log10 M0 [dyne cm]  10.7,
%the definition given by Hanks and Kanamori [1979].)
%% Define your area
%{'75','110','-10','25'};
%% Parameters to define dx, x
dx_bestDSREV=importdata('dx_bestDSREV.mat');EDGDSREV=importdata('EDGmatDSREV.mat');
dx_bestDSNORM=importdata('dx_bestDSNORM.mat');EDGDSNORM=importdata('EDGmatDSNORM.mat');
dx_bestSS=importdata('dx_bestSS.mat');EDGSS=importdata('EDGmatSS.mat');
Mwint=0:0.25:9 ;%used later for the interpolation of dx depending on the Magnitude
%% Sumatra block lower edge corner
longitudeSuma(9)=longitudeSuma(9)-0.6;%Correction Banerjee probleme de longitude
[X,Y,zone]=ll2utm(latitudeSuma,longitudeSuma);
zone=zeros(size(X))+zone;
Mw = 9.22;
Mo=7.62 * 10^22;
Modyn=Mo*1e7 ;%[dyn cm]
Mwvec = (2/3)*log10(Modyn)-10.7;
%% Geometry
dip1a=DipSuma;
depth1a=(d1Suma+d2Suma)/2*1e3;
W=(d1Suma-d2Suma)*1e3;%WidthSuma*1e3;
L=LengthSuma*1e3;
strike1a=StrikeSuma;
RakeSuma(isnan(RakeSuma)==1)=0;
%% Slip
D=[2.3 2.3 2.3 2.3 6.6 5.2 7.8 16.6 0 19.4 0 15.5 0.5 15.5 9.2];%ModelC Banerjee
% u1A (m) u1B (m) u1C (m) u1D (m) u1E (m) u1F (m) u1G (m) u1H (m) u2A (m)u2B (m) u2C (m) u3B (m) u3C (m) u3D (m)
%% Building the grid for all earthquakes
longT=-180:dxgrid:180;
latT=-90:dxgrid:90;
%Construction of the total grid point matrix
matot=zeros(numel(latT),numel(longT));
if outdata==2
    matot2=zeros(numel(latT),numel(longT));
end

for R=1:numel(WidthSuma)
    if lambdaSuma(R) == 1
        rake1a(R)= 105+RakeSuma(R) ;
    elseif lambdaSuma(R)==2
        rake1a(R)= 115 ;
    elseif lambdaSuma(R)==3;
        rake1a(R)= 90 ;
    end
end
%New loop allow to fix  a constant dx and size of the grid calculation for okubo
%Without size(grid,2) = 1.0576e+03 calcul for the first iteration TOO
%MUCH!!! let's reduce at <100 by taking the lower resolution for all plan
%assuming that there is no such big difference between the smallest plan and the biggest for sumatra
for R=1:numel(WidthSuma)
    %% Define the grid size and resolution for each earthquake
    % Change dataset depending on the Rake
    for j=Model
        if (rake1a(R)<=135 && rake1a(R)>=45)% Dip-slip |reverse fault
            EDGmat=EDGDSREV(j,end);disp ('la1');dx_best=dx_bestDSREV(j,end);
        elseif (rake1a(R) <=-45 && rake1a(R)>=-135)% Dip-slip |normal fault
            EDGmat=EDGDSNORM(j,end);disp ('la2');dx_best=dx_bestDSNORM(j,end);
        elseif (rake1a(R)>-45 && rake1a(R)<=45)%Strike slip
            EDGmat=EDGSS(j,end);disp ('la3');dx_best=dx_bestSS(j,end);
        elseif (rake1a(R)>135 && rake1a(R)<=180)%Strike slip
            EDGmat=EDGSS(j,end);disp ('la4');dx_best=dx_bestSS(j,end);
        elseif (rake1a(R)<-135 && rake1a(R)>=-180)%Strike slip
            EDGmat=EDGSS(j,end);disp ('la5');dx_best=dx_bestSS(j,end);
        end
    end
       dxint1(R)=dx_best;
        xmin1(R)=-EDGmat*4;
        xmax1(R)=EDGmat*4;
end
dxint=max(max(dxint1));
xmin=max(max(-EDGmat*4));
xmax=min(min(EDGmat*4));

for R=1:numel(WidthSuma)
    %% Okubo calculation
    disp([num2str(R),'/',num2str(numel(WidthSuma)),' Okubo calculation for Sumatra'])
    %     %% Define Est North-vector for Okubo
    E=xmin:dxint:xmax;
    N=xmin:dxint:xmax;
    %     E=xmin(R):dxint(R):xmax(R);
    %     N=xmin(R):dxint(R):xmax(R);
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
%indlocallon=find(wrapTo360(longT2(1,:))>=min(LO2) & wrapTo360(longT2(1,:))<=max(LO2));
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

% References
% models
% Leonard, M. (2010). Earthquake fault scaling: Self-consistent relating of rupture length, width, average displacement, and moment release. Bulletin of the Seismological Society of America, 100(5A), 1971-1988.
% Wells, D. L., & Coppersmith, K. J. (1994). New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. Bulletin of the seismological Society of America, 84(4), 974-1002.
% Yen, Y. T., & Ma, K. F. (2011). Source-scaling relationship for M 4.6–8.9 earthquakes, specifically for earthquakes in the collision zone of Taiwan. Bulletin of the Seismological Society of America, 101(2), 464-481.
% Gravity change
% Okubo, S. (1992). Gravity and potential changes due to shear and tensile faults in a half?space. Journal of Geophysical Research: Solid Earth, 97(B5), 7137-7144.