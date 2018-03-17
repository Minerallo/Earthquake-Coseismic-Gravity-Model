clear
close all
clc
%% exemple de grande grille
  xF = -54197:368.8190:54197;
  yF = -54197:368.8190:54197;
%% faux epicentre
los = 152.78; % deg deci
las = -60,76;
%% translation
[x0,y0,ZONE] = ll2utm(las,los);
xG = xF+x0;
yG = yF+y0;
%% Synth Z (sim Okubo)
[XG,YG] = meshgrid(xG,yG);
Z = (-50*XG.^2+3*YG.^2)*1e-13; % random
%% exemple de grille d'interpolation
res = 0.05;
[LAB, LOB] = meshgrid(las-5:res:las+5,los-5:res:los+5);
%% OPTION #1 Passage Lat/Lon avec meshgrid
[la,lo] = utm2ll(xG,yG,ZONE);
% [la,lo] = utm2ll(xG,yG,ZONE);
[LO1,LA1] = meshgrid(lo,la); 
% interp 
Zi1 = interp2(LO1,LA1,Z,LOB,LAB,'linear');
%% OPTION # 2 Passage Lat/Lon direct depuis XG/YG
[LA2,LO2] = utm2ll(XG(1:numel(XG)),YG(1:numel(YG)),ZONE);
tic
Gq=gridfit(LO2,LA2,Z(1:numel(Z)),LOB(:,1),LAB(1,:),'tilesize',150,'overlap',0.25);
Zi2 = Gq';
toc
        
%Zi2 = biharmonic_spline_interp2(LO2',LA2',Z(1:numel(Z))',LOB,LAB);
% Zi2 = griddata(LO2,LA2,Z(1:numel(Z)),LOB,LAB);

%% visu
m_proj('Miller Cylindrical','latitudes',[-5 5]+las,'longitudes',[-5 5]+los,1);
% fig points
set(figure,'Position',[1 1 900 900])
set(gcf,'PaperPositionMode','auto')
m_plot(LO1,LA1,'k.')
hold on
m_plot(LO2,LA2,'r.')
m_gshhs_l('color','k','LineWidth',1);
m_grid
title('Compare positions points')

% fig tot
set(figure,'Position',[1 1 1900 920])
set(gcf,'PaperPositionMode','auto')
subplot(2,3,1)
m_plot(LO1,LA1,'k.')
hold on
m_plot(LOB,LAB,'r.')
m_plot(LO2,LA2,'r.')
m_gshhs_l('color','k','LineWidth',1);
m_grid
title('Methode meshgrid')

subplot(2,3,4)
m_pcolor(LOB,LAB,Zi1)
hold on
m_plot(LOB,LAB,'r.')
m_plot(LO2,LA2,'r.')
m_gshhs_l('color','k','LineWidth',1);
m_grid
colorbar

subplot(2,3,2)
m_plot(LO2,LA2,'.')
hold on
m_plot(LOB,LAB,'r.')
m_gshhs_l('color','k','LineWidth',1);
m_grid
title('Methode directe')

subplot(2,3,5)
m_pcolor(LOB,LAB,Zi2)
hold on
m_plot(LOB,LAB,'r.')
m_gshhs_l('color','k','LineWidth',1);
m_grid
colorbar

subplot(2,3,6)
m_pcolor(LOB,LAB,Zi2-Zi1)
hold on
m_plot(LOB,LAB,'r.')
m_gshhs_l('color','k','LineWidth',1);
m_grid
colorbar
title('Methode directe - Methode meshgrid')
