clear
close all
clc
%% exemple de grande grille
%grille seisme
%% Input model parameters
dx_best=importdata('dx_bestDSREV.mat');EDGmat=importdata('EDGmatDSREV.mat')
Mwint=0:0.25:9 ;%used later for the interpolation of dx depending on the Magnitude
%   xF = -54197:368.8190:54197;
%   yF = -54197:368.8190:54197;
  Mwvec=0:0.2:9;
  Model=5;
  for R=1:numel(Mwvec)
      for j=Model
          xint(R)=interp1(Mwint,EDGmat(j,:),Mwvec(R),'linear');
          dxint(R)=interp1(Mwint,dx_best(j,:),Mwvec(R),'linear');
          xmin(R)=-round(xint(R));
          xmax(R)=round(xint(R));
      end
  xF = xmin(R):dxint(R):xmax(R);
  yF = xmin(R):dxint(R):xmax(R);
%% faux epicentre
los = 3.80%148.64; % deg deci
las = 78.79%-62,99;
%% translation
[x0,y0,ZONE] = ll2utm(las,los);
xG = xF+x0;
yG = yF+y0;
%% Synth Z (sim Okubo)
[XG,YG] = meshgrid(xG,yG);
  Z = (-50*XG.^2+3*YG.^2)*1e-13; % random
%% exemple de grille d'interpolation
res = 0.05;
[LAB, LOB] = meshgrid(las-40:res:las+40,los-40:res:los+40);
%% OPTION #1
[la{R},lo{R}] = utm2ll(xG,yG,ZONE);
[LO1{R},LA1{R}] = meshgrid(lo{R},la{R});
%% OPTION #2
[LA2{R},LO2{R}] = utm2ll(XG(1:numel(XG)),YG(1:numel(YG)),ZONE);
  end

  for R=1:numel(Mwvec)
%% visu
%m_proj('Miller Cylindrical','latitudes',[min(la{R}) max(la{R})] ,'longitudes',[min(lo{R}) max(lo{R})],1);
m_proj('Miller Cylindrical','latitudes',[-70 -50] ,'longitudes',[125 165],1);
% fig points
set(figure(R),'Position',[1 1 900 900])
set(gcf,'PaperPositionMode','auto')
m_plot(LO1{R},LA1{R},'k.')
hold on
m_plot(LO2{R},LA2{R},'r.')
m_coast('color','k','LineWidth',1);
m_grid
title(['Compare positions points, Magnitude',num2str(R)])
  end

  for R=1:numel(Mwvec)
deltaLO2{R}=max(LO2{R})-min(LO2{R});
deltaLO1{R}=max(max(LO1{R}))-min(min(LO1{R}));
ratioLO(R)=deltaLO1{R}*100/deltaLO2{R}

  end

  figure(111);
  plot(Mwvec(:),ratioLO(:))
% % fig tot
% set(figure,'Position',[1 1 1900 920])
% set(gcf,'PaperPositionMode','auto')
% subplot(2,3,1)
% m_plot(LO1,LA1,'k.')
% hold on
% m_plot(LOB,LAB,'r.')
% m_plot(LO2,LA2,'r.')
% m_gshhs_l('color','k','LineWidth',1);
% m_grid
% title('Methode meshgrid')
%
% subplot(2,3,4)
% m_pcolor(LOB,LAB,Zi1)
% hold on
% m_plot(LOB,LAB,'r.')
% m_plot(LO2,LA2,'r.')
% m_gshhs_l('color','k','LineWidth',1);
% m_grid
% colorbar
%
% subplot(2,3,2)
% m_plot(LO2,LA2,'.')
% hold on
% m_plot(LOB,LAB,'r.')
% m_gshhs_l('color','k','LineWidth',1);
% m_grid
% title('Methode directe')
%
% subplot(2,3,5)
% m_pcolor(LOB,LAB,Zi2)
% hold on
% m_plot(LOB,LAB,'r.')
% m_gshhs_l('color','k','LineWidth',1);
% m_grid
% colorbar
%
% subplot(2,3,6)
% m_pcolor(LOB,LAB,Zi2-Zi1)
% hold on
% m_plot(LOB,LAB,'r.')
% m_gshhs_l('color','k','LineWidth',1);
% m_grid
% colorbar
% title('Methode directe - Methode meshgrid')
