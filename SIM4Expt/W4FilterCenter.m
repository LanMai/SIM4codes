function [fG1f] = W4FilterCenter(fG1,NoisePower,co,OBJpara)
% Aim: Wiener Filtering the central frequency component
% INPUT VARIABLES
%   FiSMao: noisy central frequency component
%   OTFo: system OTF
%   co: Wiener filter constant [=1, for minimum RMS estimate]
%   SFo: scaling factor (not significant here, so set to 1)
% OUTPUT VARIABLES
%   FiSMaof: Wiener Filtered estimate of FiSMao
%   NoisePower: avg. noise power in FiSMao

w = size(fG1,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Ro = sqrt( (X-wo).^2 + (Y-wo).^2 );

% Object Power determination
Aobj = OBJpara(1);
Bobj = OBJpara(2);
Ro(wo+1,wo+1) = 1;
OBJpower = Aobj*(Ro.^Bobj);
OBJpower = OBJpower.^2;

%% Wiener Filtering
fG1f = fG1./(1 + co.*NoisePower./OBJpower);
WFilter = 1./(1 + co.*NoisePower./OBJpower);

%% for cross-checking filtered estimate visually
%{
FiSMao1 = fG1f;
pp = 3;
figure;
hold on
plot(x-wo,log(abs(fG1f(wo+pp,:))),'go--','LineWidth',2,'MarkerSize',6)
plot(x-wo,log(abs(fG1(wo+pp,:))),'o--','LineWidth',2,'MarkerSize',6)
plot(x-wo,log(abs(fG1f(wo+pp,:))),'r*--','LineWidth',2,'MarkerSize',4)
plot(x-wo,10*WFilter(wo+pp,:),'k*--','LineWidth',2,'MarkerSize',4)
title('WoFilter')
legend('FiSMaof','FiSMao','FiSMao1')
grid on
box on
%kk
%}