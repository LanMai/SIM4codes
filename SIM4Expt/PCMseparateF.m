function [kA,phaseA] = PCMseparateF(S1aTnoisy,OTFo,PSFe)
% AIM: to estimate illumination frequency and phase
% INPUT VARIABLES
%   S1aTnoisy: raw SIM image
%   OTFo: system OTF
%   PSFe: for edgetapering
% OUTPUT VARIABLES
%    kA: illumination frequency
%    phaseA: illumination phase


w = size(S1aTnoisy,1);
wo = w/2;

%% Determination of illumination frequency vectors
kA = IlluminationFreqF(S1aTnoisy,OTFo,PSFe);
magkA = sqrt(kA*kA')

%% determination of illumination phase shifts
[phaseA] = IlluminationPhaseF(S1aTnoisy,kA);
phaseA*180/pi

kA = [kA(2) kA(1)];