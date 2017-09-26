function [kA,phaseA] = PCMseparateF(S1aTnoisy,OTFo,PSFe)
% AIM: obtaining the noisy estimates of three frequency components
% INPUT VARIABLES
%   S1aTnoisy,S2aTnoisy,S3aTnoisy: 3 raw SIM images with identical 
%                                   illumination pattern orientation 
%                                   but different phase shifts
%   OTFo: system OTF
% OUTPUT VARIABLES
%   fDo,fDp,fDm: noisy estimates of separated frequency components
%   kA: (averaged) illumination frequency vector


w = size(S1aTnoisy,1);
wo = w/2;

%% Determination of illumination frequency vectors
kA = IlluminationFreqF(S1aTnoisy,OTFo,PSFe);
magkA = sqrt(kA*kA')

%% determination of illumination phase shifts
[phaseA] = IlluminationPhaseF(S1aTnoisy,kA);
phaseA*180/pi

kA = [kA(2) kA(1)];