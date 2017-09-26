function [k2fa] = IlluminationFreqF(S1aTnoisy,OTFo,PSFe)
% AIM: illumination frequency vector determination
% INPUT VARIABLES
%   S1aTnoisy: raw SIM image
%   OTFo: system OTF
% OUTPUT VARIABLE
%   k2fa: illumination frequency vector

w = size(OTFo,1);
wo = w/2;

% edge tapering raw SIM image
S1aTnoisy_et = edgetaper(S1aTnoisy,PSFe);
fS1aTnoisy_et = fftshift(fft2(S1aTnoisy_et));

% OTF cut-off freq
Kotf = OTFedgeF(OTFo);

% Approx illumination frequency vector
[k2fa,~,~] = ApproxFreqDuplex(fS1aTnoisy_et,Kotf);

fS1aTnoisy = fftshift(fft2(S1aTnoisy));
% illumination frequency vector determination by optimizing
% autocorrelation of fS1aTnoisy
OPT = 1;
PhaseKai2opt0 = @(k2fa0)PhaseKai2opt(k2fa0,fS1aTnoisy,OTFo,OPT);
options = optimset('LargeScale','off','Algorithm',...
	'active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');
% options = optimset(options,'UseParallel','always');
k2fa0 = k2fa;
% tic;
[k2fa,fval] = fminsearch(PhaseKai2opt0,k2fa0,options);
% toc
% k2a = sqrt(k2fa*k2fa')
