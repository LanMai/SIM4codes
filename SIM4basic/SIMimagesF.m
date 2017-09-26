function [Snoisy, DIoTnoisy, DIoT, Spattern] = SIMimagesF(k2,...
          DIo,PSFo,OTFo,ModFac,NoiseLevel,UsePSF,PSFe,AnglePhase,alpha)

% AIM: to generate raw sim images
% INPUT VARIABLES
%   k2: illumination frequency
%   DIo: specimen image
%   PSFo: system PSF
%   OTFo: system OTF
%   UsePSF: 1 (to blur SIM images by convloving with PSF)
%           0 (to blur SIM images by truncating its fourier content beyond OTF)
%   NoiseLevel: percentage noise level for generating gaussian noise
% OUTPUT VARIABLES
%   Snoisy: raw sim images
%   DIoTnoisy: noisy wide field image
%   DIoT: noise-free wide field image
% Spattern: illumination patterns

w = size(DIo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);

mA = 0.5; % mean illumination intensity
aA = 0.5*ModFac; % amplitude of illumination intensity above mean
aNoise = NoiseLevel/100; % corresponds to 10% noise
n = size(AnglePhase,1);
Snoisy = zeros(w,w,n);
Spattern = zeros(w,w,n);
for i = 1:n
    thetaA = AnglePhase(i,1) + alpha;
    % illumination frequency vectors
    k2a = (k2/w).*[cos(thetaA) sin(thetaA)]; 
    
    % illunination phase shifts with random errors
    psAo = AnglePhase(i,2) + 0*(0.5-rand(1,1))*pi;
    
    % illunination patterns
    sAo = mA + aA*cos(2*pi*(k2a(1,1).*(X-wo)+k2a(1,2).*(Y-wo))+psAo);   
    Spattern(:,:,i) = sAo;        
        
    % superposed Object        
    s1a = DIo.*sAo;
    
    % superposed (noise-free) Images        
    PSFsum = sum(sum(PSFo));        
    if ( UsePSF == 1 )            
        S1aT = conv2(s1a,PSFo,'same')./PSFsum;        
    else        
        s1a = edgetaper(s1a,PSFe);            
        S1aT = ifft2( fft2(s1a).*fftshift(OTFo) );            
        S1aT = real(S1aT);        
    end       
    
    % Gaussian noise generation        
    % SNR = 1/aNoise            
    % SNRdb = 20*log10(1/aNoise)        
    nS1aT = random('norm', 0, aNoise*std2(S1aT), w , w);        
    NoiseFrac = 1; %may be set to 0 to avoid noise addition         
    S1aTnoisy = S1aT + NoiseFrac*nS1aT;            
    Snoisy(:,:,i) = S1aTnoisy;
end

if ( UsePSF == 1 )
	DIoT = conv2(DIo,PSFo,'same')./PSFsum;
else
	DIo = edgetaper(DIo,PSFe);
	DIoT = ifft2( fft2(DIo).*fftshift(OTFo) );
	DIoT = real(DIoT);
end
nDIoT = random('norm', 0, aNoise*std2(DIoT), w , w);
DIoTnoisy = DIoT + NoiseFrac*nDIoT;
