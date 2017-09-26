clear all
close all
clc

w = 512;
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);

%% Generation of the PSF with Besselj.
scale = 0.63; % used to adjust PSF/OTF width
[PSFo,OTFo] = PsfOtf(w,scale); 


PSFe = fspecial('gaussian',16,2);  % for edgetapering


%% Reading input file
Io1 = imread('FractalNeuralNetwork.jpg');
Io = Io1(1:w,1:w);
DIo = double(Io);


%% Generating raw SIM Images
alpha = 1*pi/12;
k2 = 80; % illumination freq
ModFac = 0.8; % modulation factor
NoiseLevel = 10.0; % in percentage
UsePSF = 0; % 1(to blur using PSF), 0(to blur using OTF)

%% 4-frame
%%{
p1 = [0 0]; % [angularOrientation phaseShift]
p2 = [0*pi/3 3*pi/3];
p3 = [1*pi/3 0*pi/3];
p4 = [2*pi/3 0];
AnglePhase = [p1; p2; p3; p4]; 
%}

%% 5-frame
%{
p1 = [0 0]; % [angularOrientation phaseShift]
p2 = [0*pi/3 2*pi/3];
p3 = [0*pi/3 4*pi/3];
p4 = [1*pi/3 0];
p5 = [2*pi/3 0];
AnglePhase = [p1; p2; p3; p4; p5]; 
%}

%% 6-frame
%{
p1 = [0 0]; % [angularOrientation phaseShift]
p2 = [0*pi/3 2*pi/3];
p3 = [0*pi/3 4*pi/3];
p4 = [1*pi/3 0];
p5 = [1*pi/3 2*pi/3];
p6 = [2*pi/3 0];
AnglePhase = [p1; p2; p3; p4; p5; p6]; 
%}

%% 7-frame
%{
p1 = [0 0]; % [angularOrientation phaseShift]
p2 = [0*pi/3 2*pi/3];
p3 = [0*pi/3 4*pi/3];
p4 = [1*pi/3 0];
p5 = [1*pi/3 2*pi/3];
p6 = [2*pi/3 0];
p7 = [2*pi/3 2*pi/3];
AnglePhase = [p1; p2; p3; p4; p5; p6; p7]; 
%}

[Snoisy, DIoTnoisy, DIoT,Spattern0] = SIMimagesF(k2,DIo, ...
 PSFo,OTFo,ModFac,NoiseLevel,UsePSF,PSFe,AnglePhase,alpha);

 
%% Wiener Filtering of wide-field image
Dwf = WideFieldDeconvF(DIoTnoisy,OTFo,PSFe);

%%
ModFacEst = ModFac; %0.8;
n = size(AnglePhase,1);
k2a = zeros(n,2);
PhaseA = zeros(n,1);
Spattern = zeros(w,w,n);
for i = 1:n    
    S1aTnoisy = Snoisy(:,:,i);
    
    [k2a(i,:),PhaseA(i)] = PCMseparateF(S1aTnoisy,OTFo,PSFe);
    
    Spattern(:,:,i) = PatternCheckF(Spattern0(:,:,i),k2a(i,:),PhaseA(i),ModFacEst);
    
end


u = 210; % selecting width of the sub-image
uo = u/2;
OTFo = OTFresize(OTFo,u);
k2a = k2a.*(u/w);
[ MaskPetals, doubleSize ] = MaskPetalsF(OTFo,k2a);

if ( u > 180 )
    PSFe = fspecial('gaussian',16,2.0);
else    
    PSFe = fspecial('gaussian',7,0.7);
end

% coordinates for selecting the region from raw SIM images
xLeft = wo; 
yTop = wo;

%% obtaining the least square solution
[fG1, fG3]  = SIMfreqDeconvAngF(n,ModFacEst, ...
                    OTFo, Snoisy(xLeft+1:xLeft+u,yTop+1:yTop+u,:),...
                    Spattern(xLeft+1:xLeft+u,yTop+1:yTop+u,:), PSFe,...
                    MaskPetals, doubleSize, DIo(xLeft+1:xLeft+u,yTop+1:yTop+u));

    
OBJparaA = OBJ4powerPara(fG1,fG3,OTFo);
co = 1.0;
[fG1f,WFilter] = W4FilterCenter(fG1,fG3,co,OBJparaA);
G1f = real( ifft2(fftshift(fG1f)) );
   
h = 20;
figure;
imshow(G1f(h+1:u-h,h+1:u-h),[])
str = int2str(NoiseLevel);
title(['G1f, NoiseLevel = ' int2str(NoiseLevel) ]);   

bottom = min( min(min(abs(fG1))), min(min(sqrt(abs(fG3)))) );
bottom = min( bottom, min(min(abs(fG1f))) );
bottom = max( bottom, 1 );
top  = max( max(max(abs(fG1))), max(max(sqrt(abs(fG3)))) );
top  = max( top, max(max(abs(fG1f))) );
figure;
surf( log(abs(fG1)).*(log(abs(fG1))>5), 'EdgeColor','none')
colormap(jet)
axis([0 u 0 u])
box on
caxis manual % This sets the limits of the colorbar to manual for the first plot
caxis([log(bottom) log(top)]);
colorbar;
axis equal

figure;
surf( log(sqrt(abs(fG3))).*(log(sqrt(abs(fG3)))>5), 'EdgeColor','none')
colormap(jet)
axis([0 u 0 u])
box on
caxis manual % This sets the limits of the colorbar to manual for the first plot
caxis([log(bottom) log(top)]);
colorbar;
axis equal
    
figure;
surf( log(abs(fG1f)).*(log(abs(fG1f))>5), 'EdgeColor','none')
colormap(jet)
axis([0 u 0 u])
box on
caxis manual % This sets the limits of the colorbar to manual for the first plot
caxis([log(bottom) log(top)]);
colorbar;
axis equal    
   
Dwf0 = Dwf(xLeft+1:xLeft+u,yTop+1:yTop+u);
DIo0 = DIo(xLeft+1:xLeft+u,yTop+1:yTop+u);
h = 20;
figure;
imshow(Dwf0(h+1:u-h,h+1:u-h),[])
str = int2str(NoiseLevel);
title(['Dwf, NoiseLevel = ' int2str(NoiseLevel) ]);
figure;
imshow(DIo0(h+1:u-h,h+1:u-h),[])
title('DIo');

Dnoisy = DIoTnoisy(xLeft+1:xLeft+u,yTop+1:yTop+u);
figure;
imshow(Dnoisy(h+1:u-h,h+1:u-h),[])
str = int2str(NoiseLevel);
title(['DIoTnoisy, NoiseLevel = ' int2str(NoiseLevel) ]);
   
