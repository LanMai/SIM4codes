clear all
close all
clc

%% Loading system OTF file
OTFo = double(imread('OTF.tif'));
OTFo = OTFpost(OTFo); 

%% Read Expt. Raw SIM Images
aa1 = imread('sim03roi2z4.tif');
aa = aa1(:,:,2);
bb = uint16(zeros(512,512,9));
for ii=1:3
     for jj=1:3
        bb(1:512,1:512,(ii-1)*3+jj)=aa((ii-1)*512+1:ii*512,(jj-1)*512+1:jj*512,1);
     end
end

%% Pre-processing of Raw SIM Images
S1aTnoisy = PreProcessingF( bb(:,:,1) );
S2aTnoisy = PreProcessingF( bb(:,:,2) );
S3aTnoisy = PreProcessingF( bb(:,:,3) );
S1bTnoisy = PreProcessingF( bb(:,:,4) );
S2bTnoisy = PreProcessingF( bb(:,:,5) );
S3bTnoisy = PreProcessingF( bb(:,:,6) );
S1cTnoisy = PreProcessingF( bb(:,:,7) );
S2cTnoisy = PreProcessingF( bb(:,:,8) );
S3cTnoisy = PreProcessingF( bb(:,:,9) );
clear aa aa1 bb

% optional step (may be used if it produces visually better results)
S2aTnoisy = imhistmatch(S2aTnoisy,S1aTnoisy);
S3aTnoisy = imhistmatch(S3aTnoisy,S1aTnoisy);
S1bTnoisy = imhistmatch(S1bTnoisy,S1aTnoisy);
S2bTnoisy = imhistmatch(S2bTnoisy,S1aTnoisy);
S3bTnoisy = imhistmatch(S3bTnoisy,S1aTnoisy);
S1cTnoisy = imhistmatch(S1cTnoisy,S1aTnoisy);
S2cTnoisy = imhistmatch(S2cTnoisy,S1aTnoisy);
S3cTnoisy = imhistmatch(S3cTnoisy,S1aTnoisy);

S1aTnoisy = double( S1aTnoisy ); 
S2aTnoisy = double( S2aTnoisy ); 
S3aTnoisy = double( S3aTnoisy ); 
S1bTnoisy = double( S1bTnoisy ); 
S2bTnoisy = double( S2bTnoisy ); 
S3bTnoisy = double( S3bTnoisy ); 
S1cTnoisy = double( S1cTnoisy ); 
S2cTnoisy = double( S2cTnoisy ); 
S3cTnoisy = double( S3cTnoisy ); 


%% 7-frame
%{
n = 7;
ModFacEst = 0.8.*ones(n,1);
w = 512;
wo = w/2;
Snoisy = zeros(w,w,n);
Snoisy(:,:,1) = S1aTnoisy;
Snoisy(:,:,2) = S2aTnoisy;
Snoisy(:,:,3) = S3aTnoisy;
Snoisy(:,:,4) = S1bTnoisy;
Snoisy(:,:,5) = S2bTnoisy;
Snoisy(:,:,6) = S1cTnoisy;
Snoisy(:,:,7) = S2cTnoisy;
%}

%% 6-frame
%{
n = 6; 
ModFacEst = 1.0.*ones(n,1);
w = 512;
wo = w/2;
Snoisy = zeros(w,w,n);
Snoisy(:,:,1) = S1aTnoisy;
Snoisy(:,:,2) = S2aTnoisy;
Snoisy(:,:,3) = S3aTnoisy;
Snoisy(:,:,4) = S1bTnoisy;
Snoisy(:,:,5) = S2bTnoisy;
Snoisy(:,:,6) = S1cTnoisy;
%}

%% 5-frame
%%{
n = 5;
ModFacEst = 1.0.*ones(n,1);
w = size(OTFo,1);
wo = w/2;
Snoisy = zeros(w,w,n);
Snoisy(:,:,1) = S1aTnoisy;
Snoisy(:,:,2) = S2aTnoisy;
Snoisy(:,:,3) = S3aTnoisy;
Snoisy(:,:,4) = S1bTnoisy;
Snoisy(:,:,5) = S1cTnoisy;
%}

%% 4-frame
%{
n = 4;
ModFacEst = 1.0.*ones(n,1);
w = size(OTFo,1);
wo = w/2;
Snoisy = zeros(w,w,n);
Snoisy(:,:,1) = S1aTnoisy;
Snoisy(:,:,2) = S3aTnoisy;
Snoisy(:,:,3) = S1bTnoisy;
Snoisy(:,:,4) = S1cTnoisy;
%}



clear S1aTnoisy S2aTnoisy S3aTnoisy 
clear S1bTnoisy S2bTnoisy S3bTnoisy
clear S1cTnoisy S2cTnoisy S3cTnoisy

k2a = zeros(n,2);
PhaseA = zeros(n,1);
Spattern = zeros(w,w,n);
PSFe = fspecial('gaussian',14,1.7);
for i = 1:n    
    S1aTnoisy = Snoisy(:,:,i);
    
    [k2a(i,:),PhaseA(i)] = PCMseparateF(S1aTnoisy,OTFo,PSFe);
    
    
    Spattern(:,:,i) = PatternCheckF(Snoisy(:,:,i),k2a(i,:),PhaseA(i),ModFacEst(i));
    
end
% kk



u = 202;
uo = u/2;
OTFo = OTFresize(OTFo,u);
k2a = k2a.*(u/w);
[ MaskPetals, doubleSize ] = MaskPetalsF(OTFo,k2a);

PSFe = fspecial('gaussian',16,2.0); % for edgetapering

xLeft = 50;
yTop = 120;    
[fG1, fG3]  = SIMfreqDeconvAngF(n, ModFacEst, ...
                    OTFo, Snoisy(xLeft+1:xLeft+u,yTop+1:yTop+u,:),...
                    Spattern(xLeft+1:xLeft+u,yTop+1:yTop+u,:), PSFe, MaskPetals, doubleSize);
    %    
      
    %{
    figure;
    surf(log(abs(fG1)),'EdgeColor','none')
    colormap jet
    figure;
    surf(log(abs(fG3)),'EdgeColor','none')
    colormap jet
    %}
    
    %{
    fSig = fG1 - 1*sqrt(abs(fG3));    
    figure;
    surf( log(abs(fSig)).*(log(abs(fSig))>5), 'EdgeColor','none')
    colormap(jet)
    box on
    
    figure;
    hold on
    plot( [1:2*u]-u,log(abs(fG1(u+1,:))),'o-','linewidth',2)
    plot( [1:2*u]-u,0.5*log(abs(fG3(u+1,:))),'o-','linewidth',2)
    plot( [1:2*u]-u,log(abs(fSig(u+1,:))),'o-','linewidth',2)
    grid on
    box on
    
    figure;
    hold on
    plot( [1:2*u]-u,log(abs(fG1(u+1,:))),'o-','linewidth',2)
    plot( [1:2*u]-u,0.5*log(abs(fG3(u+1,:))),'o-','linewidth',2)
    plot( [1:2*u]-u,log(abs(fG1(u+1,:)))-0.5*log(abs(fG3(u+1,:))),'o-','linewidth',2)
    grid on
    box on
    %}
    
%% Determining the object power spectrum
OBJparaA = OBJ4powerPara(fG1,fG3,OTFo, doubleSize);
co = 1.0;
fG1f = W4FilterCenter(fG1,fG3,co,OBJparaA);
G1f = real( ifft2(fftshift(fG1f)) );

v = size(fG1,1);  
h = 20;
figure;
imshow(G1f(h+1:v-h,h+1:v-h),[]) 
title('SR image with artefact')
   
top  = max( max(max(abs(fG1))), max(max(sqrt(abs(fG3)))) );
top  = max( top, max(max(abs(fG1f))) );   
fG1_temp = top*(1-MaskPetals) + abs(fG1);
fG3_temp = top*(1-MaskPetals) + sqrt(abs(fG3));
fG1f_temp = top*(1-MaskPetals) + abs(fG1f);
bottom = min( min(min(fG1_temp)), min(min(sqrt(fG3_temp))) );
bottom = min( bottom, min(min(fG1f_temp)) );
clear fG1_temp fG3_temp fG1f_temp


figure;
surf( log(abs(fG1)), 'EdgeColor','none')
colormap(jet)
axis([0 v 0 v])
box on
caxis manual % This sets the limits of the colorbar to manual for the first plot
caxis([log(bottom) log(top)]);
colorbar;
axis equal
   
figure;
surf( log(sqrt(abs(fG3))) , 'EdgeColor','none')
colormap(jet)
axis([0 v 0 v])
box on
caxis manual % This sets the limits of the colorbar to manual for the first plot
caxis([log(bottom) log(top)]);
colorbar;
axis equal
    
figure;
surf( log(abs(fG1f)) , 'EdgeColor','none')
colormap(jet)
axis([0 v 0 v])
box on
caxis manual % This sets the limits of the colorbar to manual for the first plot
caxis([log(bottom) log(top)]);
colorbar;
axis equal   

%% suppressing the spurious frequency-peaks at illumination frequencies
[fG1n, NotchMask] = PeakNotchFilterF(fG1f,k2a);
G1n = real( ifft2(fftshift(fG1n)) );

figure;
imshow(G1n(h+1:v-h,h+1:v-h),[])
title('SR image with artefact suppressed')

figure;
surf( log(abs(fG1n)) , 'EdgeColor','none')
colormap(jet)
axis([0 v 0 v])
box on
caxis manual % This sets the limits of the colorbar to manual for the first plot
caxis([log(bottom) log(top)]);
colorbar;
axis equal 
    

%% Wiener Filtering of wide-field image
OTFo = double(imread('OTF.tif'));
OTFo = OTFpost(OTFo); 
PSFe = fspecial('gaussian',16,2.0);
aa1 = imread('widefield03roi2z4.tif'); % wide-field image
DIoTnoisy = uint16( aa1(:,:,2) );
clear aa1
DIoTnoisy = PreProcessingF( DIoTnoisy );
Dwf = WideFieldDeconvF(DIoTnoisy,OTFo,PSFe);

Dwf0 = Dwf(xLeft+1:xLeft+u,yTop+1:yTop+u);
DIoTnoisy0 = DIoTnoisy(xLeft+1:xLeft+u,yTop+1:yTop+u);

h = 10;
figure;
imshow(DIoTnoisy0(h+1:u-h,h+1:u-h),[]);
title('Wide-field image')
figure;
imshow(Dwf0(h+1:u-h,h+1:u-h),[]);
title('Wiener Filtered wide-field image')

