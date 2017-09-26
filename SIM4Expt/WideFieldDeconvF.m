function Dwf = WideFieldDeconvF(DIoTnoisy,OTFo,PSFe)

%% Object power parameters determination
OBJparaA = OBJpowerPara(DIoTnoisy,OTFo,PSFe);
DIoTnoisy = edgetaper(DIoTnoisy,PSFe);
fDIoTnoisy = fftshift(fft2(DIoTnoisy));

%% Wiener Filtering wide-field image
% (do not edge-taper DIoTnoisy prior to filtering) 
SFo = 1;
co = 1.0; 
[fDIoT,npDo,WFilter] = WoFilterCenter(fDIoTnoisy,OTFo,co,OBJparaA,SFo);
Dwf = real( ifft2(fftshift(fDIoT)) );
%{
figure;
imshow(Dwf,[])
figure;
imshow(DIoTnoisy,[])
%}