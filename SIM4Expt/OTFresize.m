function OTFo = OTFresize(OTFtemp,t)

w = size(OTFtemp,1);
wo = w/2;
to = t/2;
PSFtemp = real( fftshift( ifft2(fftshift(OTFtemp)) ) );
PSFo = PSFtemp(wo-to+1:wo+to,wo-to+1:wo+to);

OTFo = abs( fftshift(fft2(PSFo)) );
%{
figure;
mesh(OTFo)
figure;
mesh(OTFtemp)
% kk
%}