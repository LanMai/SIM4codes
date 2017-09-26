function [ S1aTnoisy ] = PreProcessingF( S1aTnoisy )

%{
S1aTnoisy0 = adapthisteq(S1aTnoisy);
figure;
imshow(S1aTnoisy0,[])
figure;
mesh(S1aTnoisy0)
%}

S1aTnoisy = double( S1aTnoisy ); 
minS = min(S1aTnoisy(:));
maxS = max(S1aTnoisy(:));
S1aTnoisy = 0.6*(S1aTnoisy - minS)./(maxS - minS);
S1aTnoisy = im2uint16(S1aTnoisy);
S1aTnoisy = adapthisteq(S1aTnoisy); % CLAHE
%{
figure;
mesh(S1aTnoisy)
figure;
imshow(S1aTnoisy,[])
%}


