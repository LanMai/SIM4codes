function sAo = PatternCheckF(S1aTnoisy,kA,phaseA,ModFacEst)

w = size(S1aTnoisy,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);

k2a = kA./w;
maxS = 0.6*max(S1aTnoisy(:));
sAo = 0.5 + 0.5*ModFacEst*cos(2*pi*(k2a(1,1).*(X-wo)+k2a(1,2).*(Y-wo))+phaseA);
Smix = S1aTnoisy;
pmix = 8;
pp = w/pmix;
for u = 1:pmix
    for v = 1:pmix
        if ( mod(u+v,2)==0 )
            Smix((u-1)*pp+1:u*pp,(v-1)*pp+1:v*pp) = maxS.*sAo((u-1)*pp+1:u*pp,(v-1)*pp+1:v*pp);
        end
    end
end
%% for visual verification
% figure;
% imshow(Smix,[])
% imshow(Smix(1:wo,1:wo),[])