function [ MAp, fSp, NoiseStd ] = SIMangMatrixF( ModFac, Kotf, MaskPetals, ...
    S1aTnoisy, sAo, PSFe, OTFo, doubleSize )

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Co = (X-wo) + 1i*(Y-wo);
Ro = abs(Co);
Zo = 1*(Ro<Kotf);

sAo = edgetaper(sAo,PSFe);
fsAo = fftshift(fft2(sAo));
if ( doubleSize > 0 )
    fsAo = DoubleMatSizeF(fsAo);
end
aA = 0.5*ModFac;
Threshold = 1*0.005*(aA/2);
MaskPeak = 1.*(abs(fsAo)./w^2>Threshold);
fsAo1 = fsAo.*MaskPeak;
fsAoE = zeros(w+2*Kotf,w+2*Kotf);
fsAoE(Kotf+1:Kotf+w,Kotf+1:Kotf+w) = fsAo1;

S1aTnoisy = edgetaper(S1aTnoisy,PSFe);
fS1aTnoisy = fftshift( fft2(S1aTnoisy) );
if ( doubleSize > 0 )
    fS1aTnoisy = DoubleMatSizeF(fS1aTnoisy);
end

Znoise = 1*(Ro>(Kotf+10));
hh = 10;
Znoise1 = (X>hh).*(X<(w-hh)).*(Y>hh).*(Y<(w-hh));
Znoise = Znoise.*Znoise1;
Noise0 = fS1aTnoisy.*Znoise;
Noise = sum(sum( abs(Noise0.*conj(Noise0)) ))./sum(sum(Znoise));
NoiseStd = sqrt(Noise./w^2);

Ncount = sum(sum(MaskPetals));
Zcount = sum(sum(Zo));
MAp = sparse(Zcount,Ncount);

IndexZo = zeros(Zcount,2);
IndexMaskPetals = zeros(Ncount,1);
InCounter = 0;
OutCounter = 0;
for u = 1:w^2
    
    if ( Zo(u)>0 )
        InCounter = InCounter + 1;
        IndexZo(InCounter,:) = [ u, OTFo(u) ];
    end
    
    if ( MaskPetals(u)>0 )
        OutCounter = OutCounter + 1;
        IndexMaskPetals(OutCounter) = u;
    end
    
end

fSp = fS1aTnoisy(IndexZo(:,1));

parfor u = 1:Zcount
    v = IndexZo(u,1);
    j = mod(v-1,w) + 1; % column
    i = (v-j)/w + 1; % row
    sy = i - (wo+1);
    sx = j - (wo+1);
    fsAo0 = circshift(fsAoE,[sx sy]);
    fsAo1 = fsAo0(Kotf+1:Kotf+w,Kotf+1:Kotf+w);
    
    MAp(u,:) = IndexZo(u,2).*fsAo1(IndexMaskPetals)./w^2; 
end

MAp = conj(MAp);



