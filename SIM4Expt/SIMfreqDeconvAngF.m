function [fG1, fG2] = SIMfreqDeconvAngF( n, ModFac, ...
                    OTFo, Snoisy, Spattern, PSFe, MaskPetals, doubleSize)
%

if ( doubleSize > 0 )
    OTFo = DoubleMatSizeF(OTFo);
end

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Co = (X-wo) + 1i*(Y-wo);
Ro = abs(Co);

Kotf = OTFedgeF(OTFo);
Zo = 1*(Ro<Kotf);

tic
Tbegin = toc;
for i = 1:n
    [ Mtemp{i}, Ftemp{i}, NoiseStd{i} ] = SIMangMatrixF( ModFac(i), Kotf, MaskPetals, ...
        Snoisy(:,:,i), Spattern(:,:,i), PSFe, OTFo, doubleSize ); 
    [100*i toc]
end
Tend = toc;
dT = Tend - Tbegin
% kk

M = Mtemp{1};
F = Ftemp{1};
for i = 2:n
   M = cat(1, M, Mtemp{i} );
   F = cat(1, F, Ftemp{i} );
end
clear Mtemp Ftemp

Mbegin = toc;
G1 = (M'*M)\(M'*F);
Mend = toc;
Mdt = Mend - Mbegin

Zcount = sum(sum(Zo));
IndexZo = zeros(Zcount,1);
InCounter = 0;
for u = 1:w^2
    if ( Zo(u)>0 )
        InCounter = InCounter + 1;
        IndexZo(InCounter,:) =  u;
    end
end

NoiseIterations = 512;
NoiseTemplate = zeros(n*Zcount,NoiseIterations);
for j = 1:NoiseIterations
    for i = 1:n
        ImageNoise = fftshift(fft2( random('norm', 0, NoiseStd{i}, w , w) ));
        NoiseTemplate((i-1)*Zcount+1:i*Zcount, j ) = ImageNoise(IndexZo);
    end
end

Mbegin = toc;
G2temp = (M'*M)\(M'*NoiseTemplate);
G2temp = G2temp.*conj(G2temp);
G2 = sum(G2temp,2)./NoiseIterations;
Mend = toc;
Mdt = Mend - Mbegin 

fG1 = zeros(w,w);
fG2 = zeros(w,w);

countC = 0;
for t = 1:w^2
	if (MaskPetals(t)>0)
        j = mod(t-1,w) + 1; % column
        i = (t-j)/w + 1; % row
        countC = countC + 1; 
        fG1(j,i) = G1(countC);
        fG2(j,i) = G2(countC);
    end
end 
