function [ MaskPetals, doubleSize ] = MaskPetalsF(OTFo,k2a)

w = size(OTFo,1);
wo = w/2;

Kotf = OTFedgeF(OTFo);
kTemp = sqrt( k2a(1,:)*k2a(1,:)' ) + Kotf;
doubleSize = 0;
if ( kTemp+20 > wo )
    doubleSize = 1;
    w = 2*w;
    wo = w/2;
end

x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Co = (X-wo) + 1i*(Y-wo);
Ro = abs(Co);

Zo = 1*(Ro<Kotf);
MaskPetals = Zo;
n = size(k2a,1);
for i = 1:n
    ka = round(k2a(i,:));
    Ca = ka(1) + 1i*ka(2);
    Rap = abs(Co-Ca);
    Ram = abs(Co+Ca);
    Zap = 1*(Rap<Kotf);
    Zam = 1*(Ram<Kotf);
    MaskPetals = 1 - (1-MaskPetals).*(1-Zap).*(1-Zam);
end
%{
figure;
mesh(Zo)
figure;
mesh(MaskPetals)
[ Amask ] = AppoMask( MaskPetals, 0.4 );
figure;
mesh(Amask)
kk
%}