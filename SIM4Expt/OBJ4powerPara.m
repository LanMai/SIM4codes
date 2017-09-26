function OBJparaA = OBJ4powerPara(fG1,fG3,OTFo,doubleSize)
% AIM: determination of object power parameters Aobj and Bobj
% INPUT VARIABLES
%   fDIoTnoisy: FT of central frequency component
%   OTFo: system OTFo
% OUTPUT VARIABLES
%   OBJparaA: [Aobj Bobj]

if ( doubleSize > 0 )
    OTFo = DoubleMatSizeF(OTFo);
end

w = size(OTFo,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Cv = (X-wo) + 1i*(Y-wo);
Ro = abs(Cv);

%% object power parameters through optimization
OBJparaOpt0 = @(OBJpara0)OBJ4paraOpt(OBJpara0,fG1,fG3,OTFo);
options = optimset('LargeScale','off','Algorithm',...
	'active-set','MaxFunEvals',500,'MaxIter',500,'Display','notify');

% obtaining crude initial guesses for Aobj and Bobj 
Kotf = OTFedgeF(OTFo);
Zm = (Ro>0.3*Kotf).*(Ro<0.4*Kotf);
Aobj = sum(sum(abs(fG1.*Zm)))./sum(sum(Zm));
Bobj = -0.5;
OBJpara0 = [Aobj Bobj];

% optimization step
[OBJparaA,fval] = fminsearch(OBJparaOpt0,OBJpara0,options);
