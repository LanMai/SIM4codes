function OTFtemp = DoubleMatSizeF(OTFo)

w = size(OTFo,1);
wo = w/2;

t = 2*w;
to = t/2;
OTFtemp = zeros(t,t);
OTFtemp(to-wo+1:to+wo,to-wo+1:to+wo) = OTFo;



