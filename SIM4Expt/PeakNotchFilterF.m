function [fG1n,NotchMask] = PeakNotchFilterF(fG1,k2a)

ko = zeros(3,2);
ko(1,:) = k2a(1,:);

n = size(k2a,1);
count = 1;
for i = 2:n
    Cv = k2a(i,1) + 1i.*k2a(i,2);
    Co = ko(count,1) + 1i.*ko(count,2);
    dAng = 180/pi*angle( Cv*conj(Co) );
    % [i dAng]
    if ( abs(dAng) > 5 )
        count = count + 1;
        ko(count,:) = k2a(i,:);
    end
end

w = size(fG1,1);
wo = w/2;
x = linspace(0,w-1,w);
y = linspace(0,w-1,w);
[X,Y] = meshgrid(x,y);
Co = (X-wo) + 1i*(Y-wo);
Ro = abs(Co);

alpha = 0.07;
beta = 1.0;
So = 1 - exp(-alpha*(Ro.^beta));
%{
figure;
mesh(So)
%}

%{
figure;
hold on
plot( x-wo,So(wo+1,:),'o-','linewidth',2)
%}

NotchMask = ones(w,w);
for i = 1:3
    Cv = ko(i,1) + 1i.*ko(i,2);
    Rp = abs(Co - Cv);
    Rm = abs(Co + Cv);
    Sp = 1 - exp(-alpha*(Rp.^beta));
    Sm = 1 - exp(-alpha*(Rm.^beta));
    NotchMask = NotchMask.*Sp.*Sm;
end
%{
figure;
mesh(NotchMask)

Znotch = 1.*(NotchMask>0.97);
figure;
mesh(Znotch)
%}

fG1n = fG1.*NotchMask;

        
    