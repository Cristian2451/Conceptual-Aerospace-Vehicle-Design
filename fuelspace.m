clear
clc

% parameters to edit
fspar = 0.15;
rspar = 0.6;
s = 14.1;
Croot = 4.027;
Ctip = 1.611;
sections = 1000;
FusWidth = 2.786;


% import aerofoil coordinates
coords = readmatrix('aerofoil.txt');
coords(length(coords)+1,1) = coords(1,1);
coords(length(coords)+1,2) = coords(1,2);

% increase density of points
for i = 1:3
    coords = interp(coords);
end

% initialise index of new array
j = 1;
fuelbox(1,2) = 0;

for i = 1:length(coords)
    if (coords(i,1)>fspar)&&(coords(i,1)<rspar)
        fuelbox(j,1) = coords(i,1);
        fuelbox(j,2) = coords(i,2);
        j = j + 1;
    end
end

% complete the circle
fuelbox(j,1) = fuelbox(1,1);
fuelbox(j,2) = fuelbox(1,2);

% find area of fuelbox
A = 0;
x = 0.5;
y = 0;
for i = 1:(length(fuelbox)-1)
    A = A + 0.5*abs( x*(fuelbox(i,2)-fuelbox(i+1,2)) + fuelbox(i,1)*(fuelbox(i+1,2)-y) + fuelbox(i+1,1)*(y-fuelbox(i,2)) );
end

% A = (1/2) |x1(y2 − y3) + x2(y3 − y1) + x3(y1 − y2)|

figure
plot(coords(:,1),coords(:,2),'x')
hold on
plot(fuelbox(:,1),fuelbox(:,2),'-')
axis equal
hold off

% initialise volume
Vtot = 0;
Vsec(sections) = 0;
seff = s - 1;
ds = seff/sections;

for i = 1:(sections-1)
    Vsec(i) = ds*A*(chorddist( (i+1)*ds, s, Croot, Ctip ))^2;
    Vtot = Vtot + Vsec(i);
end

% wingbox
Vwb = FusWidth*A*(chorddist( 0, s, Croot, Ctip ))^2;

% total wing fuel volume
% both wings, 85% usable volume between spars, 4% lost to tank structure,
% 5% lost to expansion space

Veff = (2*0.85*0.96*0.95*Vtot) + Vwb;
Fmasswing = 775*Veff;

msec = 0.85*0.96*0.95*775*Vsec;
x = linspace(FusWidth,s,sections);
figure
plot(x,msec)
