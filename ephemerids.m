function ef=ephemerids(n,t)
switch n
    case 1,
        theta0=0;
        OMEGA=0;
    case 2,
        theta0=-45;
        OMEGA=0;
    case 3,
        theta0=-90;
        OMEGA=0;
    case 4,
        theta0=-135;
        OMEGA=0;
    case 5,
        theta0=-180;
        OMEGA=0;
    case 6,
        theta0=-225;
        OMEGA=0;
    case 7,
        theta0=-270;
        OMEGA=0;
    case 8,
        theta0=-315;
        OMEGA=0;
    case 9,
        theta0=15;
        OMEGA=120;
    case 10,
        theta0=-30;
        OMEGA=120;
    case 11,
        theta0=-75;
        OMEGA=120;
    case 12,
        theta0=-120;
        OMEGA=120;
    case 13,
        theta0=-165;
        OMEGA=120;
    case 14,
        theta0=-210;
        OMEGA=120;
    case 15,
        theta0=-255;
        OMEGA=120;
    case 16,
        theta0=-300;
        OMEGA=120;
    case 17,
        theta0=30;
        OMEGA=240;
    case 18,
        theta0=-15;
        OMEGA=240;
    case 19,
        theta0=-60;
        OMEGA=240;
    case 20,
        theta0=-105;
        OMEGA=240;
    case 21,
        theta0=-150;
        OMEGA=240;
    case 22,
        theta0=-195;
        OMEGA=240;
    case 23,
        theta0=-240;
        OMEGA=240;
    case 24,
        theta0=-285;
        OMEGA=240;
    otherwise, disp('Incorrect # of Sat.')
end
theta0=theta0*pi/180;
OMEGA=OMEGA*pi/180;
% ci=cos(1.131);
% si=sin(1.131);
ci=cos(deg2rad(55));  % for GPS
si=sin(deg2rad(55)); 
Tsat=t*1.552448385e-4;
Tear=-t*7292115e-11;
r=25478136;
V1=-r*7292115e-11; V2=r*1.552448385e-4;

X=r*(cos(theta0+Tsat).*cos(OMEGA+Tear)-sin(theta0+Tsat).*sin(OMEGA+Tear)*ci);
Y=r*(cos(theta0+Tsat).*sin(OMEGA+Tear)+sin(theta0+Tsat).*cos(OMEGA+Tear)*ci);
Z=r*sin(theta0+Tsat)*si;

VX=-V2*(sin(theta0+Tsat).*cos(OMEGA+Tear)+cos(theta0+Tsat).*sin(OMEGA+Tear)*ci)...
    -V1*(cos(theta0+Tsat).*sin(OMEGA+Tear)+sin(theta0+Tsat).*cos(OMEGA+Tear)*ci);
VY=-V2*(sin(theta0+Tsat).*sin(OMEGA+Tear)-cos(theta0+Tsat).*cos(OMEGA+Tear)*ci)...
    +V1*(cos(theta0+Tsat).*cos(OMEGA+Tear)-sin(theta0+Tsat).*sin(OMEGA+Tear)*ci);
VZ=V2*cos(theta0+Tsat)*si;
ef=[X; Y; Z; VX; VY; VZ];