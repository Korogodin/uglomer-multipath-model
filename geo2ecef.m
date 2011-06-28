function res=geo2ecef(sv)
R=6367443.5;
x=(R+sv(3))*cos(sv(1))*cos(sv(2));
y=(R+sv(3))*cos(sv(2))*sin(sv(1));
z=(R+sv(3))*sin(sv(2));
Cg2e=[-sin(sv(1))  -sin(sv(2))*cos(sv(1))  cos(sv(1))*cos(sv(2));
       cos(sv(1))  -sin(sv(2))*sin(sv(1))  cos(sv(2))*sin(sv(1));
       0            cos(sv(2))             sin(sv(2))];
V=Cg2e*sv(4:6);
A=Cg2e*sv(7:9);
res=[x;
     y;
     z;
     V;
     A];
