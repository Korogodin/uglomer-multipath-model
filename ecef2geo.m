function res=ecef2geo(sv)
R=6367443.5;
L=atan2(sv(2),sv(1));
B=atan2(sv(3), norm(sv(1:2)));
h=norm(sv(1:3))-R;
Cg2e=[-sin(L)  -sin(B)*cos(L)  cos(L)*cos(B);
       cos(L)  -sin(B)*sin(L)  cos(B)*sin(L);
       0            cos(B)           sin(B)];
V=Cg2e'*sv(4:6);
A=Cg2e'*sv(7:9);
res=[L;
     B;
     h;
     V;
     A];
