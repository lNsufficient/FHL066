data %innehåller värden för att hitta noderna vettigt
n1 = [L;L*tand(30);h_2]; %y-koordinaten fås genom en tänkt bisektris från 
%nod 7 till nod 1
n2 = [L*3/2; L*sind(60);h_1];
n3 = [L*cosd(60); L*sind(60);h_1];
n4 = [L;0;h_1];
n5 = [2*L;0;0];
n6 = [L; 2*L*sind(60);0];
n7 = [0;0;0];
coord = [n1'; n2'; n3'; n4'; n5'; n6'; n7'];
coord0=coord;
Enod = [1 1 2;
        2 1 3;
        3 1 4;
        4 2 3;
        5 3 4;
        6 2 4;
        7 2 6;
        8 3 6;
        9 3 7;
        10 4 7;
        11 4 5;
        12 2 5];
    
    Edof=[Enod(:,1), node_dof(Enod(:,2)),node_dof(Enod(:,3))];
nelm = length(Enod);
nnod = length(coord);
ndof = nnod*3;
a=zeros(ndof,1);

top_dof=3; %z-förskj. i nod 1

P_end=-5e-2; %slutgiltig
nbr_steps=1000;

f=zeros(ndof,1);
df=f;
df(top_dof)=P_end/nbr_steps;

K=zeros(ndof);

bc=[];
for i=5:7
    bc = [bc; node_dof(i)', zeros(3,1)];
end

[Ex0,Ey0,Ez0]=coordxtr(Edof,coord0,node_dof((1:nnod)'),2);
eldraw3(Ex0,Ey0,Ez0,[1 4 1]);
