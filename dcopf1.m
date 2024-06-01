%function dcopf
mpc=loadcase('case30');
Pl=mpc.bus(:,3);
Pgmin=mpc.gen(:,10);
Pgmax=mpc.gen(:,9);
a=mpc.gencost(:,end-2);
b=mpc.gencost(:,end-1);
c=mpc.gencost(:,end-1);
Ng=length(Pgmax);
mpc.branch(:,6)=25;
Fmax=mpc.branch(:,6);


H = makePTDF(mpc,1);


Pg=sdpvar(Ng,1);


obj=a'*Pg+ones(Ng,1)'*b;
obj=0;
for i=1:Ng
   obj=obj+a(i)*Pg(i)^2+b(i)*Pg(i)+c(i);
end
con=[sum(Pg)==sum(Pl)];
con=con+[Pgmin<=Pg<=Pgmax];
con=con+[-Fmax<=H(:,mpc.gen(:,1))*Pg-H*Pl<=Fmax];
ops = sdpsettings('solver','gurobi');
solvesdp(con,obj,ops)
double(obj)
% double(Pg)

PP=double(Pg)

FF=double(H(:,mpc.gen(:,1))*Pg-H*Pl);
FFbind=find(abs(FF+Fmax)<=0.0000001|abs(FF-Fmax)<=0.0000001)
PPmar=find((PP~=Pgmin)&(PP~=Pgmax))






