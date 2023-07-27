function dcopf2
mpc=loadcase('case118');
Pl=mpc.bus(:,3);

Pgmin=mpc.gen(:,10);
Pgmax=mpc.gen(:,9);
mpc.branch(:,6)=200;
Fmax=mpc.branch(:,6);

c=mpc.gencost(:,end-2);
a=mpc.gencost(:,end-1);
b=mpc.gencost(:,end);

Ng=length(Pgmax);
Nl=length(Pl);
Nline=length(Fmax);

B=makeBdc(mpc);

tap = ones(Nline, 1);                              %% default tap ratio = 1
i = find(mpc.branch(:, 9));                       %% indices of non-zero tap ratios
tap(i) = mpc.branch(i, 9);                        %% assign non-zero tap ratios

%% build connection matrix Cft = Cf - Ct for line and from - to buses
f = mpc.branch(:, 1);                           %% list of "from" buses
t = mpc.branch(:, 2);                           %% list of "to" buses
i = [(1:Nline)'; (1:Nline)'];                         %% double set of row indices
S = sparse(i, [f;t], [ones(Nline, 1); -ones(Nline, 1)], Nline, Nl);    %% connection matrix


X=sparse((1:Nline)',(1:Nline)',1./mpc.branch(:,4)./tap,Nline, Nline)
A=sparse( mpc.gen(:,1), (1:Ng)',ones(Ng, 1), Nl, Ng);

% X=diag(1./mpc.branch(:,4)./tap);



theta=sdpvar(Nl,1);
Pg=sdpvar(Ng,1);
obj=a'*Pg+ones(Ng,1)'*b;
con=[theta(1)==0]+[B*theta==A*Pg-Pl+S'*X*(mpc.branch(:,10))];
con=con+[Pgmin<=Pg<=Pgmax];
con=con+[-Fmax<=X*(S*theta-mpc.branch(:,10))<=Fmax];

ops = sdpsettings('solver','gurobi');
solvesdp(con,obj,ops)

double(obj)
% double(Pg)
