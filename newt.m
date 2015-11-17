nodes %denna kör automatiskt data.m



for n = 1:nbr_steps
    f = f + df;
    %u = u
    K=0*K;
    for j = 1:nelm
        index=Edof(j,2:); %de noder stången gränsar till
        ec=coord(index,:)';
        ed=a(node_dof(index)); 
        [es,~]=bar3gs(ec,ep,ed);
        Ke = bar3ge(ec,ep,ed,es);
        K(index,index)=K(index,index)+Ke;
    end
    while res > TOL
        
    end
        
end
