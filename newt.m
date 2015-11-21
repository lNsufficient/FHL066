nodes %denna kör automatiskt data.m

TOL=1e-4;

plot_a = zeros(nbr_steps, 1);
plot_f = zeros(nbr_steps, 1);

for n = 1:nbr_steps
    f = f + df;
    %u = u
    K=0*K;
   
    %Skaffa hela K, för alla element.
    for j = 1:nelm
        index=Edof(j,2:end); %de noder stången gränsar till
        ec=coord(index,:)';
        % ATT GÖRA
        % plocka ut dof med node_dof
        % 6 st dof, inte bara 2 noder
% därför blir det helt in i helvete fel
%
%
%
        ed=a(node_dof(index)); 
        [es,~]=bar3gs(ec,ep,ed);
        Ke = bar3ge(ec,ep,ed,es);
        K(index,index)=K(index,index)+Ke;
    end
    
    
    res=TOL+1;
    while res > TOL
        r = f;
        for j = 1:nelm
            index=Edof(j,2:end); %de noder stången gränsar till
            ec=coord(index,:)';
            [es, ~] = bar3gs(ec,ep,ed);
            ef=bar3gf(ec,ed,es);
            r(index) = r(index) - ef;
        end
        da = solveq(K, r, bc);
        a=a+da;
        r(bc(:,1))=0;
        res = norm(r);
        for k=1:3
        coord(:,k) = coord0(:,k)+a(k:3:(end+k-3));
        end
    end
    plot_f(n) = f(top_dof);
    plot_a(n) = a(top_dof);
end
plot(f, a)