nodes %denna k√∂r automatiskt data.m

TOL=norm(df)*1e-3;

plot_a = zeros(nbr_steps, 1);
plot_f = zeros(nbr_steps, 1);

lambda=0;
S=0;

P = zeros(ndof,1);
P(top_dof) = P_end;
G = P*0;
f_int=P*0;

TOL_0 = TOL; %what should be considered a zero for pq
%TOL - what should be considered a zero for res

l=1e-5;

psy=0;

for n = 1:nbr_steps
    a_i = a;
    lambda_i = lambda;
    S_i = S;
    
    res = TOL+1;
    i = 1;
    while (res > TOL)
        
        for j = 1:nelm
            index_dof=Edof(j,2:end); %de frihg st√•ngen gr√§nsar till
            index_nod=Enod(j,2:end); % de noder st gr t
            ec=coord0(index_nod,:)';
            
            ed=a(index_dof);
            [es,~]=bar3gs(ec,ep,ed);
            Ke = bar3ge(ec,ep,ed,es);
            K(index_dof,index_dof)=K(index_dof,index_dof)+Ke;
            ef=bar3gf(ec,ed,es);
            f_int(index_dof) = f_int(index_dof) + ef;
        end
        
        da_G = solveq(K,-G,bc);
        da_P = solveq(K,P,bc);
        
        delta_a=a_i-a;
        delta_lambda=lambda_i-lambda;
        
        [a1, a2, a3, a4, a5] = calc_a(P, delta_a, delta_lambda, da_P, da_G, l, psy);
        
        if (abs(a2^2-a3*2*a1)<TOL_0) %Detta motsvarar att bÂda rˆtterna
            %‰r samma.
            d_lambda = -2/(2*a1);
        elseif (a2^2-2*a3*a1 < 0) %H‰r r‰cker det att kolla om det ‰r 
            %mindre ‰n noll, -TOL_0 t‰cks av if-satsen
            disp('Complex Solutions')
        else
            d_lambda=sign(delta_a'*da_P)*l/(sqrt(da_P'*da_P+psy*P'*P)); 
        end
        
        lambda_i = lambda_i + d_lambda;
         
        a_i = a_i + da_G + da_P*d_lambda;
        
        G=f_int-lambda_i*P;
        i=i+1;
        
        if mod(i,100)==0
            keyboard;
        end
    end
    
    plot_f(n) = f(top_dof);
    plot_a(n) = a(top_dof);
    
end
    