perturb_switch = 0;
%0. ingen st�rn
%1. st�rn i geometri
%2. st�rn i lasten

nodes %denna kÃ¶r automatiskt data.m


%plot_a = zeros(nbr_steps, 1);
%plot_f = zeros(nbr_steps, 1);

plot_a = [];
plot_f = [];

lambda=0;
S=0;

P = zeros(ndof,1);
P_end=-.15;
if perturb_switch == 2
    P(top_dof) = P_end*cosd(12);
    P(top_dof-1) = P_end*sind(12);
else
    P(top_dof) = P_end;
end
G = P*0;
f_int=P*0;


%TOL=norm(P_end/nbr_steps)*1e-2;
l=4e-4;
l_0=l;
TOL = l*1e-3;

TOL_0 = TOL; %what should be considered a zero for pq
%TOL - what should be considered a zero for res


psy=1/44;
Lamb = [0 0 0];
n=0;

old_a = a; 

USE_HOOKE = 0;
continuous_plot = 0;

while lambda < 1
    n=n+1   
    
    a_i = a;
    lambda_i = lambda;
    S_i = S;
    
    res = TOL+1;
    i = 1;
    
    
    G=G*0;
    
    while (res > TOL)
        
        K=K*0;
        f_int=f_int*0;
        
        
        for j = 1:nelm
            index_dof=Edof(j,2:end); %de frihg stÃ¥ngen grÃ¤nsar till
            index_nod=Enod(j,2:end); % de noder st gr t
            ec=coord0(index_nod,:)';
            
            ed=a_i(index_dof);
            if USE_HOOKE
                [es,~]=bar3gs(ec,ep,ed);
                Ke = bar3ge(ec,ep,ed,es);
            else
                [es,eps]=bar3gs_log(ec,ep,ed);
                Ke = bar3ge_log(ec,ep,ed,es,eps);
            end
            K(index_dof,index_dof)=K(index_dof,index_dof)+Ke;
            ef=bar3gf(ec,ed,es);
            f_int(index_dof) = f_int(index_dof) + ef;
        end
        
        da_G = solveq(K,-G,bc);
        da_P = solveq(K,P,bc);
        
        %[n,norm(da_G)]
        
        delta_a=a_i-a;
        delta_lambda=lambda_i-lambda; 
        [a1, a2, a3, a4, a5] = calc_a(P, delta_a, delta_lambda, da_P, da_G, l, psy);
        
        if i > 1
            %             if (abs(a2^2-a3*2*a1)<TOL_0) %Detta motsvarar att båda rötterna
            %                 %är samma.
            %                 d_lambda = -a2/(2*a1);
            %                 Lamb = Lamb + [1 0 0];
            %             elseif ((a2/(2*a1)^2-a3/a1) < 0) %Här räcker det att kolla om det är
            %                 %mindre än noll, -TOL_0 täcks av if-satsen
            %                 disp('Complex Solutions')
            %                 Lamb = Lamb + [0 1 0];
            %             else

            
            %s=sign(delta_a'*da_P);
            p=a2/a1;        q=a3/a1;
            d_lambda_test=-p/2+sqrt((p/2)^2-q)*[1; -1];
            if (a4+a5*d_lambda_test(1) > a4+a5*d_lambda_test(2))
                d_lambda_test=d_lambda_test(1);
            else
                d_lambda_test=d_lambda_test(2);
            end

            %                 Lamb = Lamb + [0 0 1];
            %             end
            
            temp=1;
        else
            temp=2;
            delta_a_n = a - old_a;
            s=sign(delta_a_n'*da_P);
            if (n == 1)
                s = 1;
            end
            d_lambda_test=s*sqrt(-a3/a1);
            
%              if (a4+a5*d_lambda_test < a4-a5*d_lambda_test)
%                  d_lambda_test=-d_lambda_test;
%              end
            
            s=1;
        end   
        if ~isreal(d_lambda_test)
            keyboard
            l = l/2;
            continue;
        else
            d_lambda = d_lambda_test;
            l = min(l*2,l_0);
        end
            
        lambda_i = lambda_i + d_lambda; 
        
        a_i = a_i + da_G + da_P*d_lambda;
        
        G=f_int-lambda_i*P;
        G(bc(:,1)) = 0;
        res=norm(G);
        i=i+1;
        
        if mod(i,100)==0
        %    keyboard;
            break;
        end
        
        
    end
    old_a = a;
    a=a_i;

    if n==30
    oldLambda = lambda-lambda_i;
    end
    
    
    if (n > 30 &&((lambda_i-lambda)/oldLambda > 2))
        keyboard
    end
    %oldLambdaChange = lambdaChange;
    lambda=lambda_i;
   
    plot_f(n) = lambda*P_end;
    plot_a(n) = a(top_dof);
    
    if mod(n,5000)==0
        keyboard;
    end
    
    if continuous_plot
        for k=1:3
                coord(:,k) = coord0(:,k)+a(k:3:(end+k-3));
        end
        [Ex,Ey,Ez]=coordxtr(Edof,coord,node_dof((1:nnod)'),2);
        clf;
        eldraw3(Ex,Ey,Ez,[1 4 1]);
        pause;
    end
end


%%
%plot(abs(plot_a), abs(plot_f))
d_lambda = (plot_f(2:end) - plot_f(1:end-1))/P_end;
plot(real(d_lambda))
plot(plot_a, plot_f)
xlabel('förskjutning / meter')
ylabel('kraft / Newton')
A = zeros(length(plot_a), 3);
F = zeros(length(plot_f), 3);
A(:,perturb_switch+1) = plot_a;
F(:,perturb_switch+1) = plot_f;

 for k=1:3
            coord(:,k) = coord0(:,k)+a(k:3:(end+k-3));
 end
[Ex,Ey,Ez]=coordxtr(Edof,coord,node_dof((1:nnod)'),2);
eldraw3(Ex,Ey,Ez,[1 4 1]);
