function [x, somedata]=dmrg_solve2(A, y, tolerance)

% Inner parameters
max_full_size=20000;
prec_compr=1e-3;
prec_tol=1e-1;
prec_iters=10;

dropsweeps=1;
ddpow = 0.1; % stepsize for d-power in truncations
min_dpow = 1; % Minimal d-power for truncation
ddrank = 1; % stepsize for additional rank
min_drank = 1; % Minimal additional rank
d_pow_check = 0; % d-power for checking the convergence
bot_conv = 0.1; % bottom convergence factor - if better, we can decrease dpow and drank
top_conv = 0.99; % top convergence factor - if worse, we have to increase dpow and drank

bs_treshold = 0.0001*0; % Treshold from the previous residual to consider a local system as "bad"
trunc_to_true = 2; % Truncation error to true residual treshold

use_self_prec=false;
nswp=10;
nrestart=40;
gmres_iters=2;
% local_prec = 'als';
local_prec = 'selfprec';
local_format = 'full';
% local_format = 'tt';
rmax=1000;
tol=tolerance;
verb=1;
kickrank = 2;
x0=y;
x=x0;
P = core(tt_eye(tt_size(y), d));;






d=size(A,1);



x{1}=reshape(x{1}, size(x{1},1), 1, size(x{1},2));
y{1}=reshape(y{1}, size(y{1},1), 1, size(y{1},2));
A{1}=reshape(A{1}, size(A{1},1),size(A{1},2), 1, size(A{1},3)); %Bydlocode (@)
P{1}=reshape(P{1}, size(P{1},1),size(P{1},2), 1, size(P{1},3));

phA = cell(d,1);
phy = cell(d,1);
dx_old = ones(d,1);
dx = zeros(d,1);
% artificial rank additions
drank = ones(d,1)*min_drank;
% d-power for stronger compression tolerance./(d.^dpows)
dpows = ones(d,1)*min_dpow;

%  chkvec = tt_random(tt_size(y), max(size(y)), kickrank);
%  chkvec{1}=reshape(chkvec{1}, size(chkvec{1},1), 1, size(chkvec{1},2));
%  phAchk = cell(d,1);
%  phychk = cell(d,1);


somedata = cell(4,1); % swp, conds
somedata{2} = zeros(d, nswp);

sol_hist = cell(3,1);
sol_hist{1}=x;
max_res_old = 0;
last_sweep = false;
for swp=1:nswp
%     z = x;
    % 1-to-d orthogonalization
    rvx = 1; rnewx=1; phAold=1; phyold=1;

    for i=1:d-1
        cre = x{i};
        n1 = size(cre,1); rx1 = size(cre,2); rx2 = size(cre,3);
        cre = reshape(permute(cre, [2 1 3]), rx1, n1*rx2);
        cre = rvx*cre; % size rnew,n1,rx2
        rx1=rnewx;
        cre = reshape(cre, rx1*n1, rx2);
        [q,rvx]=qr(cre,0); % size rx1*n1,r2new - r2new,rx2
        rnewx = min(rx1*n1, rx2);
        x{i}=permute(reshape(q, rx1, n1, rnewx), [2 1 3]);

        % Now, update phi. phA=X' PA X, phY = X' PY
        a1 = A{i};
        n1=size(a1,1); m1=size(a1,2); ra1=size(a1,3); ra2=size(a1,4);
        p1 = P{i};
        k1=size(p1,1); rp1=size(p1,3); rp2=size(p1,4);
        y1 = y{i};
        ry1=size(y1,2); ry2=size(y1,3);
        x1 = x{i};

        rxm1=size(x1,2); rxm2=size(x1,3); rxn1=rxm1; rxn2=rxm2;
        phAold = reshape(phAold, rxn1*rp1*ra1, rxm1);
        x1 = reshape(permute(x1, [2 1 3]), rxm1, m1*rxm2);
        phAold=phAold*x1; % size rxn1*rp1*ra1*m1*rxm2
        phAold=reshape(phAold, rxn1, rp1, ra1, m1, rxm2);
        phAold=reshape(permute(phAold, [1 2 5 4 3]), rxn1*rp1*rxm2, m1*ra1);
        a1 = reshape(permute(a1, [2 3 1 4]), m1*ra1, n1*ra2);
        phAold=phAold*a1; % size rxn1*rp1*rxm2*n1*ra2
        phAold=reshape(phAold, rxn1, rp1, rxm2, n1, ra2);
        phAold=reshape(permute(phAold, [1 3 5 4 2]), rxn1*rxm2*ra2, n1*rp1);
        p1 = reshape(permute(p1, [2 3 1 4]), n1*rp1, k1*rp2);
        phAold=phAold*p1; % size rxn1*rxm2*ra2*k1*rp2
        phAold=reshape(phAold, rxn1, rxm2, ra2, k1, rp2);
        phAold=reshape(permute(phAold, [2 5 3 4 1]), rxm2*rp2*ra2, k1*rxn1);
        x1 = reshape(x{i}, k1*rxn1, rxn2);
        phAold=phAold*conj(x1); % size rxm2*rp2*ra2*rxn2 <--- complex conjugate!
        phAold = reshape(phAold, rxm2, rp2, ra2, rxn2);
        phAold = permute(phAold, [4 2 3 1]); % we need rxn2,rp2,ra2,rxm2
        phA{i}=phAold;

        phyold = reshape(phyold, rxn1*rp1, ry1);
        y1 = reshape(permute(y1, [2 1 3]), ry1, n1*ry2);
        phyold=phyold*y1; % size rxn1*rp1*n1*ry2
        phyold = reshape(phyold, rxn1, rp1, n1, ry2);
        phyold=reshape(permute(phyold, [1 4 3 2]), rxn1*ry2, n1*rp1);
        p1=reshape(permute(P{i}, [2 3 1 4]), n1*rp1, k1*rp2);
        phyold=phyold*p1; % size rxn1*ry2*k1*rp2
        phyold=reshape(phyold, rxn1, ry2, k1, rp2);
        phyold=reshape(permute(phyold, [4 2 3 1]), rp2*ry2, k1*rxn1);
        x1=reshape(x{i}, k1*rxn1, rxn2);
        phyold=phyold*conj(x1); % size rp2*ry2*rxn2 <--- complex conjugate!
        phyold=permute(reshape(phyold, rp2, ry2, rxn2), [3 1 2]);
        phy{i}=phyold;
    end;
    % convolve rv with the last cre
    cre = x{d};
    n1 = size(cre,1); rx1 = size(cre,2); rx2 = size(cre,3);
    cre = reshape(permute(cre, [2 1 3]), rx1, n1*rx2);
    cre = rvx*cre; % size rnew,n1,rx2
    x{d}=permute(reshape(cre, rnewx, n1, rx2), [2 1 3]);

    % Now, start the d-to-1 DMRG iteration
    max_res = 0;
    phAold=1; phyold=1;

    for i=d:-1:2
        a2=A{i}; a1=A{i-1}; ra1=size(a1,3); ra2=size(a1,4); ra3=size(a2,4);
        n1 = size(a1,1); m1=size(a1,2); n2=size(a2,1); m2=size(a2,2);
        p2=P{i}; p1=P{i-1}; rp1=size(p1,3); rp2=size(p1,4); rp3=size(p2,4);
        k1 = size(p1,1); k2=size(p2,1);

        y1=y{i-1}; y2=y{i}; ry1=size(y1,2); ry2=size(y1,3); ry3=size(y2,3);
        x1=x{i-1}; x2=x{i}; rx1=size(x1,2); rx2=size(x1,3); rx3=size(x2,3);

        % Compute RHS: phy{i-2}*P1*y1*y2*P2*phyold
        if (i>2)
            rhs1 = phy{i-2};
        else
            rhs1=1;
        end;
        rhs1 = reshape(rhs1, rx1*rp1, ry1);
        y1 = reshape(permute(y1, [2 1 3]), ry1, n1*ry2);
        rhs1 = rhs1*y1; % size rx1*rp1*n1*ry2
        rhs1 = reshape(rhs1, rx1, rp1, n1, ry2);
        rhs1=reshape(permute(rhs1, [1 4 3 2]), rx1*ry2, n1*rp1);
        p1 = reshape(permute(p1, [2 3 1 4]), n1*rp1, k1*rp2);
        rhs1=rhs1*p1; % size rx1*ry2*k1*rp2
        rhs1=reshape(rhs1, rx1, ry2, k1, rp2);
        rhs1=reshape(permute(rhs1, [1 3 4 2]), rx1*k1, rp2*ry2);

        y2=reshape(permute(y2, [2 1 3]), ry2*n2, ry3);
        phyold2 = reshape(phyold, rx3*rp3, ry3);
        rhs2 = y2*(phyold2.'); % size ry2*n2, rx3*rp3
        rhs2 = reshape(rhs2, ry2, n2, rx3, rp3);
        rhs2 = permute(rhs2, [1 3 2 4]);
        rhs2 = reshape(rhs2, ry2*rx3, n2*rp3);
        p2 = reshape(permute(p2, [2 4 1 3]), n2*rp3, k2*rp2);
        rhs2 = rhs2*p2; % size ry2*rx3, k2*rp2
        rhs2 = reshape(rhs2, ry2, rx3, k2, rp2);
        rhs2 = permute(rhs2, [4 1 3 2]);
        rhs2 = reshape(rhs2, rp2*ry2, k2*rx3);

        if (strcmp(local_format, 'full'))
            rhs = rhs1*rhs2;
            rhs = reshape(rhs, rx1*k1*k2*rx3, 1);
        else
            rhs = cell(2,1);
            rhs{1} = rhs1;
            rhs{2} = rhs2.';
        end;

        rxn1=rx1; rxn3=rx3;
        rxm1=rx1; rxm2=rx2; rxm3=rx3;
        if (i>2)
            B = phA{i-2};
        else
            B=1;
        end;
              
        if (i>2)
            B = phA{i-2};
        else
            B=1;
        end;        
        
        B = reshape(B, rxn1, rp1, ra1, rxm1);
        B = reshape(permute(B, [1 4 2 3]), rxn1*rxm1*rp1, ra1);
        a1 = reshape(permute(A{i-1}, [3 2 1 4]), ra1, m1*n1*ra2);
        B = B*a1; % size rxn1*rxm1*rp1*m1*n1*ra2
        B = reshape(B, rxn1,rxm1,rp1,m1,n1,ra2);
        B = permute(B, [1 2 4 6 3 5]);
        B = reshape(B, rxn1*rxm1*m1*ra2, rp1*n1);
        p1 = reshape(permute(P{i-1}, [3 2 1 4]), rp1*n1, k1*rp2);
        B = B*p1; % size rxn1*rxm1*m1*ra2*k1*rp2
        B = reshape(B, rxn1,rxm1,m1,ra2,k1,rp2);
        B = permute(B, [1 5 2 3 6 4]);
        B = reshape(B, rxn1*k1*rxm1*m1, rp2*ra2);
        % This is the first term of tensor-structured matrix B \otimes B2
        %  Now, the second
        B2 = permute(phAold, [1 2 4 3]);
        B2 = reshape(B2, rxn3*rp3*rxm3, ra3);
        a2 = reshape(A{i}, n2*m2*ra2, ra3);
        B2 = B2*(a2.'); % size rxn3*rp3*rxm3*n2*m2*ra2
        B2 = reshape(B2, rxn3, rp3, rxm3, n2, m2, ra2);
        B2 = permute(B2, [1 3 5 6 4 2]);
        B2 = reshape(B2, rxn3*rxm3*m2*ra2, n2*rp3);
        p2 = reshape(permute(P{i}, [2 4 1 3]), n2*rp3, k2*rp2);
        B2 = B2*p2; % size rxn3*rxm3*m2*ra2*k2*rp2
        B2 = reshape(B2, rxn3, rxm3, m2, ra2, k2, rp2);
        B2 = permute(B2, [5 1 3 2 6 4]);
        B2 = reshape(B2, k2*rxn3*m2*rxm3, rp2*ra2);
        
        rB=rp2*ra2;
        MatVec='bfun2';
        if (((rxn1*k1*k2*rxn3<max_full_size))||(rB>max(rxn1*k1, rxn3*k2)))&&(strcmp(local_format, 'full'))
            MatVec='full';
            if (rxn1*k1*k2*rxn3>max_full_size)
                MatVec='half-full';
            end;
            B = B*(B2.'); % size rxn1*k1*rxm1*m1*k2*rxn3*m2*rxm3
            B = reshape(B, rxn1,k1,rxm1,m1,k2,rxn3,m2,rxm3);
            B = permute(B, [1 2 5 6 3 4 7 8]);
            B = reshape(B, rxn1*k1*k2*rxn3, rxm1*m1*m2*rxm3);
        else
            B1 = reshape(B, rxn1*k1, rxm1*m1, rB);
            B2 = reshape(B2, k2*rxn3, m2*rxm3, rB);
            B=cell(2,1);
            B{1}=B1;
            B{2}=B2;
        end;

        % Form previous solution
        x1 = reshape(permute(x{i-1}, [2 1 3]), rxm1*m1, rxm2);
        x2 = reshape(permute(x{i}, [2 1 3]), rxm2, m2*rxm3);
        if (strcmp(local_format, 'full'))
            sol_prev = x1*x2;
            sol_prev = reshape(sol_prev, rxm1*m1*m2*rxm3, 1);
        else
            sol_prev = cell(2,1);
            sol_prev{1} = x1;
            sol_prev{2} = x2.';
        end;

        real_tol = (tol/(d^dpows(i)))/trunc_to_true;

        if (strcmp(local_format, 'tt'))
            mv = @(vec,tolerance,mr)bfun3(B,vec,tolerance,mr);
        else
            if (strcmp(MatVec, 'bfun2'))
                mv=@(vec)bfun2(B, vec, rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
                mv_t=@(vec)bfun2_t(B, vec, rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
                
            else
                mv = @(vec)(B*vec);
                mv_t = @(vec)(B'*vec);
                
            end;
        end;

        % Check the previous residual
        if (strcmp(local_format, 'tt'))
            res_prev = mv(sol_prev, [], []);
            normf = exp(0.5*tt_dot2(rhs, rhs));
            res_prev = tt_dist3(res_prev, rhs)/normf;
        else
            res_prev = mv(sol_prev);
            normf = norm(rhs);
            res_prev = norm(res_prev - rhs)/normf;
        end;

        % We will solve the system only if res_prev>0.1*max_res_prev
        if (~last_sweep)&&(res_prev>bs_treshold*max_res_old)
            if (strcmp(MatVec,'full'))
                disp("FULL LINEAR SYSTEM")
                somedata{3}(i, swp) = res_prev;
                sol = B \ rhs;
                
                res=B*sol;
                res_true = norm(res-rhs)/norm(rhs);
            else
                if (strcmp(local_format, 'full'))
                    [sol_new,flg] = gmres(mv, rhs, nrestart, real_tol, 2, [], [], sol_prev);
                    res_new=norm(mv(sol_new)-rhs)/normf;
                    conv_factor=(res_new/res_prev);
                    if (res_new*(conv_factor)>real_tol && use_self_prec && strcmp(MatVec, 'bfun2')) % we need a prec.
                        if (strcmp(local_prec, 'selfprec'))
                            iB=tt_minres_selfprec(B, prec_tol, prec_compr, prec_iters, 'right');

                            resid = rhs-mv(sol_new);
                            [dsol,flg] = gmres(@(vec)bfun2(B, bfun2(iB,vec,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3),...
                                rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3), resid, nrestart, real_tol/res_new, gmres_iters);
                            dsol = bfun2(iB,dsol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3);
                            sol = sol_new+dsol;

                        end;
                        if (strcmp(local_prec, 'als'))
                            sol = als_solve_rx_2(B, rhs, real_tol, [], sol_new);
                        end;
                    else
                        [sol,flg] = gmres(mv, rhs, nrestart, real_tol, gmres_iters, [], [], sol_new);
                    end;
                    if (flg>0)
                        fprintf('-warn- gmres did not converge at block %d\n', i);
                    end;                    

                    res=mv(sol);
                    res_true = norm(res-rhs)/normf;
                else
                    sol_new = tt_gmres(mv, rhs, real_tol, 2, nrestart, real_tol, real_tol, [], [], [], sol_prev);
                    res_new=tt_dist3(mv(sol_new,[],[]),rhs)/normf;
                    conv_factor=(res_new/res_prev);
                    if (res_new*conv_factor>real_tol && use_self_prec && strcmp(MatVec, 'bfun2')) % we need a prec.
                            iB=tt_minres_selfprec(B, prec_tol, prec_compr, prec_iters, 'right');
                            sol = tt_gmres(@(vec,tolerance,mr)bfun3(B, vec, tolerance, mr), rhs, real_tol, gmres_iters, nrestart, real_tol, real_tol, @(vec,tolerance,mr)bfun3(iB, vec, tolerance, mr), [], [], sol_new);
                    else
                        sol = tt_gmres(mv, rhs, real_tol, gmres_iters, nrestart, real_tol, real_tol, [], [], [], sol_new);
                    end;
                    res=mv(sol,[],[]);
                    res_true = tt_dist3(res,rhs)/normf;
                end;
            end;

            if (strcmp(local_format, 'full'))
                dx(i) = norm(sol-sol_prev,'fro')/norm(sol_prev,'fro');
            else
                dx(i) = tt_dist3(sol, sol_prev)/exp(0.5*tt_dot2(sol,sol));
            end;
        else
            res_true = res_prev;
            dx(i)=0;
	    sol = sol_prev;
        end;

        if ((res_true>res_prev/trunc_to_true))&&(res_true>real_tol)&&(~last_sweep)
	    fprintf('--warn-- the residual damp by gmres was smaller than in the truncation\n');
%             keyboard;
            sol = sol_prev;
            res_true = res_prev;
        end;

        if (res_prev>max_res)
            max_res = res_prev;
        end;

        
        if (strcmp(local_format, 'full'))
            nrmsol = norm(sol, 'fro');
        else
            nrmsol = exp(0.5*tt_dot2(sol,sol));
        end
        if (nrmsol==0)
            dx(i)=0;
        end;

        if (swp==1)
            dx_old(i)=dx(i);
        end;

        % The new core does not converge - increase rank
        if (dx(i)/dx_old(i)>top_conv)&&(dx(i)>tolerance/(d^d_pow_check))
            drank(i)=drank(i)+ddrank;
            dpows(i)=dpows(i)+ddpow;
        end;
        % The new core converges well - try to decrease rank
        if (dx(i)/dx_old(i)<bot_conv)||(dx(i)<tolerance/(d^d_pow_check))
            drank(i)=max(drank(i)-ddrank, min_drank);
            dpows(i)=max(dpows(i)-ddpow, min_dpow);
        end;
        % perform simple compression for the last sweep
        if (last_sweep)
            dpows(i)=min(0.5, min_dpow);
        end;

        if (res_prev>bs_treshold*max_res_old)&&(strcmp(local_format, 'full'))
            if (mod(swp,dropsweeps)~=0)&&(swp>1)&&(~last_sweep)
                [u,s,v]=svd(sol-reshape(sol_prev,[rxm1*m1,m2*rxm3]),'econ');
            else
                if (~last_sweep)
                    sol=reshape(sol,[rxm1*m1,m2*rxm3]);
                    [u,s,v]=svd(sol,'econ');
                else
                    [x2,rv]=qr(x2.', 0);
                    x1 = x1*(rv.');
                    [u,s,v]=svd(x1, 'econ');
                    v = x2*v;
                end;
            end;
            s = diag(s);
            flm=norm(s);
            %Truncation block. We have to make it smarter by binary search
            r0 = 1; rM = min(size(s,1),rmax); r = round((r0+rM)/2);
            while (rM-r0>2)
                er0=norm(s(r+1:numel(s)));
                if (mod(swp,dropsweeps)~=0)&&(swp>1)&&(~last_sweep)
                    sol = sol_prev+reshape(u(:,1:r)*diag(s(1:r))*(v(:,1:r))',rxm1*m1*m2*rxm3, 1);
                else
                    sol = reshape(u(:,1:r)*diag(s(1:r))*(v(:,1:r))',rxm1*m1*m2*rxm3, 1);
                end;
                if (strcmp(MatVec,'full')||strcmp(MatVec,'half-full'))
                    resid = norm(B*sol-rhs)/norm(rhs);
                else
                    resid = norm(bfun2(B,sol,rxm1,m1,m2,rxm3,rxn1,k1,k2,rxn3)-rhs)/norm(rhs);
                end;
                
                if ((resid<max(res_true*trunc_to_true, tolerance/(d^dpows(i)))) ) %Value of the rank is OK
                    rM = r-1;
                    r = round((r0+rM)/2);
                else %Is not OK.
                    r0 = r;
                    r = round((r0+rM)/2);
                end;
            end
            r = r0;
            % Line search - if the rank is underestimated
            cursol = cell(2,1);
            cursol{1}=u(:,1:r);
            cursol{2}=conj(v(:,1:r))*diag(s(1:r));
            if (strcmp(MatVec,'full')||strcmp(MatVec,'half-full'))
                resid = B*full(tt_tensor(cursol),rxm1*m1*m2*rxm3)-rhs;
            else
                resid = full(tt_tensor(tt_mv(B,cursol)),rxm1*m1*m2*rxm3)-rhs;
            end;
            while (r<min(size(s,1), rmax))
                r=r+1;
                er0=norm(s(r+1:numel(s)));
                cursol{1}=u(:,r);
                cursol{2}=conj(v(:,r))*s(r);
                if (strcmp(MatVec,'full')||strcmp(MatVec,'half-full'))
                    resid = B*full(tt_tensor(cursol),rxm1*m1*m2*rxm3)+resid;
                    
                else
                    resid = full(tt_tensor(tt_mv(B,cursol)),rxm1*m1*m2*rxm3)+resid;

                end;
                normres = norm(resid)/norm(rhs);
                
                if ((normres<max(res_true*trunc_to_true, tolerance/(d^dpows(i)))) ) %Value of the rank is OK
                    break;
                end;
            end;

            if (~last_sweep)
                r = r+drank(i); % we want even larger ranks
            end;

            v = conj(v);
        else
            if (strcmp(local_format, 'tt'))
                x1 = sol{1};
                x2 = sol{2}.';
            end;
            % We do not have to decimate the whole supercore,
            % only one of factors, as we have the previous solution
            [v,rv]=qr(x2.',0); % size m2*rxm3, rxm2' - rxm2',rxm2
            r = size(v,2);
            u = x1*rv.';
            s = ones(r,1);
        end;
        r = min(r, max(size(s))); % but not too large
        r = min(r,rmax);

        v = v(:,1:r);
        u = u(:,1:r)*diag(s(1:r));


        if (mod(swp,dropsweeps)~=0)&&(swp>1)&&(~last_sweep)
            u = [x1, u];
            v = [x2.', v];
            [v,rv]=qr(v,0);
            u = u*(rv.');
            r = size(v,2);
        else
            if (~last_sweep)
                vr=randn(size(v,1),kickrank);
%                 v=reort(v,vr);
                [v,rr]=qr([v,vr], 0);
%                 radd=size(v,2)-r;
                radd=kickrank;
                if ( radd > 0 )
                    ur=zeros(size(u,1),radd);
                    u=[u,ur]*(rr.');
                end
                r = size(u,2);
%                 r=r+radd;
            end;
        end;

        x{i}=permute(reshape(v, m2, rxm3, r), [1 3 2]);
        x{i-1}=permute(reshape(u, rxm1, m1, r), [2 1 3]);
        rxm2=r; rxn2=r;



        phAold = reshape(phAold, rxn3*rp3*ra3, rxm3);
        x2 = reshape(x{i}, m2*rxm2, rxm3);
        phAold = phAold*(x2.'); % size rxn3*rp3*ra3*m2*rxm2
        phAold = reshape(phAold, rxn3, rp3, ra3, m2, rxm2);
        phAold = permute(phAold, [1 2 5 4 3]);
        phAold = reshape(phAold, rxn3*rp3*rxm2, m2*ra3);
        a2 = reshape(permute(A{i}, [2 4 1 3]), m2*ra3, n2*ra2);
        phAold=phAold*a2; % size rxn3*rp3*rxm2*n2*ra2
        phAold=reshape(phAold, rxn3, rp3, rxm2, n2, ra2);
        phAold = permute(phAold, [1 3 5 4 2]);
        phAold = reshape(phAold, rxn3*rxm2*ra2, n2*rp3);
        p2 = reshape(permute(P{i}, [2 4 1 3]), n2*rp3, k2*rp2);
        phAold=phAold*p2; % size rxn3*rxm2*ra2*k2*rp2
        phAold=reshape(phAold, rxn3,rxm2,ra2,k2,rp2);
        phAold = permute(phAold, [2 3 5 1 4]);
        phAold = reshape(phAold, rxm2*ra2*rp2, rxn3*k2);
        x2 = reshape(permute(x{i}, [3 1 2]), rxn3*k2, rxn2);
        phAold = phAold*conj(x2); % size rxm2*ra2*rp2*rxn2 <-- cplx conjugate!
        phAold = permute(reshape(phAold, rxm2, ra2, rp2, rxn2), [4 3 2 1]);

        phyold = reshape(phyold, rxn3*rp3, ry3);
        y2 = reshape(y{i}, n2*ry2, ry3);
        phyold = phyold*(y2.'); % size rxn3*rp3*n2*ry2
        phyold = reshape(phyold, rxn3, rp3, n2, ry2);
        phyold = permute(phyold, [1 4 3 2]);
        phyold = reshape(phyold, rxn3*ry2, n2*rp3);
        p2 = reshape(permute(P{i}, [2 4 1 3]), n2*rp3, k2*rp2);
        phyold = phyold*p2; % size rxn3*ry2*k2*rp2
        phyold = reshape(phyold, rxn3, ry2, k2, rp2);
        phyold = permute(phyold, [4 2 1 3]);
        phyold = reshape(phyold, rp2*ry2, rxn3*k2);
        x2 = reshape(permute(x{i}, [3 1 2]), rxn3*k2, rxn2);
        phyold = phyold*conj(x2); % size rp2*ry2*rxn2 <-- cplx conjugate!
        phyold = permute(reshape(phyold, rp2, ry2, rxn2), [3 1 2]);
    end;

    max_res_old = max_res;
    if (mod(swp,100)==0)
        max_res_old = 0;
    end;

    if (last_sweep)
        break;
    end;
    if (max_res<tol/(d^d_pow_check))||(swp==nswp-1)
        last_sweep=true;
%         break;
    end;

    dx_old = dx;
end;

x{1}=reshape(x{1}, size(x{1},1), size(x{1},3));

% x = tt_compr2(x, tolerance, rmax);

if (nargout>1)
    somedata{1} = swp;
end;

end

function [y]=bfun2(B, x, rxm1, m1, m2, rxm3, rxn1, k1, k2, rxn3)
% Computes (B{1} \otimes B{2})x
% B{1} is of sizes rxn1*k1, rxm1*m1, rB
% B{2} is of sizes k2*rxn3, m2*rxm3, rB
rB=size(B{1},3);
x = reshape(x, rxm1*m1, m2*rxm3);
B1 = permute(B{1}, [3 1 2]);
B1 = reshape(B1, rB*rxn1*k1, rxm1*m1);
y = B1*x; % size rB*rxn1*k1,m2*rxm3
y = reshape(y, rB, rxn1*k1, m2*rxm3);
y = permute(y, [3 1 2]);
y = reshape(y, m2*rxm3*rB, rxn1*k1);
B2 = reshape(B{2}, k2*rxn3, m2*rxm3*rB);
y = B2*y; % size k2*rxn3,rxn1*k1
y = reshape(y.', rxn1*k1*k2*rxn3, 1);
end


function [y]=bfun2_t(B, x, rxm1, m1, m2, rxm3, rxn1, k1, k2, rxn3)
% Computes (B{1}' \otimes B{2}')x
% B{1} is of sizes rxn1*k1, rxm1*m1, rB
% B{2} is of sizes k2*rxn3, m2*rxm3, rB
rB=size(B{1},3);
x = reshape(x, rxm1*m1, m2*rxm3);
B1 = permute(B{1}, [3 2 1]);  % HERE
B1 = reshape(B1, rB*rxn1*k1, rxm1*m1);
y = B1*x; % size rB*rxn1*k1,m2*rxm3
y = reshape(y, rB, rxn1*k1, m2*rxm3);
y = permute(y, [3 1 2]);
y = reshape(y, m2*rxm3*rB, rxn1*k1);
B2 = reshape(permute(B{2}, [2,1,3]), k2*rxn3, m2*rxm3*rB);
y = B2*y; % size k2*rxn3,rxn1*k1
y = reshape(y.', rxn1*k1*k2*rxn3, 1);
end


function [y] = bfun3(A, x, tolerance, mr)
% For the 2d--TT MatVec
y = tt_mv(A, x);
if (nargin<4)
    mr = [];
end;
if (nargin>2)&&(~isempty(tolerance))
    y = tt_compr2(y, tolerance, mr);
end;

end

