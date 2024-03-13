function sol = evaluate_sol(some_norm,dim,some_tmp_sol,previus_sol)
    % Recover the scales
    % Distribute norms equally...
    some_other_norm = exp(sum(log(some_norm))/dim);
    % ... and plug them into x
    for i=1:dim
        some_tmp_sol{i} = some_tmp_sol{i}*some_other_norm;
    end;

    sol = cell2core(previus_sol, some_tmp_sol);
return
end