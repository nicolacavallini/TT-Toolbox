function res = evaluate_residual(matrix,rhs,sol)
    res = mtimes(matrix,sol);
    res = minus(rhs,res);
return
end