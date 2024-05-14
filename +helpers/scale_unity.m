function [Cs, rhss] = scale_unity(C,rhs)
% Intuitive scaling loosely based on the 1-Frobenius norm of C
cnorm = full(sum(sum(abs(C)))/nnz(C));
if cnorm > 0
    Cs = C/cnorm;
    rhss = rhs/cnorm;
    fs = "Scaled matrix to unity, cnorm = %f \n";
    fprintf(fs, cnorm)
end
end