function [C, r] = sqrt_truncate_r(S, tau)

%@S: n x n Similarity matrix
%@r: rank tuning parameter
%info: finds square root of S using eigen decompostion and truncates to
%rank r

    %[U,Z] = eig(S);
    [U, Z, V] = svd(S);

    z = diag(Z);
    n = size(U,2);

    for r = 1:n
        C = U(:,1:r)*diag(sqrt(z(1:r)));
        objfunc = (norm(S-C*C','fro')/norm(S,'fro'));
        if objfunc <= tau
            break
        end
    end
end
