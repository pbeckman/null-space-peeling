function Q = col_basis(M, k)
    if nargin < 2
        k = size(M,2);
    end
    [Qf, ~, ~] = svd(M); % using svd instead of QR orders the basis by significance
    Q = Qf(:,1:k);
end