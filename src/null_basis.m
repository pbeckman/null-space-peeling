function Q = null_basis(M, k)
    [Qt, ~, ~] = svd(M'); % using svd instead of QR orders the basis by significance
    Q = Qt(:,end-k+1:end);
end