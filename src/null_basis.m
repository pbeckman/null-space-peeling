function Q = null_basis(M, k)
    [Qt, ~] = qr(M');
    Q = Qt(:,end-k+1:end);
end