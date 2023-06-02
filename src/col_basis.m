function Q = col_basis(M, k)
    [Qf, ~] = qr(M);
    Q = Qf(:,1:k);
end