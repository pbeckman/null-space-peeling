function left = plot_factor(fig, F, marg, sc, left)
    S = sparse(F);
    S(S==0) = NaN;
    width  = sc*size(S,2);
    subplot('Position', [left, marg, width, 1-marg], 'Parent', fig);
    imagesc(S)
    colormap(flipud(gray))
    left = left + width + marg;
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    pbaspect([size(S,2) size(S,1) 1])
end

%% Plot factors (this may be very slow for matrices of size n > 1000)

fig = figure(2);
clf

marg = 0.01;
left = marg;
sc   = (1-3*marg) / (2*sum(cellfun(@(F) size(F,1), A.U)) + sum(cellfun(@(F) size(F,1), A.D)));
for l=A.tree.lvl:-1:1
    left = plot_factor(fig, A.U{l+1}, marg, sc, left);
end
left = plot_factor(fig, A.D{1}, marg, sc, left);
for l=1:A.tree.lvl
    left = plot_factor(fig, A.V{l+1}', marg, sc, left);
    left = plot_factor(fig, A.D{l+1},  marg, sc, left);
end