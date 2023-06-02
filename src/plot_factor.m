function left = plot_factor(fig, F, marg, sc, left)
    S = sparse(F);
    S(S==0) = NaN;
    width  = sc*size(S,2);
    subplot('Position', [left, marg, width, 1-marg], 'Parent', fig);
    imagesc(S)
    colormap(flipud(purples))
    left = left + width + marg;
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    pbaspect([size(S,2) size(S,1) 1])
end