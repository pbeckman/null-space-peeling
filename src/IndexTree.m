classdef IndexTree
    properties
        lvl
        idx
    end
    methods
        function tree = IndexTree(n, lvl)
            tree.lvl = lvl;
            tree.idx = cell(lvl+1,1);
            tree.idx{1} = {[1 n]};
            for l=1:lvl
                idxpar = tree.idx{l};
                idxl   = cell(2^l,1);
                for b=1:2:2^l
                    bpar = idxpar{(b-1)/2+1};
                    idxl{b}   = [bpar(1), bpar(1) + floor((bpar(2)-bpar(1))/2)];
                    idxl{b+1} = [bpar(1) + floor((bpar(2)-bpar(1))/2) + 1, bpar(2)];
                end
                tree.idx{l+1} = idxl;
            end
        end
        function a = ancestor(~, al, b, cl)
            % index of ancestor at level al for child box b at level cl
            a = fix((b-1) / 2^(cl-al)) + 1;
        end
    end
end