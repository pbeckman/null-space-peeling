function idx_out = combine_idx(idx_in)
    idx_out = cell(length(idx_in)/2, 1);
    for b=1:length(idx_out)
        i1 = idx_in{2*(b-1) + 1};
        i2 = idx_in{2*(b-1) + 2};
        idx_out{b} = [i1(1), i2(2)];
    end
end