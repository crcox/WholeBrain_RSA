function c = drop_index_and_adjust(full_index, val_to_drop)
    c = full_index(val_to_drop~=full_index);
    z = c > val_to_drop;
    c(z) = c(z) - 1;
end