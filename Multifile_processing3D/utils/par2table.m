function table_out = par2table(table,par)
    %creating a new table from the par structure
    new = struct2table(par, 'AsArray', true);
    %appending the original table
    table_out = [table; new];
end 