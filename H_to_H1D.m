function H_1D = H_to_H1D(H,max_grid)
if ~exist('max_grid','var')
    max_grid=max(H);
end
Ndigits=ceil(log10(max(10,abs(max_grid+1))));
Npowers=flip(cumsum(flip(Ndigits)));
Npowers=[Npowers(2:end) 0];
H_sparse_powers=H.*10.^Npowers;
H_1D=sum(H_sparse_powers,2);

end
