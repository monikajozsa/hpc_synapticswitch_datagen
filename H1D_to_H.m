function H = H1D_to_H(H_1D,Nspecies)

H=zeros(length(H_1D),Nspecies);
digits_all=numel(num2str(max(floor(H_1D))));
digits_1=numel(num2str(max(floor(H_1D))))-numel(num2str(min(floor(H_1D))))+1;
H(:,1)=floor(H_1D/10^(digits_all-digits_1));
H_1D=H_1D-H(:,1)*10^(digits_all-digits_1);
digits_2=numel(num2str(max(floor(H_1D))))-numel(num2str(min(floor(H_1D))))+1;
H(:,2)=floor(H_1D/10^(digits_all-digits_1-digits_2));
H(:,3)=H_1D-H(:,2)*10^(digits_all-digits_1-digits_2);

end