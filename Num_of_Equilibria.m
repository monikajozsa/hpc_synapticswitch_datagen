function [Neq ,Jsubs]= Num_of_Equilibria(rate_constants,xbar)
    Knum_E=rate_constants.Knum_E;
    Kdenum_E=rate_constants.Kdenum_E;
    Knum_I=rate_constants.Knum_I;
    Kdenum_I=rate_constants.Kdenum_I;
    beta=rate_constants.betas;
    if size(rate_constants.Knum_E,1)==2
        syms x1 x2
        eqns = [
            sum((Knum_E(2,1)*x2./xbar(2)*xbar(1))./(Kdenum_E(2,1)+x2/xbar(2)) + (Knum_I(2,1)*x2./xbar(2))./(Kdenum_I(2,1)+x2/xbar(2)))==beta(1)*x1, ...
            sum((Knum_E(1,2)*x1./xbar(1)*xbar(2))./(Kdenum_E(1,2)+x1/xbar(1)) + (Knum_I(1,2)*x1./xbar(1))./(Kdenum_I(1,2)+x1/xbar(1)))==beta(2)*x2];
        S = solve(eqns,[x1 x2]);
        Neq=[length(find(S.x1>0)) length(find(S.x2>0))];
        
        diff_funct = [(Knum_E(2,1)*x2./xbar(2)*xbar(1))./(Kdenum_E(2,1)+x2/xbar(2))-beta(1)*x1, ...
            (Knum_E(1,2)*x1./xbar(1)*xbar(2))./(Kdenum_E(1,2)+x1/xbar(1))-beta(2)*x2];
        if all(Neq==1)
            J=jacobian(diff_funct);
            Jsubs=subs(J,[x1 x2],xbar'); % we calculate this for stability check eig(Jsubs)
        else
            Jsubs=0;
        end
    end    
    if size(rate_constants.Knum_E,1)==3
        syms x1 x2 x3
        eqns = [
            sum((Knum_E([2 3],1).*[x2;x3]./xbar([2 3])*xbar(1))./(Kdenum_E([2 3],1)+[x2;x3]./xbar([2 3])) + (Knum_I([2 3],1)*xbar(1))./(Kdenum_I([2 3],1)+[x2;x3]./xbar([2 3])))==beta(1)*x1, ...
            sum((Knum_E([1 3],2).*[x1;x3]./xbar([1 3])*xbar(2))./(Kdenum_E([1 3],2)+[x1;x3]./xbar([1 3])) + (Knum_I([1 3],2)*xbar(2))./(Kdenum_I([1 3],2)+[x1;x3]./xbar([1 3])))==beta(2)*x2, ...
            sum((Knum_E([1 2],3).*[x1;x2]./xbar([1 2])*xbar(3))./(Kdenum_E([1 2],3)+[x1;x2]./xbar([1 2])) + (Knum_I([1 2],3)*xbar(3))./(Kdenum_I([1 2],3)+[x1;x2]./xbar([1 2])))==beta(3)*x3];
        S = solve(eqns,[x1 x2 x3],'ReturnConditions',true);
        Neq=[length(find(S.x1>0)) length(find(S.x2>0)) length(find(S.x3>0))];
        
        diff_funct=[sum((Knum_E([2 3],1).*[x2;x3]./xbar([2 3])*xbar(1))./(Kdenum_E([2 3],1)+[x2;x3]./xbar([2 3])) + (Knum_I([2 3],1)*xbar(1))./(Kdenum_I([2 3],1)+[x2;x3]./xbar([2 3])))-beta(1)*x1, ...
            sum((Knum_E([1 3],2).*[x1;x3]./xbar([1 3])*xbar(2))./(Kdenum_E([1 3],2)+[x1;x3]./xbar([1 3])) + (Knum_I([1 3],2)*xbar(2))./(Kdenum_I([1 3],2)+[x1;x3]./xbar([1 3])))-beta(2)*x2, ...
            sum((Knum_E([1 2],3).*[x1;x2]./xbar([1 2])*xbar(3))./(Kdenum_E([1 2],3)+[x1;x2]./xbar([1 2])) + (Knum_I([1 2],3)*xbar(3))./(Kdenum_I([1 2],3)+[x1;x2]./xbar([1 2])))-beta(3)*x3];
        if all(Neq==1)
            J=jacobian(diff_funct);
            Jsubs=subs(J,[x1 x2 x3],xbar'); % we calculate this for stability check eig(Jsubs)
        else
            Jsubs=0;
        end
    end
end