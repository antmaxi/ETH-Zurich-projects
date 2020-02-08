function C = chi2_cost(s1,s2, n_out, tangent_flag)

    n1 = size(s1,1);
    n2 = size(s2,1);
    if n1~=n2
        disp("different sizes");
    end
    C = zeros(n1, n2);
    for i = 1:n1-n_out
        S1 = s1(i, :,:);
        for j = 1:n2-n_out
            if (tangent_flag)
                S2 = s2(j,:,:);
                minus2 = (S1 - S2).^2;
                plus = S1 + S2;
                C1 = 0.5*sum(minus2(plus ~= 0)./plus(plus ~= 0));
                S2 = fliplr(s2(j,:,:));
                minus2 = (S1 - S2).^2;
                plus = S1 + S2;
                C2 = 0.5*sum(minus2(plus ~= 0)./plus(plus ~= 0));
                C(i,j) = min(C1,C2);
            else
                S2 = s2(j,:,:);
                minus2 = (S1 - S2).^2;
                plus = S1 + S2;
                C(i,j) = 0.5*sum(minus2(plus ~= 0)./plus(plus ~= 0));
            end
        end
    end
    out_cost = 55; % ouliers cost if used
    if n_out
        C(:,n2-n_out+1:n2) = out_cost; 
        C(n1-n_out+1:n1,:) = out_cost;
    end
    %disp(C);
end