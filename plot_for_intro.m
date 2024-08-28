n = 10;
singularvalues1 = 1./(1:n);
ks = [10, 20, 40, 80, 160, 320, 640, 1280];
T = 1000;

p = 3;
figure('Position', [100, 100, 600, 500])

variance1 = [];
ourbound = [];
boundKV = [];

for k = ks
    disp(k)
    approx = zeros(1, T);
    estvar = zeros(1, T);
    for t = 1:T
        Y = diag(singularvalues1) * [randn(n,k)];
        M = Y'*Y;
        X = triu(M, 1);
        approx(t) = nchoosek(k, p)\trace(M*X^(p-1));
        
        if p==2 && includevarvar==1
            % varvar
            avg = norm(Y, 'fro')^2 / k;
            estvar(t) = sum((vecnorm(Y).^2 - ones(1, k)*avg).^2)/(2*(k-1));
        end
    end
    variance1 = [variance1, var(approx)];
    
    ourbound = [ourbound, newbound(k, p, singularvalues1)];

    boundKV = [boundKV, 2^(12*p)*p^(6*p)*3^p*max(n^(p-2)/k^p, ...
        max(n^(1/2-1/p)/k, 1/k))];

end

loglog(ks, variance1 / spnorm(singularvalues1, 2*p)^2, '-o', 'linewidth', 2, ...
    'DisplayName','true variance')
hold on
% loglog(ks, ourbound, '--x', 'linewidth', 2, 'DisplayName', 'our bound')
loglog(ks, boundKV, ':o', 'linewidth', 2, 'DisplayName','bound from [6]')

legend
xlabel('k')
ylabel('variance')

set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(gcf, 'example1', 'epsc')


function x = newbound(k, p, sv)
    y = 1;
    for i = 0:p-1
        y = y*(k-p-i)/(k-i);
    end
    x = y*spnorm(sv,2*p)^2;

    y = p^2 / (k-p+1);
    for i = 0:p-2
        y = y * (k-p-i)/(k-i);
    end
    x = x + y*(spnorm(sv,2*p)^2 + 2*spnorm(sv,4*p));

    y = 1/factorial(2);
    for i = 0:(p-3)
        y = y*(k-p-i)/(k-i);
    end
    for i = 0:1 % from 0 to r-1
        y = y*(p-i)^2/(k-p+1+i);
    end
    x = x + y*(6*spnorm(sv,4*p) + 3*spnorm(sv,4)*spnorm(sv,4*p-4));

    for r = 3:p
        y = 1/factorial(r);
        for i = 0:p-r-1
            y = y*(k-p-i)/(k-i);
        end
        for i = 0:r-1
            y = y*(p-i)^2/(k-p+1+i);
        end

        for ell = 1:r
            if ell == 1
                coeff = 3^r;
            elseif ell==2
                coeff = 3^(r-ell) * (2^(r+1)-2);
            else
                coeff = 3^(r-ell)*ell^r;
            end
            x = x + coeff*y*spnorm(sv,4)^(ell-1)*spnorm(sv,4*p-4*(ell-1));
        end
    end
    x = x - spnorm(sv,2*p)^2;
end