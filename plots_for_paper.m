n = 100;
% singularvalues = 0.8.^(0:n-1);
% ps = [4, 6, 8];
% titlesave = "zeroeight";
% includeKVbound = 0;

singularvalues = 1./(1:n).^2;
ps = [3, 5, 7];
titlesave = "1overi";
includeKVbound = 0;

% singularvalues = ones(n, 1);
% ps = [2, 4, 6];
% titlesave = "identity";
% includeKVbound = 1;


ks = [16, 20, 40, 80, 160, 320, 500];
T = 10000;

figure('Position', [100, 100, 1600, 500])

includextrace = 0;
includevarvar = 0;

for i = 1:length(ps)
    p = ps(i);
    variance = [];
    first = [];
    second = [];
    bound = [];
    bound2 = [];
    boundKV = [];
    varvar = [];
    varxtrace = [];

    for k = ks
        disp(k)

        if (k < 100)
            T = 100000;
        else
            T = 10000;
        end
        approx = zeros(1, T);
        estvar = zeros(1, T);
        for t = 1:T
            Y = diag(singularvalues) * [randn(n,k)];
            M = Y'*Y;
            X = triu(M, 1);
            approx(t) = nchoosek(k, p)\trace(M*X^(p-1));
            
            if p==2 && includevarvar==1
                % varvar
                avg = norm(Y, 'fro')^2 / k;
                estvar(t) = sum((vecnorm(Y).^2 - ones(1, k)*avg).^2)/(2*(k-1));
            end
        end
        variance = [variance, var(approx)];

        if p==2 && includevarvar==1
            varvar = [varvar, var(estvar)];
        end
        
        first = [first, 2*p^2*(sum(singularvalues.^(4*p)))/k];
        second = [second, second_order_estimate(k, p, singularvalues)];
        bound = [bound, newbound(k, p, singularvalues)];
        boundKV = [boundKV, 2^(12*p)*p^(6*p)*3^p*max(n^(p-2)/k^p, ...
            max(n^(1/2-1/p)/k, 1/k))* spnorm(singularvalues, 2*p)^2];
        
    end

    subplot(1, 3, i)
    plot_our_bounds(p, ks, variance, first, second, bound, boundKV,...
        includeKVbound, varvar, includevarvar, varxtrace, includextrace)
    hold off

end

set(findall(gcf,'-property','FontSize'),'FontSize',14)
saveas(gcf, titlesave, 'epsc')

function x = second_order_estimate(k, p, sv)
    x = 2*p^2*spnorm(sv, 4*p)/k;
    y = spnorm(sv, 4*p) + 3/2*spnorm(sv,4)*spnorm(sv, 4*p-4)-1/2*spnorm(sv,2*p)^2;
    x = x + y*p^2*(p-1)^2/(k^2);
end

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

function x = c(r, ell)
    if ell == 1
        x = 3^r;
    elseif ell == 2
        x = 3^(r-2)*(2^(r+1)-1);
    else
        x = 3^(r-ell) * ell^r;
    end
end



function plot_our_bounds(p, ks, variance, first, second, bound1, ...
    boundKV, includeKV, varvar, includevarvar, varxtrace, includextrace)
    loglog(ks, variance, '-s', 'linewidth', 2, 'DisplayName','variance')
    hold on
    loglog(ks, first, '-o', 'linewidth', 2, 'DisplayName','first-order est.')
    loglog(ks, second, '-x', 'linewidth', 2, 'DisplayName','second-order est.')
    loglog(ks, bound1, ':p', 'linewidth', 2, 'DisplayName','upper bound')
    if includeKV == 1
        loglog(ks, boundKV, ':d', 'linewidth', 2, 'DisplayName','bound from [6]')
    end
    if includevarvar == 1 && p==2
        loglog(ks, varvar, '-o', 'linewidth', 2, 'DisplayName', 'varvar')
    end
    if includextrace == 1
        loglog(ks, varxtrace, '-o', 'linewidth', 2, 'DisplayName', 'XTrace')
    end

    legend
    xlabel('k')
    ylabel('variance')
    title(['p=', num2str(p)]);
end









