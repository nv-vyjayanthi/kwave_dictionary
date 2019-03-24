% do iterative recon sparsity of  prec

function prec = deconviter(y,prec, A, L, Niter, la, alg )

switch alg
    case 'ista'
        figure
        for iter =1:Niter
            yi = idct(dct(prec).*A);
            prec = prec - (2/L)*idct( dct(yi-y) .* A);
            prec = max(0, soft(prec, la/L));
        end
        
    case 'fista'
        q=prec;
        t=1;
        for iter =1:Niter
            tprev = t;
            prev = prec; 
            yi = idct(dct(q).*A);
            q = q - (2/L)*idct( dct(yi-y) .* A);
            q = max( 0, soft(q, la/L));
            prec = q;
            t = (1+sqrt(1+4*t^2))/2;
            q = prec + (tprev-1)/t*(prec-prev);
        end
        
        case 'fastnet'
        q=prec;
        t=1;
        for iter =1:Niter
            tprev = t;
            prev = prec; 
            yi = idct(dct(q).*A);
            q =  q - (2/L)*idct( dct(yi-y) .* A);
            q = max(0, proxnet(q, la/L, 0.05));
            prec = q;
            t = (1+sqrt(1+4*t^2))/2;
            q = prec + (tprev-1)/t*(prec-prev);
        end
end



% soft thresholding 
function x = soft(x,alpha)
x = sign(x).*max(abs(x)-alpha,0);
return;

% elastic net thresholding 
function x = proxnet(x,alpha,c)
x = x/(1+c*alpha);
alpha=alpha/(1+c*alpha);
x = soft(x,alpha);
return;






