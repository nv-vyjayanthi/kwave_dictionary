% do 1D iterative reconstruction
% use FISTA to minimize
% 1/2 * || A*PREC - Y ||_2^2 + la1 * ||PREC||_21 + la2/2 * ||PREC||_2^2
% 1/2 * || A*PREC - Y ||_2^2 + la1 * ||PREC||_1 + la2/2 * ||PREC||_2^2 
%
% Y     ... data 
% Prec  ... recovered pressure
% A     ... PSF 
% L     ... Lipschitz constant
% Niter ... iteration number in fista
% la1   ...
% mu    ... 
% regularization   ... type of regularization term
% ... 'fista'   R = la1 * ||PREC||_1 + la2/2 * ||PREC||_2^2 ... sparse
% ... 'b-fista' R = la1 * ||PREC||_21 + la2/2 * ||PREC||_2^2 ... blocksparse  
% 
% last modified : 25.02.2016
%

function Prec = reconiter1D(Y,Prec, A, L, Niter, la1, la2, alg )

[n,m] = size(Y);
switch alg
    case 'fista'
        figure
        for iter =1:Niter
            Yi = idct(dct(Prec).*A);
            Prec = Prec - (2/L)*idct( dct(Yi-Y) .* A);
            Prec = max(0,Prec);
            Prec = proxnet(Prec,la1,la2);
        end
       
        
        case 'Bfista'
        Q=Prec;
        t=1;
        for iter =1:Niter
            tprev = t;
            Pprev = Prec; 
            Yi = idct(dct(Q).*A);
            Q =  Q - (2/L)*idct( dct(Yi-Y) .* A);
            %Q = max(0,Q);
            Q = netblock(Q,m, la1/L, la2/L);
            Prec = Q;
            t = (1+sqrt(1+4*t^2))/2;
            Q = Prec + (tprev-1)/t*(Prec-Pprev);
        end
end



% elastic net thresholding for fista
function X = proxnet(X,la1,la2)
X   = X/(1+la2);
la1 = la1/(1+la2);
X   = sign(X).*max(abs(X)-la1,0);
return;


% block elastic net for block fista
function X = netblock(X,m,la1,la2)
X   = X/(1+la2);
la1 = la1/(1+la2);

% rq = sum(X.^2,2);  rq(~(rq==0)) = max( ( 1 - la1^2./rq(~(rq==0))),0);
% above version does something different; similar to james-stein thresholding 

rq = sqrt(sum(X.^2,2));  rq(~(rq==0)) = max( ( 1 - la1./rq(~(rq==0))),0);

X = X .* ( rq * ones(1,m) );
return;





