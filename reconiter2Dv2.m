% do 2D iterative reconstruction
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
% ... 'rb-fista' R = la1 * ||PREC||_21 + la2/2 * ||PREC||_2^2 
% ... uses JOINT sparsity in specke patterns and rows of image  
% 
% last modified : 25.02.2016
%

% =================================
function Prec = reconiter2D( Y,Prec, A, L, Niter, la1, la2, regularization, plotit )

[Nx,~,N] = size(Y); % notice: size Y = [Nx,Nx,N]

Q=Prec;
t=1;

if plotit; figure; end 
for iter =1:Niter
    if (~plotit && mod(iter,20) == 1);  disp([' iter ' num2str(iter), ' von ', num2str(Niter) ]); end;
    tprev = t;
    Pprev = Prec;
    Yi = convolve3(Q,Nx,N,A);
    Q =  Q - (2/L)*convolve3(Yi-Y,Nx,N,A);
    % type of regularization
    % Q = max(0,Q);  % positivity constraint
    switch regularization
        case 'fista'
            Q = proxnet(Q,la1/L, la2/L);
        case 'b-fista'
            Q = netblock(Q,Nx,N,la1/L, la2/L );
        case 'rb-fista'
            Q = rownetblock(Q,Nx,N, la1/L, la2/L );
    end
    
    Prec = Q;
    t = (1+sqrt(1+4*t^2))/2;
    Q = Prec + (tprev-1)/t*(Prec-Pprev);
    
    if plotit
        subplot(121), imagesc(sum(Prec,3)); colormap(hot);
        title([regularization, ' (mean): iter ' num2str(iter), ' von ', num2str(Niter) ]); pause(0),drawnow    
        subplot(122), imagesc(sqrt(var(Prec,0,3))); colormap(hot);
        title([regularization, ' (var): iter  ' num2str(iter), ' von ', num2str(Niter) ]); pause(0),drawnow
    end
end
% =================================        


% =================================
% convolution with A
% applied to 1st+2nd component  
function F = convolve3(F,Nx,N,A)
f = zeros(Nx);    
     for j=1:N
        f(:,:) = F(:,:,j); 
        f = convolve(f,A);
        F(:,:,j) = f;
    end
return;
%  =================================


% =================================
% block elastic net thresholding of X
% each block is a vector X(i,j,:)
function X = netblock(X,Nx,N,al1,al2)
Y = reshape(X,Nx^2,N);  
Y  = Y/(1+al2);
al1 = al1/(1+al2);
% rq = sum(Y.^2,2)/N; rqp=rq;
% rq(~(rq==0)) = max( ( 1 - al^2./rq(~(rq==0))),0);
rq = sqrt(sum(Y.^2,2)); rqp=rq;
rq(~(rq==0)) = max( ( 1 - al1./rq(~(rq==0))),0);

Y = Y .* ( rq * ones(1,N) );
X = reshape(Y,Nx,Nx,N); 
%if plotit; subplot(133), imagesc(rqp); title('squared norms of blocks'); pause(0),drawnow; end
return;
% =================================

% =================================

% =================================
% block1d elastic net thresholding of X 
% each block is a vector X(i,:,:)
% this requires that phantom is quite
% constant along 2d dimension
function X = rownetblock(X,Nx,N,al,c)
X = permute(X,[2 1 3]);
Y = reshape(X,Nx,Nx*N);
Y  = Y/(1+c*al);
al = al/(1+c*al);
rq = sum(Y.^2,2)/(Nx*N); rqp=rq;
rq(~(rq==0)) = max( ( 1 - al^2./rq(~(rq==0))),0);
Y = Y .* ( rq * ones(1,Nx*N) );
X = reshape(Y,Nx,Nx,N);
X = permute(X,[2 1 3]);
% if plotit; subplot(133), plot(1:Nx,rqp); title('squared norms of blocks'); pause(0),drawnow; end
return;
% =================================


function X = proxnet(X,al1,al2)
al1 =  al1/100;
X   = X/(1+al2);
al1 = al1/(1+al2);
X   = sign(X).*max(abs(X)-al1,0);
%if plotit; subplot(133), plot(X);  end
return;


