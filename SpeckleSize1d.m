% function to calculate speckle size in 1d
function [speckle_FWHM,V] = SpeckleSize1d(speckle_1d)

% input: speckle_1d whose first dimension is the speckle pattern in 1d and
% the second dimension is the number of speckle patterns

% output: FWHM of 1d speckle with Gaussian fitting

N = size(speckle_1d,2);
E1=speckle_1d(:,1);                  %creates an array with the values of column 1
s1=size(xcov(E1));          %finds the size of the autocovariance array
D1=zeros(s1);               %creates an empty array of size s1
D1=double(D1);              %typecasts the values as doubles
for j=1:N
C1=speckle_1d(:,j);
D1=imadd(D1,xcov(C1,'coeff'));      %sums the xcov arrays for all rows into D1
end
V=(D1/max(D1));             %the finished vertical product, V, is normalized

helper2 = 1:size(V);

gauss1 = fittype('gauss1');                                      %sets up the Gaussian curve fitting
excludeLowV = excludedata(helper2',V,'range',[.2,1]);
optionsV = fitoptions(gauss1);
optionsV.Exclude = excludeLowV;

[VFit, VFitStats] = fit(helper2',V,gauss1, optionsV);

%The SpeckleSize code technically does not find the FWHM and 1/e2 widths of
%the Gaussian fits, but rather the widths at which the Gaussians fall to 
%values of .5 and 1/e^2.  This provides a better idea of the speckles� size
%even when the Gaussian�s amplitude is not unity, as expected for a perfect
%fit of the normalized data, but can produce unexpected trends...


%FWHM values (Full width when the fit = .5)  
speckle_FWHM = (2*(VFit.c1)*sqrt(-log(.5/(VFit.a1))));


