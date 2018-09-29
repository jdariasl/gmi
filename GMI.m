% function [mDataI,mDataIu,sMix] = GMI(mData,iM)
%
% Imputates missing variables using a Gaussian Mixture Model (GMM). It is based on the NetLab toolbox.
%
% Inputs:
%
%        mData: represents the data whose expectation is maximized, with
%               each row corresponding to a sample. Missing values must be
%               represented by NaN.
%
%
%        iM: Number of Components in Mixture
%        
%        sCovarType: The mixture model type defines the covariance structure 
%                   of each component  Gaussian:
%                   'spherical' = single variance parameter for each component: 
%                                 stored as a vector
%               	'diag' = diagonal matrix for each component: stored as 
%                            rows of a matrix. Default option.
%                	'full' = full matrix for each component: stored as 3d array
%
% Ouputs:
%
%       sMix is  a MatLab structure containing the GMI model.
%
%       mDataI es the data with the missing values inputed by the GMM
%       model. mDataIu is the unnormalized version of mDataI.
%
% See: GMItest
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Copyright (c) 2016, Julián David Arias Londoño All rights reserved. %%%
%%%%%%%%%%%%%%%%%%% Department of Systems Engineering %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Universidad de Antioquia, Colombia %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  2016  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [mDataI,mDataIu,sMix,media,desvia] = GMI(mData,iM,sCovarType)

if nargin < 3, sCovarType = 'diag'; end

%--------------------------------------------------------------------------
[iNx,iNv] = size(mData);
%--------------------------------------------------------------------------
%sMixR.iNx = iNx;
%--------------------------------------------------------------------------
R = isnan(mData);
media = zeros(1,iNv);
desvia = zeros(1,iNv);
for i = 1:iNv
    media(i) = mean(mData(R(:,i)~=1,i));
    desvia(i) = std(mData(R(:,i)~=1,i));
end
mDataN = NormalizaMissingValues(mData,media,desvia,R);
%--------------------------------------------------------------------------
%sMixR.vMeanX = media;
%sMixR.vSigmaX = desvia;
%--------------------- Initialization ------------------------------------- 
mDataNt = mDataN;
mDataNt(R) = 0;%Initial imputation using mean
sMix = gmm(iNv, iM, sCovarType);
vOptions = foptions;
vOptions(14) = 20; % maximum number of iterations
vOptions(1) = -1; % Switch off all messages, including warning
sMix = gmminit(sMix, mDataNt, vOptions);
%---------------------- Adjusting -----------------------------------------
vOptions(14) = 40; % maximum number of iterations
vOptions(5) = 1; % Reset covariance matrices when singular
vOptions(1) = 1; % Switch messages on; 
try
    sMix = gmmem_imputation(sMix, mDataN, vOptions);
    mDataI = GMITest(sMix,mDataN);
catch
    disp('Means instead of EM');
    mDataI = mDataNt;
end
mDataIu = (mDataI.*repmat(desvia,iNx,1)) + (repmat(media,iNx,1));
mDataIu(~R) = mData(~R);