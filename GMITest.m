function xo = GMITest(sMix,xo)

%iNx = size(xo,1);
R = isnan(xo);
SwMV = find(sum(R,2)); % Samples with missing values
iM = sMix.ncentres;


sCovarType = sMix.covar_type;
%--------------------------------------------------------------------------
sMixR.sMixOriginal = sMix;
%--------------------------------------------------------------------------
for k = 1:length(SwMV)
    i = SwMV(k);
    %--------------------------------------------------------------------------
    iNx = sum(R(i,:)==0);
    iNy = sum(R(i,:)==1);
    %----------------------------------------------------------------------
    sMixR.iNx = iNx;
    sMixR.iNy = iNy;
    %----------------------------------------------------------------------
    indx = R(i,:)~=1;
    indy = R(i,:)==1;
    sMixR.mCentresX = sMix.centres(:,indx);
    sMixR.mCentresY = sMix.centres(:,indy);
    %--------------------------------------------------------------------------
    %---------------------- Mixture of mData features -------------------------
    sMixX = gmm(iNx, iM, sCovarType);
    sMixX.centres = sMixR.mCentresX;
    sMixX.priors = sMix.priors;
    %--------------------------------------------------------------------------
    %--------------------- Conditional Covarinces -----------------------------
    switch sCovarType
        case 'diag'
            sMixX.covars =  sMix.covars(:,indx);
            for j = 1:iM
                sMixR.mCovars(j).mCovarsXX = diag(sMix.covars(j,indx));
                sMixR.mCovars(j).mCovarsYY = diag(sMix.covars(j,indy));
                sMixR.mCovars(j).mCovarsXY = zeros(iNx,iNy);
                sMixR.mCovars(j).mCovarsYX = zeros(iNy,iNx);
            end
            
        case 'full'
            sMixX.covars =  sMix.covars(indx,indx,:);
            for j = 1:iM
                sMixR.mCovars(j).mCovarsXX = sMix.covars(indx,indx,j);
                sMixR.mCovars(j).mCovarsYY = sMix.covars(indy,indy,j);
                sMixR.mCovars(j).mCovarsXY = sMix.covars(indx,indy,j);
                sMixR.mCovars(j).mCovarsYX = sMix.covars(indy,indx,j);
            end
            
        case 'spherical'
            sMixX.covars = sMix.covars;
            for j = 1:iM
                sMixR.mCovars(j).mCovarsXX = diag(ones(1,iNx)*sMix.covars(j));
                sMixR.mCovars(j).mCovarsYY = diag(ones(1,iNy)*sMix.covars(j));
                sMixR.mCovars(j).mCovarsXY = zeros(iNx,iNy);
                sMixR.mCovars(j).mCovarsYX = zeros(iNy,iNx);
            end
    end
    sMixR.sMixX = sMixX;  
    sMixR.bFlagNorm = 0;
    %----------------------------------------------------------------------
    [imputations,~] = GMRTest(sMixR,xo(i,indx));
    xo(i,indy) = imputations;
    %----------------------------------------------------------------------
end
%--------------------------------------------------------------------------