%clear all;
% This code depends on the NetLab toolbox. It is free available from:
% http://www.aston.ac.uk/eas/research/groups/ncrg/resources/netlab/downloads/
%--------------------------------------------------------------------------
%---------------------- Example one ---------------------------------------
ndata = 500;
data = randn(ndata, 2);
prior = [0.3 0.5 0.2];
% Mixture model swaps clusters 1 and 3
datap = [0.2 0.5 0.3];
datac = [0 2; 0 0; 2 3.5];
datacov = repmat(eye(2), [1 1 3]);
data1 = data(1:prior(1)*ndata,:);
data2 = data(prior(1)*ndata+1:(prior(2)+prior(1))*ndata, :);
data3 = data((prior(1)+prior(2))*ndata +1:ndata, :);

% First cluster has axis aligned variance and centre (2, 3.5)
data1(:, 1) = data1(:, 1)*0.1 + 2.0;
data1(:, 2) = data1(:, 2)*0.8 + 3.5;
datacov(:, :, 3) = [0.1*0.1 0; 0 0.8*0.8];

% Second cluster has variance axes rotated by 30 degrees and centre (0, 0)
rotn = [cos(pi/6) -sin(pi/6); sin(pi/6) cos(pi/6)];
data2(:,1) = data2(:, 1)*0.2;
data2 = data2*rotn;
datacov(:, :, 2) = rotn' * [0.04 0; 0 1] * rotn;

% Third cluster is at (0,2)
data3(:, 2) = data3(:, 2)*0.1;
data3 = data3 + repmat([0 2], prior(3)*ndata, 1);

% Put the dataset together again
data = [data1; data2; data3];
figure;
plot(data(:,1),data(:,2),'.');
%--------------------------------------------------------------------------
indNAN = randperm(size(data,1));
indVar = (rand(1,size(data,1))>0.5) + 1;
cont = 0;
for MValPer = 0.05:0.05:0.5 %Percentage of missing values
    cont = cont + 1;
    for M = 2:10 % Number of Gaussians in mixture
        
        %----------------------------------------------------------------------
        indtem = indNAN(1:ceil(ndata*MValPer));
        dataNaN = data;
        ind1 = indVar(indtem) == 1;
        dataNaN(ind1,1) = NaN;
        ind2 = indVar(indtem) == 2;
        dataNaN(ind2,2) = NaN;
        %----------------------------------------------------------------------
        [~,mDataIu,~] = GMI(dataNaN,M,'full');
        %----------------------------------------------------------------------
        %--------------- MSE estimation ---------------------------------------
        tem = sum((data(ind1,1) - mDataIu(ind1,1)).^2);
        tem = tem + sum((data(ind2,2) - mDataIu(ind2,2)).^2);
        Error(cont,M-1) = sqrt(tem/length(indtem));
    end
end
%--------------------------------------------------------------------------
MValPer = 0.05:0.05:0.5;
M = 2:10;
figure(),surf(MValPer,M,Error');
colorbar;
