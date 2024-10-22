function data_resized = getRandomCurves(covfunc,arg1,arg2,X)
% covfunc: string, "matern" or "squared_exponential"
% arg1 and arg2 must be strings


filename = sprintf('../../../Data/gaussian_process_realisations/curve_%s_p_%s_ell_%s.mat',covfunc,arg1,arg2);
load(filename,'data','x');


minx = min(X);
maxx = max(X);

x = x/max(x)*(maxx-minx) + minx;

[m,~] = size(data);
data_resized = zeros(m,length(X));

parfor i = 1:m

    data_resized(i,:) = interp1(x,data(i,:),X);

end

data_resized = data_resized*10^(-8);