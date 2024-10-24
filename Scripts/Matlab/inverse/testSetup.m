function L = testSetup(X, Y, true_val)

beta = 10^19;
L = -1/2*beta*sum(abs((true_val-Y)).^2,'all')/length(true_val);
%Lprev = -sum(abs((Y_true-Y_last)*10^7).^2,'all')*10;
%Lmean = -sum(abs((Y_true-(mean(F_array) + 2*10^(-8)))*10^7).^2,'all')*10;