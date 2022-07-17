function statis = statisData(data)
%% This code computes some statistics of scRNA-seq data.
% Input:  data: a vector of scRNA-seq data
% Output: statis: [mean,noise,sk,kt,bc,ff].

data_mean = mean(data);
data_var = var(data,1);
data_cv2 = data_var / data_mean^2;
data_fano = data_var / data_mean;
data_sk = skewness(data) + 1;
data_kt = kurtosis(data);
data_bc = ((data_sk-1)^2 + 1)/data_kt;

statis = [data_mean,data_cv2,data_fano,data_sk,data_kt,data_bc];
end
