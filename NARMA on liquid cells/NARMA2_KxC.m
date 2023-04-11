function [err_train, err_test] = NARMA2_KxC(Ntrain, Ntest, K, C, flag, pmax)
Nn = 1000;

uk = rand(Nn+1,1)/2;
yk = zeros(Nn+1,1);
for i = 3:Nn+1
    yk(i) = 0.4*yk(i-1)+0.4*yk(i-1)*yk(i-2)+0.6*(uk(i-1))^3+0.1;
end


dt = 0.1; %ms
D = 1/dt;
J = 1;
pw = 25;  %ms
res_t = 10;  %ms
%detector_int = 0.1;  %ms
%pmax = 0.08;%mW
poffset = 0.0; %mW

segm = round(ceil(Nn/J));



L = round(segm*(pw+res_t)*D);

% K = 24;
% C = 24;


inp_w = rand(K,1);

inp = uk*inp_w';
inp = (inp-min(inp,[],'all'))/max(inp,[],'all');

inpm = zeros(K,L);

for k = 1:K
    for i = 1:segm
        for j = 1:J
            inpm(k,round((i-1)*(pw+res_t)*D+(j-1)*(pw*D/J))+1:round((i-1)*(pw+res_t)*D+j*pw*D/J)) = sqrt(pmax*inp((i-1)*J+j,k)+poffset);
        end
    end
end



out_M = zeros(L,K);
hvec = 0.8*ones((K-1)*C/2+mod(C,2)/2);
hmatrix = zeros(L,(K-1)*C/2+mod(C,2)/2);
Pmatrix = zeros(L,(K-1)*C/2+mod(C,2)/2);
theta = pi/3*ones((K-1)*C/2+mod(C,2)/2,1);

for i = 1:L
    inp_ = inpm(:,i);
    [temp, hnvec,P] = RC_KxC(inp_,theta,hvec,dt,K,C,flag);
    out_M(i,1:K) = (abs(temp(1:K))).^2;
    hvec = hnvec;
    hmatrix(i,:) = hvec(:);
    Pmatrix(i,:) = P(:);
    %disp(i);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Training %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


train_start = 20;
out_train = out_M(train_start*(pw+res_t)*D+1:(train_start+Ntrain)*(pw+res_t)*D,:);
input_train = zeros(Ntrain,K*pw*D/J);

segm = round(ceil(Ntrain/J));
for i = 1:segm
    for j = 1:J
        for k = 1:K
            input_train((i-1)*J+j,(k-1)*pw*D/J+1:k*pw*D/J) = out_train(round((i-1)*(pw+res_t)*D+(j-1)*(pw*D/J))+1:round((i-1)*(pw+res_t)*D+j*pw*D/J),k);
        end
    end
end

ytrain = yk(train_start+2:train_start+Ntrain+1);

kc = 0.000001;

w = ridge(ytrain,input_train,kc,0);

yo = w(1) + input_train*w(2:end);

err_train = immse(ytrain,yo)/immse(ytrain,zeros(Ntrain,1));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Testing %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


test_start = 700;
out_test = out_M(test_start*(pw+res_t)*D+1:(test_start+Ntest)*(pw+res_t)*D,:);
input_test = zeros(Ntest,2*pw*D/J);

segm = round(ceil(Ntest/J));
for i = 1:segm
    for j = 1:J
        for k = 1:K
            input_test((i-1)*J+j,(k-1)*pw*D/J+1:k*pw*D/J) = out_test(round((i-1)*(pw+res_t)*D+(j-1)*(pw*D/J))+1:round((i-1)*(pw+res_t)*D+j*pw*D/J),k);
        end
    end
end

ytest = yk(test_start+2:test_start+Ntest+1);

yo = w(1) + input_test*w(2:end);

err_test = immse(ytest,yo)/immse(ytest,zeros(Ntest,1));

    
    
end