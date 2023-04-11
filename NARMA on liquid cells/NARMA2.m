function [err_train, err_test] = NARMA2(Ntrain, Ntest)
Nn = 2000;

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
pmax = 0.15;%mW
poffset = 0.0; %mW

segm = round(ceil(Nn/J));


L = round(segm*(pw+res_t)*D);

K = 2;


%inp_w = rand(K,1);
inp_w = zeros(K,1);
inp_w(1) = 1;

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


temp = zeros(2,1);
out_M = zeros(L,2);
hvec = 0.8;
hmatrix = zeros(1,L);
Pmatrix = zeros(1,L);
omvec = zeros(2,L);
theta = zeros(1,2);
theta(1) = 0.6155;
%theta(1) = pi/4;
theta(2) = pi/3;

for i = 1:L
    [temp(1), temp(2),hnvec,P,o1,o2] = RC_1coup(inpm(1,i),inpm(2,i),theta,hvec,dt);
    out_M(i,1:2) = (abs(temp(1:2))).^2;
    hvec = hnvec;
    hmatrix(i) = hvec;
    Pmatrix(i) = P;
    omvec(1,i) = abs(o1)^2;
    omvec(2,i) = abs(o2)^2;
    %disp(i);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Training %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ntrain = 400;
train_start = 20;
out_train = out_M(train_start*(pw+res_t)*D+1:(train_start+Ntrain)*(pw+res_t)*D,:);
input_train = zeros(Ntrain,2*pw*D/J);

segm = round(ceil(Ntrain/J));
for i = 1:segm
    for j = 1:J
        for k = 1:2
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

%Ntest = 100;
test_start = 1500;
out_test = out_M(test_start*(pw+res_t)*D+1:(test_start+Ntest)*(pw+res_t)*D,:);
input_test = zeros(Ntest,2*pw*D/J);

segm = round(ceil(Ntest/J));
for i = 1:segm
    for j = 1:J
        for k = 1:2
            input_test((i-1)*J+j,(k-1)*pw*D/J+1:k*pw*D/J) = out_test(round((i-1)*(pw+res_t)*D+(j-1)*(pw*D/J))+1:round((i-1)*(pw+res_t)*D+j*pw*D/J),k);
        end
    end
end

ytest = yk(test_start+2:test_start+Ntest+1);

yo = w(1) + input_test*w(2:end);

err_test = immse(ytest,yo)/immse(ytest,zeros(Ntest,1));


    
    
end