function [err_train, err_test, Conf] = RC_implementation(K, C, arr_req)
    
    tic
    [images_train, labels_train, count_train, images_test, labels_test, count_test] =  get_data(arr_req); %%% Get training and Testing data
    
    
    %% PCA
    N_train = sum(count_train);
    N_test = sum(count_test);

    M = size(arr_req,2);
    count_train_cum = cumsum(count_train);
    y_train = zeros(N_train,M);
    for i = 1:M
        y_train(count_train_cum(i)+1:count_train_cum(i+1),i) = 1;
    end

    count_test_cum = cumsum(count_test);
    y_test = zeros(N_test,M);
    for i = 1:M
        y_test(count_test_cum(i)+1:count_test_cum(i+1),i) = 1;
    end

    N = N_train;

    inpm = images_train;
    inp = zeros(N,784);
    for i=1:N
        temp1 = inpm(:,:,i);
        temp1 = temp1';
        temp = reshape(temp1,[1,784]);
        inp(i,:) = temp;
    end

    XTrain = inp;
    [coeff,scoreTrain,~,~,explained,mu] = pca(XTrain);

    idx = find(cumsum(explained)>95,1);
    %idx = round(idx,-2);
    %idx = floor(idx/8)*8;
    %idx = 27;
    disp(['# PC = ', num2str(idx)])
    inp = scoreTrain(:,1:idx);
    for i = 1:N
        inp(i,:) = inp(i,:) - min(inp(i,:));
        inp(i,:) = inp(i,:)/max(inp(i,:));
    end
    inp_train = inp;
    %y_train = labels_train;


    N = N_test;
    start = 1;

    inpm = images_test;
    inp = zeros(N,784);
    for i=1:N
        temp1 = inpm(:,:,i);
        temp1 = temp1';
        temp = reshape(temp1,[1,784]);
        inp(i,:) = temp;
    end
    XTest = inp;
    scoreTest = (XTest-mu)*coeff;
    inp = scoreTest(:,1:idx);
    for i = 1:N
        inp(i,:) = inp(i,:) - min(inp(i,:));
        inp(i,:) = inp(i,:)/max(inp(i,:));
    end
    inp_test = inp;
    %y_test = labels_test;

    %% RC variables
    dt = 0.1; %ms
    D = 1/dt;
    J = 25;
    pw = 25;  %ms
    res_t = 10;  %ms
    detector_int = 0.1;  %ms
    pmax = 0.08;%mW
    poffset = 0.0; %mW
    theta = pi/4*ones((K-1)*C/2+mod(C,2)/2,1);
    
    %% Data preparation for Train-RC
    
    N = N_train;
    segm = round(ceil(idx/J));
    inp =zeros(N,J*segm);
    inp(1:N,1:idx) = inp_train;

    L = round(idx*(pw+res_t)*D/J);
    L = K*round(ceil(L/K));
    inpm = zeros(N, L);
    for j = 1:N
        for i = 1:segm
            for k = 1:J
                inpm(j,round((i-1)*(pw+res_t)*D+(k-1)*(pw*D/J))+1:round((i-1)*(pw+res_t)*D+k*pw*D/J)) = sqrt(pmax*inp(j,(i-1)*J+k)+poffset);
            end
        end
    end

    %% Train-RC!
    out_M = zeros(N,L/K,K);
    hvec = 0.8*ones((K-1)*C/2+mod(C,2)/2);
    hmatrix = zeros(N,L/K,(K-1)*C/2+mod(C,2)/2);
    Pmatrix = zeros(N,L/K,(K-1)*C/2+mod(C,2)/2);

    for j = 1:N
        for i = 1:L/K
            inp_ = zeros(K,1);
            for k = 1:K
                inp_(k) = inpm(j,i+(k-1)*L/K);
            end
            [temp, hnvec,P] = RC_KxC(inp_,theta,hvec,dt,K,C);
            out_M(j,i,1:K) = (abs(temp(1:K))).^2;
            hvec = hnvec;
            hmatrix(j,i,:) = hvec(:);
            Pmatrix(j,i,:) = P(:);
            %disp(i);
        end
    end
    
    %% Gather training output
    %disp(L*dt/(detector_int))
    out_Mo = reshape(out_M,[N round(L*dt/(detector_int))]);
    temp = out_Mo;
    Ki = 10;
    out_size = round(floor(L*dt/(detector_int)/Ki));
    out_Mo = zeros(N,out_size);
    for j = 1:N
        for i = 1:out_size
            out_Mo(j,i) = mean(temp(j,(i-1)*Ki+1:i*Ki));
        end 
    end 
    out_Mo = cat(2,out_Mo,inp);
    
    %% Train!
    Cw = size(out_Mo,2)+1;
    y = y_train;
    N = N_train;
    kc = 0.001;
    w_matrix = zeros(Cw,M);

    for i = 1:M
        w_matrix(:,i) = ridge(y(:,i),out_Mo,kc,0);
    end


    yo = ones(N,1)*w_matrix(1,:) + out_Mo*w_matrix(2:end,:);
    [My,ymo] = max(yo,[],2);
    yth = ymo;
    for i = 1:N
        for j = 1:M
            if ymo(i) == j
                yth(i) = arr_req(j);
            end
        end
    end

    err_ct_train = 0;
    y = labels_train;
    for i=1:N
        if y(i)~=yth(i)
            %disp(i);
            err_ct_train = err_ct_train+1;
        end
    end
    
    %% Data preparation for Test-RC
    
    N = N_test;
    segm = round(ceil(idx/J));
    inp =zeros(N,J*segm);
    inp(1:N,1:idx) = inp_test;

    L = round(idx*(pw+res_t)*D/J);
    L = K*round(ceil(L/K));
    inpm = zeros(N, L);
    for j = 1:N
        for i = 1:segm
            for k = 1:J
                inpm(j,round((i-1)*(pw+res_t)*D+(k-1)*(pw*D/J))+1:round((i-1)*(pw+res_t)*D+k*pw*D/J)) = sqrt(pmax*inp(j,(i-1)*J+k)+poffset);
            end
        end
    end
    
    %% Test-RC
    out_M = zeros(N,L/K,K);
    hvec = 0.8*ones((K-1)*C/2+mod(C,2)/2);
    hmatrix = zeros(N,L/K,(K-1)*C/2+mod(C,2)/2);
    Pmatrix = zeros(N,L/K,(K-1)*C/2+mod(C,2)/2);

    for j = 1:N
        for i = 1:L/K
            inp_ = zeros(K,1);
            for k = 1:K
                inp_(k) = inpm(j,i+(k-1)*L/K);
            end
            [temp, hnvec,P] = RC_KxC(inp_,theta,hvec,dt,K,C);
            out_M(j,i,1:K) = (abs(temp(1:K))).^2;
            hvec = hnvec;
            hmatrix(j,i,:) = hvec(:);
            Pmatrix(j,i,:) = P(:);
            %disp(i);
        end
    end
    
    %% Gather testing data
    out_Mo = reshape(out_M,[N round(L*dt/(detector_int))]);
    temp = out_Mo;
    Ki = 10;
    out_size = round(floor(L*dt/(detector_int)/Ki));
    out_Mo = zeros(N,out_size);
    for j = 1:N
        for i = 1:out_size
            out_Mo(j,i) = mean(temp(j,(i-1)*Ki+1:i*Ki));
        end 
    end 
    out_Mo = cat(2,out_Mo,inp);
    
    %% Test!
    out_MoT = out_Mo;
    y = labels_test;
    N = N_test;
    yo = ones(N,1)*w_matrix(1,:) + out_MoT*w_matrix(2:end,:);
    [My,ymo] = max(yo,[],2);
    yth = ymo;
    for i = 1:N
        for j = 1:M
            if ymo(i) == j
                yth(i) = arr_req(j);
            end
        end
    end
    err_ct_test = 0;
    for i=1:N
        if y(i)~=yth(i)
            %disp(i);
            err_ct_test = err_ct_test+1;
        end
    end

    %% Function returns
    err_train = err_ct_train*100/N_train;
    err_test = err_ct_test*100/N_test;
    Conf = confusionmat(y-1,yth-1);
    toc

end