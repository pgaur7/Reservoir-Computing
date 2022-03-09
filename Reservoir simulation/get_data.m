function [images_train, label_train, count_train, images_test, label_test, count_test] =  get_data(arr_req)
    arr_size = size(arr_req,2);
    read_data;
    image_ct = zeros(10,1);
    id_matrix = zeros(60000,10);

    for i = 1:60000
        image_ct(labels(i)+1) = image_ct(labels(i)+1) + 1;
        id_matrix(image_ct(labels(i)+1),labels(i)+1) = i;
    end
    out_ct = 0;
    ct_vec = zeros(arr_size+1,1);
    for i = 1:arr_size
        out_ct = out_ct + image_ct(arr_req(i)+1);
        ct_vec(i+1) = image_ct(arr_req(i)+1);
    end
    ct_cum = cumsum(ct_vec);

    out_labels = zeros(out_ct,1);
    out_images = zeros(28,28,out_ct);
    for j = 1:arr_size
        out_labels(ct_cum(j)+1:ct_cum(j+1)) = arr_req(j);
        out_images(:,:,ct_cum(j)+1:ct_cum(j+1)) = images(:,:,id_matrix(1:image_ct(arr_req(j)+1),arr_req(j)+1));

    end
    
    images_train = out_images;
    label_train = out_labels;
    count_train = ct_vec;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    image_test_ct = zeros(10,1);
    id_test_matrix = zeros(10000,10);

    for i = 1:10000
        image_test_ct(labels_test(i)+1) = image_test_ct(labels_test(i)+1) + 1;
        id_test_matrix(image_test_ct(labels_test(i)+1),labels_test(i)+1) = i;
    end

    out_test_ct = 0;
    ct_test_vec = zeros(arr_size+1,1);
    for i = 1:arr_size
        out_test_ct = out_test_ct + image_test_ct(arr_req(i)+1);
        ct_test_vec(i+1) = image_test_ct(arr_req(i)+1);
    end
    ct_test_cum = cumsum(ct_test_vec);

    out_test_labels = zeros(out_test_ct,1);
    out_test_images = zeros(28,28,out_test_ct);
    for j = 1:arr_size
        out_test_labels(ct_test_cum(j)+1:ct_test_cum(j+1)) = arr_req(j);
        out_test_images(:,:,ct_test_cum(j)+1:ct_test_cum(j+1)) = images_test(:,:,id_test_matrix(1:image_test_ct(arr_req(j)+1),arr_req(j)+1));

    end
    images_test = out_test_images;
    label_test = out_test_labels;
    count_test = ct_test_vec;
end