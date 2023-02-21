for paramno  = 1:1
    if paramno == 2
        continue
    end
    noise_index = 1;
    method = 5;
    for trial = 1:50
        heterosig = false;
        heterol = false;
        heteroseps = true;
        driver;
        fprintf('trial: %d\n', trial);
    end
%     for trial = 1:50
%         heterosig = false;
%         heterol = true;
%         driver
%         fprintf('ltrial: %d\n', trial);
%     end
%     save('./result/result_l', 'result_l');
%     for trial = 1:50
%         tic
%         heterosig = true;
%         heterol = true;
%         driver
%         fprintf('lsigtrial: %d\n', trial);
%         toc
%     end
%     save('./result/result_lsig', 'result_lsig');
end
