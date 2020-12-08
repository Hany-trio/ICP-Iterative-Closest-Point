%% 主函数 main
    %20201206 NICP 考虑法线约束的ICP算法
    clc
    clear
    
%% 读取点云数据

    cloud1 = pcread("rabbit_Segment_202011142150.pcd");
    cloud2 = pcread("rabbit_Segment_2020111421501.pcd");   
    
    cloud1 = cloud1.Location;
    cloud2 = cloud2.Location;
    
%     pcshow(cloud1);

%% ICP
    
    n_ICP(cloud1,cloud2)
    
    
    