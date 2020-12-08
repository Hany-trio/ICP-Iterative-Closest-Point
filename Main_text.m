%% Main 2 the second Main program
    % just solve the rotation param
    % the file just as a test

    clc
    clear
    
%% 读取点云数据

    cloud1 = pcread("rabbit_Segment_202011142150.pcd");
%     cloud2 = pcread("rabbit_Segment_202011142150_onlyRO.pcd");   
    cloud2 = pcread("rabbit_Segment_202011142150.pcd");

    cloud1 = cloud1.Location;
    cloud2 = cloud2.Location;
   
    
    R = [0.9814  ,  0.0093  ,  0.1917;
         0.0740  ,  0.9032  , -0.4228;
         -0.1770 ,  0.4292  ,  0.8857 ];
    
     
     cloud2 = ( R * cloud2')';
    
     close all
     pcshow(cloud1 , [0,1,0] , 'MarkerSize' , 20);
     hold on 
     pcshow(cloud2 , [1,0,0] , 'MarkerSize' , 20);
     
     
     
     RO_ICP(cloud1,cloud2)
    
     
     
     
    function [T] = RO_ICP(Refer_cloud,Match_cloud)
    
      % 输入 ： 点云 1 2 文件； 返回： T 变换矩阵
    % Refer_cloud 参考点云 不做变换操作
    % Match_cloud 匹配点云 做变换操作
    
    
    %% 计算初始平移系数
    
    [ Refer_cloud_l , ~ ] = size(Refer_cloud);  %读取点云大小
    [ Match_cloud_l , ~ ] = size(Match_cloud);
    
    CoG_Refer_cloud = sum(Refer_cloud) / Refer_cloud_l;  %计算点云重心
    CoG_Match_cloud = sum(Match_cloud) / Match_cloud_l;    
    
    Initial_Translation_Param = CoG_Refer_cloud - CoG_Match_cloud;  %初始平移量
    Initial_Rotation_Param = diag(ones(1,3));                       %初始旋转量
    
   
  %% 先进行初始变换 
  
    
    
    R = Initial_Rotation_Param;   % 对平移以及旋转进行初始化
    t = [0,0,0]'; 
    T(1:3,1:3) = R;
    T(1:3,4) = t;
    T(4,1:4) = [0,0,0,1];
    
    
   	%% 计算邻近关系 取欧式距离最小值作为对应点 
    
	Refer_kdtreeobj = KDTreeSearcher(Refer_cloud, 'distance' , 'euclidean' );    %建立 Refercloud的kdtree
    
    
    
    
    %% 迭代计算
    
    correspondence = zeros(Match_cloud_l,1);    
    
    
    while(true)   %迭代计算开始
        
    H = zeros(3);  % 定义 Hessian 矩阵
    b = zeros(3,1); % 定义 J * e 
    
        for i = 1 : Match_cloud_l   %确定点对应关系

            [ idx, ~ ] = knnsearch( Refer_kdtreeobj , Match_cloud(i,1:3) , 'dist' , 'euclidean' , 'k' , 1); 
            correspondence(i,1) = idx; 

        end
  
        
        for i = 1 : Match_cloud_l
            
            J =  - ToAntisymmetric_Mat( Match_cloud(i,1:3)' ) ; 
            
            H = H + ( J' * J ); % Hessian
            
            b = b + J' * ( Match_cloud(i,1:3)' - Refer_cloud( correspondence(i,1) , 1:3)' );
            

        end
        
        
        delt = H \ -b;         
        
        fa = delt(1:3,1);
        
        norm_delt = norm(fa);
        a = fa / norm_delt;
        
%         J_ = ( sin( norm_delt ) / norm_delt ) * eye(3) + ( 1 - sin(norm_delt) / norm_delt )*( a * a' )+(( 1 - cos(norm_delt) ) / norm_delt ) * ToAntisymmetric_Mat(a);
        deltR = cos( norm_delt ) * eye(3) + ( 1 - cos( norm_delt ) ) * ( a * a' ) + sin( norm_delt ) * ToAntisymmetric_Mat(a);
 
        Match_cloud = ( deltR * Match_cloud' )';    
        
        close all
        pcshow(Match_cloud , [0,1,0] , 'MarkerSize' , 20);
        hold on 
        pcshow(Refer_cloud , [1,0,0] , 'MarkerSize' , 20);
        
    end
    
    

        
end
    




function [cloud] = Trans(cloud,T)

    [l , w] = size(cloud);
    
    if ( w ~= 3 )
        
        cloud =cloud';
        [l , ~] = size(cloud);
        
    end
    
    
    cloud(:,4) = ones(l,1);
    
    cloud = ( T * cloud' )';
    
    cloud = cloud( :, 1:3 );



end


function [Mat] = ToAntisymmetric_Mat(vector)
    
    [ w , ~ ] = size(vector);
    
    if(w ~= 3)
        
        vector = vector';
        
    end
       
    a1 = vector(1,1);
    a2 = vector(2,1);
    a3 = vector(3,1);
    
    
    Mat = [ 0  , -a3 ,  a2;
            a3 ,  0  , -a1;
           -a2 ,  a1 ,  0];

end