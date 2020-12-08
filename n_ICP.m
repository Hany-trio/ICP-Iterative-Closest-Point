function [T] = n_ICP(Refer_cloud,Match_cloud)
    
    % ���� �� ���� 1 2 �ļ��� ���أ� T �任����
    % Refer_cloud �ο����� �����任����
    % Match_cloud ƥ����� ���任����
    
    
    %% �����ʼƽ��ϵ��
    
    [ Refer_cloud_l , ~ ] = size(Refer_cloud);  %��ȡ���ƴ�С
    [ Match_cloud_l , ~ ] = size(Match_cloud);
    
    CoG_Refer_cloud = sum(Refer_cloud) / Refer_cloud_l;  %�����������
    CoG_Match_cloud = sum(Match_cloud) / Match_cloud_l;    
    
    Initial_Translation_Param = CoG_Refer_cloud - CoG_Match_cloud;  %��ʼƽ����
    Initial_Rotation_Param = diag(ones(1,3));                       %��ʼ��ת��
    
   
  %% �Ƚ��г�ʼ�任 
  
    
    
    R = Initial_Rotation_Param;   % ��ƽ���Լ���ת���г�ʼ��
    t = [0,0,0]'; 
    T(1:3,1:3) = R;
    T(1:3,4) = t;
    T(4,1:4) = [0,0,0,1];
    
    
   	%% �����ڽ���ϵ ȡŷʽ������Сֵ��Ϊ��Ӧ�� 
    
	Refer_kdtreeobj = KDTreeSearcher(Refer_cloud, 'distance' , 'euclidean' );    %���� Refercloud��kdtree
    
    
    
    
    %% ��������
    
    H = zeros(6);  % ���� Hessian ����
    b = zeros(6,1); % ���� J * e 
    
    correspondence = zeros(Match_cloud_l,1);    
    
    
    while(true)   %�������㿪ʼ
        
    H = zeros(6);  % ���� Hessian ����
    b = zeros(6,1); % ���� J * e 
    
        for i = 1 : Match_cloud_l   %ȷ�����Ӧ��ϵ

            [ idx, ~ ] = knnsearch( Refer_kdtreeobj , Match_cloud(i,1:3) , 'dist' , 'euclidean' , 'k' , 1); 
            correspondence(i,1) = idx; 

        end

%         ps = Match_cloud - CoG_Match_cloud;   % ȥ���Ļ�����
%         pt = Refer_cloud - CoG_Refer_cloud;

        
        
        for i = 1 : Match_cloud_l
            
%             J = [ eye(3) , - ToAntisymmetric_Mat( R * Match_cloud(i,1:3)' + t) ];   % �ſɱȾ���
            J = [ eye(3) , - ToAntisymmetric_Mat( Match_cloud(i,1:3)' ) ]; 
            
            H = H + ( J' * J ); % Hessian
            
            b = b + J' * ( Match_cloud(i,1:3)' - Refer_cloud( correspondence(i,1) , 1:3)' );
            

        end
% 
%         R = T( 1:3 , 1:3);
%         t = T( 1:3 , 4  );
        
        
        delt = H \ -b;         
        
        t_ = delt(1:3,1);
        
        fa = delt(4:6,1);
        
        norm_delt = norm(fa);
        a = fa / norm_delt;
        
        J_ = ( sin( norm_delt ) / norm_delt ) * eye(3) + ( 1 - sin(norm_delt) / norm_delt )*( a * a' )+(( 1 - cos(norm_delt) ) / norm_delt ) * ToAntisymmetric_Mat(a);
        deltR = cos( norm_delt ) * eye(3) + ( 1 - cos( norm_delt ) ) * ( a * a' ) + sin( norm_delt ) * ToAntisymmetric_Mat(a);
        
        
        deltT = [deltR     ,  J_ * t_;
                 zeros(1,3) , 1      ];

        T = deltT ;

        
        Match_cloud = Trans( Match_cloud , T );
        
        
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
