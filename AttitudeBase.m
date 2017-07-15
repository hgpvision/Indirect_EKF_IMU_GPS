classdef AttitudeBase
    %该类封装了基本的姿态变化关系函数及相关矩阵计算等
    
    methods
                
        function cnb = a2cnb(o,atti)
            %输入：姿态角(pitch, roll, yaw)
            %输出：导航系到体系的姿态变换矩阵(navigation frame to body frame)
            sp = sin(atti(1));cp = cos(atti(1));    %pitch
            sr = sin(atti(2));cr = cos(atti(2));    %roll
            sy = sin(atti(3));cy = cos(atti(3));    %yaw
            cnb = zeros(3);
            cnb(1,1) = cr*cy - sr*sp*sy;
            cnb(1,2) = cr*sy + sr*sp*cy;
            cnb(1,3) = -sr*cp;
            cnb(2,1) = -cp*sy;
            cnb(2,2) = cp*cy;
            cnb(2,3) = sp;
            cnb(3,1) = sr*cy + cr*sp*sy;
            cnb(3,2) = sr*sy - cr*sp*cy;
            cnb(3,3) = cr*cp;
        end
        
        function a = cnb2atti(o,cnb )
            %a2cnb()的反函数
            gama = -atan(cnb(1,3)/cnb(3,3));
            theta = asin(cnb(2,3));
            ctheta = cos(theta);
            psi = acos(cnb(2,2)/ctheta);
            if -cnb(2,1)/ctheta < 0
                psi = -psi;
            end
            a = [theta gama psi];      
        end
        
        function cnb = quat2cnb( o,quat )
            %输入：姿态四元数
            %输出：导航系到体系的姿态变换矩阵(navigation frame to body frame)
            cnb(1,1) = quat(1)*quat(1) + quat(2)*quat(2) - quat(3)*quat(3) - quat(4)*quat(4);
            cnb(1,2) = 2*(quat(2)*quat(3)+quat(1)*quat(4));
            cnb(1,3) = 2*(quat(2)*quat(4)-quat(1)*quat(3));
            cnb(2,1) = 2*(quat(2)*quat(3)-quat(1)*quat(4));
            cnb(2,2) = quat(1)*quat(1)-quat(2)*quat(2)+quat(3)*quat(3)-quat(4)*quat(4);
            cnb(2,3) = 2*(quat(3)*quat(4)+quat(1)*quat(2));
            cnb(3,1) = 2*(quat(2)*quat(4)+quat(1)*quat(3));
            cnb(3,2) = 2*(quat(3)*quat(4)-quat(1)*quat(2));
            cnb(3,3) = quat(1)*quat(1) - quat(2)*quat(2) - quat(3)*quat(3) + quat(4)*quat(4);
        end
        
        function A = w2Datti( o,atti )
            %输入：姿态角
            %输出：姿态角速率与角速度的关系矩阵
            % A = [-sin(atti(2)), 0, cos(atti(2));
            %     cos(atti(2))*cos(atti(1)),0,sin(atti(2))*cos(atti(1));
            %     sin(atti(1))*sin(atti(2)), cos(atti(1)), -sin(atti(1))*cos(atti(2))]/cos(atti(1));
            sp = sin(atti(1));cp = cos(atti(1));
            sr = sin(atti(2));cr = cos(atti(2));
            A = [cp*cr, 0, sr*cp;
                sp*sr, cp, -cr*sp;
                -sr, 0, cr]/cp;
        end
        
        function T = Datti2w(o,atti)
            %输入：姿态角
            %输出：角速度与姿态角速率的关系矩阵
            sp = sin(atti(1));cp = cos(atti(1));
            sr = sin(atti(2));cr = cos(atti(2));
            T = [cr, 0 -sr*cp;
                0, 1,sp;
                sr, 0, cp*cr];       
        end
        
        function deltheta = w2omega( o,w )
            %输入：角速度的变化量
            %输出：四元数变化率与四元数的关系矩阵
            delthetax = w(1);
            delthetay = w(2);
            delthetaz = w(3);
            deltheta(1,1) = 0;
            deltheta(1,2) = -delthetax;
            deltheta(1,3) = -delthetay;
            deltheta(1,4) = -delthetaz;
            deltheta(2,1) = delthetax;
            deltheta(2,2) = 0;
            deltheta(2,3) = delthetaz;
            deltheta(2,4) = -delthetay;
            deltheta(3,1) = delthetay;
            deltheta(3,2) = -delthetaz;
            deltheta(3,3) = 0;
            deltheta(3,4) = delthetax;
            deltheta(4,1) = delthetaz;
            deltheta(4,2) = delthetay;
            deltheta(4,3) = -delthetax;
            deltheta(4,4) = 0;
            
        end
        
        function skew_a=a2skew_a(o,a)
            %a矢量的反对称阵计算
            skew_a=[0 -a(3) a(2)
                    a(3) 0 -a(1)
                    -a(2) a(1) 0];
        end
        
        function pq=QuatMulitMat(o,p,q)
            %p、q四元数乘法计算
            matp=[p(1) -p(2) -p(3) -p(4)
                p(2) p(1) -p(4) p(3)
                p(3) p(4) p(1) -p(2)
                p(4) -p(3) p(2) p(1)];
            pq=matp*q;
        end
        
        function pq=QuatMulitMati(o,p,q)
            %p、q四元数乘法计算
            matiq=[q(1) -q(2) -q(3) -q(4)
                q(2) q(1) q(4) -q(3)
                q(3) -q(4) q(1) q(2)
                q(4) q(3) -q(2) q(1)];
            pq=matiq*p;
        end
    end  
end
