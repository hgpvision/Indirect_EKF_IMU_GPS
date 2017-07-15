%%本程序为GPS+INS在MAV导航上的融合程序惯导解算（滤波融合）部分，采用的是间接卡尔曼滤波
%%其中惯导用来进行状态预测，GPS用来滤波矫正（即量测方程中只有GPS量测量，惯导测量值直接用在状态预测方程中）
%%使用状态预测方程及惯导测量值进行状态预测，求出偏差状态预测方程的有关系数
%%每次使用偏差状态预测方程预测得到偏差，如果有GPS测量值，则进行量测更新，并对状态进行偏差矫正
%%没有GPS测量值，则没有量测更新，也没有状态偏差矫正，想当于没有进行滤波  

%%参考文献：
%%[1]Monocular Vision for Long-term Micro Aerial State Estimation:A
%%Compendium, Stephan Weiss,2013
%%[2]：卡尔曼滤波与组合导航原理，秦永元，P49，用到了离散系统噪声方差计算公式
%%[3]：Indirect Kalman Filter for 3D Attitude Estimation, Trawny, 2005 
%%2015/12/30
%%*****************************************************************************************************%%

classdef InsSolver
    properties
        attiCalculator;
        Qc;
        Rc;
    end
    methods
        function o = InsSolver(Qc0,Rc0)
            o.attiCalculator = AttitudeBase();
            o.Qc=Qc0;
            o.Rc=Rc0;
        end
        
        
        function dy = imuDynamics( o,t,state,x )
            %四元数的imu动力学方程
            acc = x(1:3);
            gyro = x(4:6);
            v = state(4:6);
            quat = state(7:10);
            ba = state(11:13);
            bw = state(14:16);
            
%             na = zeros(3,1);
%             nw = zeros(3,1);
%             na = randn(3,1)*0.1;
%             nw = randn(3,1)*0.1;
            
            g = [0 0 -9.8]';
            cnb = o.attiCalculator.quat2cnb(quat);
            omega = o.attiCalculator.QuatMulitMati(quat,[0;gyro - bw]);
            dy = [v;cnb'*(acc - ba)+g;omega/2;zeros(6,1)];
            
        end
        
        function [state,predErrorState] = imu2state(o,acc_nose,gyro_noise,gps_noise,state0,errorstate0,tspan,step,isODE45)
            %acc和gyro分别为加表陀螺仪数据
            %state0：        初始状态
            %errorstate0：   偏差初始状态
            %Qc0：           误差方差阵
            %tspan：         时间
            %step：          步长
            %isODE45：       使用ODE45或者单步积分
            if nargin < 9
                isODE45 = 1;
            end
            
            N = length(tspan);
            state = zeros(N,length(state0));
            predErrorState=zeros(N,length(errorstate0));
        
            predP=cell(N,1);    %偏差预测误差方差
            
            state(1,:) = state0';
            predErrorState(1,:)= errorstate0';
            predP{1}=zeros(15); %初偏差预测误差方差
            
            for i = 1:N - 1
                state(i+1,:) = o.statePrediction(state(i,:)',...
                                                [tspan(i) tspan(i+1)],...
                                                [acc_nose(i,:) gyro_noise(i,:)]',...
                                                isODE45)';
                                            
                predQ=o.QuatNormalize(state(i+1,7:10));         %四元数归一化（进入本次迭代归一化一次）

                predCbn=o.attiCalculator.quat2cnb(predQ);

                [Fd,Gc,Fc]=o.getExponentMatFd(predCbn,step,acc_nose(i+1,:)-state(i+1,14:16),gyro_noise(i+1,:)-state(i+1,11:13));%acc取i行还是i+1行？

%                 predErrorState(i+1,:)=predErrorState(i,:)*Fd';  %误差状态预测
                predErrorState(i+1,:)=zeros(1,15);
                
                Qd=o.getPredCovarianceMatQd(Gc,Fc,o.Qc,step);
                
                predP{i+1}=Fd*predP{i}*Fd'+Qd;                  %预测偏差方差矩阵计算
                
%                 state(i+1,7:10)=o.QuatNormalize(state(i+1,7:10));%四元数归一化（进入下次迭代之前归一化一次，动力学方程中没有归一化）
                
                if gps_noise(i+1,1)==1
                    %如若想看看不进行滤波的结果，可以设置始终不进入该if条件语句，运行即可
                    %如果有GPS测量值，则进行两侧更新以及状态偏差矫正
                    %其中GPS测量值对姿态是不可观的？不对四元数更新（否则姿态反而发散，还对位置和速度有一定影响）
                                        
                    [postErrorState,postP]=o.MeasurementUpdate((gps_noise(i+1,2:4)-state(i+1,1:3)),predP{i+1},predErrorState(i+1,:));
                    
                    predErrorState(i+1,:)=postErrorState;   %用新得到的后验偏差误差方差阵更新预测过程中得到的先验方差阵，以供后面迭代使用
                    predP{i+1}=postP;
                    
                    state(i+1,1:6)=state(i+1,1:6)+predErrorState(i+1,1:6);
                    
%                     %----------
%                     %四元数归一化方式1：再归一化偏差四元数，再进行四元数相乘矫正
%                     delta_q=0.5*predErrorState(i+1,7:9);    %由角度误差更新得到的四元数矢量部分误差更新
%                     temp=delta_q*delta_q';
%                     if temp>1
%                         delta_q_hat=[1,delta_q]/sqrt(1+temp);
%                     else
%                         delta_q_hat=[sqrt(1-temp),delta_q];
%                     end
%                     state(i+1,7:10)=o.attiCalculator.QuatMulitMat(state(i+1,7:10),delta_q_hat');
%                     %----------
%                     
                    %----------
                    %四元数归一化方式2：先四元数相乘，再归一化
                    state(i+1,7:10)=o.attiCalculator.QuatMulitMat(state(i+1,7:10),[1,0.5*predErrorState(i+1,7:9)]');
                    state(i+1,7:10)=o.QuatNormalize(state(i+1,7:10));
                    %----------
                    
                    state(i+1,11:16)=state(i+1,11:16)+predErrorState(i+1,10:15);    %这个是否需要进行偏差矫正以后视情况而定，涉及姿态的就不用了
%                 else
%                     %没有GPS量测值时，不用进行偏差矫正，不然结果发散
%                     state(i+1,1:6)=state(i+1,1:6)+predErrorState(i+1,1:6);
%                     state(i+1,11:16)=state(i+1,11:16)+predErrorState(i+1,10:15);
                end
            end
        end
        
        function state = statePrediction(o,state0,tspan,imudata,isODE45)
            %一步状态预测
            
            if isODE45
                %ode45速度很慢
                %这里的每一步就是状态预测
                [t,y] = ode45(@o.imuDynamics,tspan,state0,[],imudata);
                state = y(end,:);
            else
                step = tspan(2) - tspan(1);
                dy = o.imuDynamics([],state0,imudata);
                state = state0 + step*dy;
            end
        end
        
        function [Fd,Gc,Fc]=getExponentMatFd(o,predCbn,step,a_hat,omega_hat)
            %计算误差状态方程（偏差量测方程）的离散化雅克比矩阵Fd，以及离散化噪声雅克比矩阵Gc（包含线性化和离散化两个过程）
            %preCbn：    预测状态方程（预测量测方程）得到的体系相对导航系的旋转矩阵
            %step：      步长
            %a_hat：     a_hat=acc-ba_hat，加表线性加速度测量值减去预测漂移值，这里ba_hat暂定恒为0
            %omega_hat： omega_hat=gyro-bw_hat，陀螺仪角速度测量值减去预测漂移值，这里bw_hat暂定恒为0
            
            skew_a_hat=o.attiCalculator.a2skew_a(a_hat);
            skew_omega_hat=o.attiCalculator.a2skew_a(omega_hat);
            
            A=-predCbn'*skew_a_hat*(step^2/2*eye(3)-step^3/6*skew_omega_hat+step^4/24*skew_omega_hat^2);
            B=-predCbn'*skew_a_hat*(-step^3/6*eye(3)+step^4/24*skew_omega_hat-step^5/120*skew_omega_hat^2);
            C=-predCbn'*skew_a_hat*(step*eye(3)-step^2/2*skew_omega_hat+step^3/6*skew_omega_hat^2);
            D=-A;
            E=eye(3)-step*skew_omega_hat+step^2/2*skew_omega_hat^2;
            F=-step*eye(3)+step^2/2*skew_omega_hat-step^3/6*skew_omega_hat^2;
            
            Fd=[eye(3) step*eye(3) A B -predCbn'*step*step/2
                zeros(3) eye(3) C D -predCbn'*step
                zeros(3) zeros(3) E F zeros(3)
                zeros(3) zeros(3) zeros(3) eye(3) zeros(3)
                zeros(3) zeros(3) zeros(3) zeros(3) eye(3)];
            Gc=[zeros(3) zeros(3) zeros(3) zeros(3)
                -predCbn' zeros(3) zeros(3) zeros(3)
                zeros(3) zeros(3) -eye(3) zeros(3)
                zeros(3) zeros(3) zeros(3) eye(3)
                zeros(3) eye(3) zeros(3) zeros(3)];
            Fc=[zeros(3) eye(3) zeros(3) zeros(3) zeros(3)
                zeros(3) zeros(3) -predCbn'*skew_a_hat zeros(3) -predCbn'
                zeros(3) zeros(3) -skew_omega_hat -eye(3) zeros(3)
                zeros(3) zeros(3) zeros(3) zeros(3) zeros(3)
                zeros(3) zeros(3) zeros(3) zeros(3) zeros(3)];
        end
        
        function Qd=getPredCovarianceMatQd(o,Gc,Fc,Qc,step)
            %计算离散系统预测噪声方差阵Qd
            %Fc：    偏差预测状态方程雅克比矩阵
            %Gc：    偏差预测状态方程过程噪声项雅克比矩阵
            %Qc：    过程噪声方差
            %step：  步长

            Q_hat=Gc*Qc*Gc';
            Qd=step*Q_hat+step^2*(0.5*Fc*Q_hat+Q_hat*0.5*Fc');
        end
        
        function [postErrorState,postP]=MeasurementUpdate(o,error_gps_noise,predP,predErrorState)
            %量测更新函数
            %error_gps_noise：   GPS测量值偏差（即测量值与估计值之差，此处使用间接卡尔曼滤波）
            %predP：             先验偏差误差方差
            %predErrorState：    先验偏差估计值
            
            H=[eye(3) zeros(3,12)];                 %量测矩阵
            K=predP*H'*(H*predP*H'+o.Rc)^-1;        %kalman增益
            
            postErrorState=predErrorState+(error_gps_noise-predErrorState*H')*K';%后验偏差状态估计
            
            postP=predP-K*H*predP;                  %后验偏差误差方差阵
        end
        
        function Q=QuatNormalize(o,q)
            %四元数归一化
            %q：输入四元数
            
            Q=q/norm(q);
        end
    end
    
end