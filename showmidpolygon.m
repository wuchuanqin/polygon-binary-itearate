%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% showmidpolygon(3, 2, 200, 1, 300, 0);
% showmidpolygon(4, 2, 200, 1, 300, 0);
% showmidpolygon(5, 2, 200, 1, 300, 0);
% showmidpolygon(3, 2, 200, 1, 300, 2);
% showmidpolygon(4, 2, 200, 1, 300, 2);
% showmidpolygon(5, 2, 200, 1, 300, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%N: 正多行顶点数
%I: 迭代次数
%L: 正多边形原始边长
%K: 比例系数
%A: 坐标轴范围
%m: 绘图模式
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showmidpolygon(N, I, L, K, A, m)
    if (m<0) || (m>I)
        return;
    end;
    
    [P, IofP]=midpolygon(N, I, L, K);
             
    hold on;
    plot([P(IofP(1,1):IofP(1,2)), P(IofP(1,1))], '-*b');
    axis([0-A, A, 0-A, A]);
    axis equal;
        
    if m == 0
        for ki=1:I+1              
            plot([P(IofP(ki,1):IofP(ki,1)+N-1), P(IofP(ki,1))], '-*r');
            ki2=IofP(ki, 1)+N;
            while ki2<IofP(ki, 2)
                plot([P(ki2:ki2+3-1), P(ki2)], '-*r');
                ki2=ki2+3;
            end;
        end;
    else        
        ki=m+1;
        plot([P(IofP(ki, 1):IofP(ki, 1)+N-1), P(IofP(ki, 1))], '-*r');
        ki2=IofP(ki, 1)+N+N*3;
        while ki2<IofP(ki, 2)
            plot([P(ki2:ki2+3-1), P(ki2)], '-*r');
            ki2=ki2+4*3;
        end;        
    end;
    hold off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, IofP] = midpolygon(N, I, L, K)
    if (N<1) || (N>8) || (I<1) || (I>5) || (L<1) || (L>1000)
       P=0+0*1j;
       return;
    end;

    p0=L*exp(1j*2*pi*(0:N-1)/N);    

    c0=1;
    ck=0;
    sum=0          
    IofP=zeros(I+1, 2);    
    IofP(1,:)=[1, N];
    for kt=1:I
        ck=(N+1)+(c0-1)*4;
        IofP(kt+1,:)=[IofP(kt,2)+1, IofP(kt,2)+(N+(ck-1)*3)];
        c0=ck;
        sum=sum+(N+(ck-1)*3);
    end;        
    pm=zeros(1, N+sum);   
    pm(1:N)=p0;
    
    c0=1;
    ck=0;    
    knl=1;
    kns=N+1;
    cp=zeros(1, I);    
    for ki=1:I        
        %%%%%%%%迭代正>3多边形迭代一次%%%%%%%%        
        %构造当前次需要迭代的正>3多角形顶点
        pn_t=pm(knl:knl+N-1);

        %进行一次迭代
        pn=midvector(pn_t,N,K);          
        
        %存储当前次迭代所得正>3多边形顶点        
        pm(kns:kns+N-1)=pn;  

        %存储当前次正>3多边形迭代所得正=3多边形顶点        
        kns_t=kns+N;
        for kn=1:N            
            pm(kns_t:kns_t+2)=[pn(kn), pn_t(mod(kn,N)+1), pn(mod(kn,N)+1)];                        
            kns_t=kns_t+3;
        end;

        %如果需要进行正=3多边形迭代，那么ki3_s为存储其迭代所得顶点的在pm中的开始索引值
        k3s=kns_t;        
        %如果需要进行正=3多边形迭代，那么ki3_l为加载所需迭代顶点在pm中的开始索引值
        k3l=knl+N;

        if ki > 1           
            %%%%%%%%迭代正=3多边形%%%%%%%%
            k3l_t=k3l;
            k3s_t=k3s;            
            while(k3l_t<kns)
                %构造当前次需要迭代的正=3三角形顶点
                p3_t=pm(k3l_t:k3l_t+2);

                %进行一次迭代
                p3=midvector(p3_t, 3, K);

                %存储当前次正=3多边形迭代所得第一个正=3多边形顶点
                pm(k3s_t:k3s_t+2)=p3;

                %存储当前次正=3多边形迭代所得其余正=3多边形顶点
                k3s_t=k3s_t+3;
                for k3=1:3
                    pm(k3s_t:k3s_t+2)=[p3(k3), p3_t(mod(k3,3)+1), p3(mod(k3,3)+1)];
                    k3s_t=k3s_t+3;
                end;

                k3l_t=k3l_t+3;
            end;
        end;

        ck=(N+1)+(c0-1)*4;
        c0=ck;

        knl=kns;        
        kns=kns+(N+(ck-1)*3);
    end

    P=pm(1:kns-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out_p]=midvector(in_p, N, k)
    out_p=zeros(1,N);
    for i=1:N
        out_p(i)=in_p(i)+in_p(mod(i,N)+1);
    end;        
    out_p=out_p*k/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
