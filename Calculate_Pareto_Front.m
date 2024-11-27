function pareto_front = Calculate_Pareto_Front(obj_fun, num_n, num_obj)

switch obj_fun
    case 'ZDT1'
        f(:,1)    = (0:1/(num_n-1):1)';
        f(:,2)    = 1 - f(:,1).^0.5;
        pareto_front = f;
    case 'ZDT2'
        f(:,1)    = (0:1/(num_n-1):1)';
        f(:,2)    = 1-f(:,1).^2;
        pareto_front = f;
    case 'ZDT3'
        f(:,1)    = (0:1/(num_n-1):1)';
        f(:,2)    = 1-f(:,1).^0.5-f(:,1).*sin(10*pi*f(:,1));
        f = f(Front_Rank(f,1)==1,:);
        pareto_front = f;
    case 'ZDT4'
        f(:,1)    = (0:1/(num_n-1):1)';
        f(:,2)    = 1-f(:,1).^0.5;
        pareto_front = f;
    case 'ZDT6'
        minf1     = 0.280775;
        f(:,1)    = (minf1:(1-minf1)/(num_n-1):1)';
        f(:,2)    = 1-f(:,1).^2;
        pareto_front = f;
    case 'IDTLZ1'
        f = UniformPoint(num_n,num_obj)/2;
        pareto_front = (1-f)/2;
    case 'IDTLZ2'
        f = UniformPoint(num_n,num_obj);
        pareto_front = 1 - f./repmat(sqrt(sum(f.^2,2)),1,num_obj);
    case 'DTLZ1'
        f = UniformPoint(num_n,num_obj)/2;
        pareto_front = f;
    case 'DTLZ2'
        f = UniformPoint(num_n,num_obj);
        pareto_front = f./repmat(sqrt(sum(f.^2,2)),1,num_obj);
    case 'DTLZ3'
        f = UniformPoint(num_n,num_obj);
        pareto_front = f./repmat(sqrt(sum(f.^2,2)),1,num_obj);
    case 'DTLZ4'
        f = UniformPoint(num_n,num_obj);
        pareto_front = f./repmat(sqrt(sum(f.^2,2)),1,num_obj);
    case 'DTLZ5'
        f = [0:1/(num_n-1):1;1:-1/(num_n-1):0]';
        f = f./repmat(sqrt(sum(f.^2,2)),1,size(f,2));
        f = [f(:,ones(1,num_obj-2)),f];
        pareto_front = f./sqrt(2).^repmat([num_obj-2,num_obj-2:-1:0],size(f,1),1);
    case 'DTLZ6'
        f = [0:1/(num_n-1):1;1:-1/(num_n-1):0]';
        f = f./repmat(sqrt(sum(f.^2,2)),1,size(f,2));
        f = [f(:,ones(1,num_obj-2)),f];
        pareto_front = f./sqrt(2).^repmat([num_obj-2,num_obj-2:-1:0],size(f,1),1);
    case 'DTLZ7'
        interval     = [0,0.251412,0.631627,0.859401];
        median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
        %         X            = ReplicatePoint(10000,num_obj-1);
        if num_obj  > 2
            sample_num = (ceil(num_n^(1/(num_obj - 1))))^(num_obj - 1);
            Gap       = 0:1/(sample_num^(1/(num_obj - 1))-1):1;
            eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:(num_obj - 1))))
            eval(sprintf('X=[%s];',sprintf('c%d(:),',1:(num_obj - 1))))
        else
            X = (0:1/(num_n-1):1)';
        end
        X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
        X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
        pareto_front            = [X,2*(num_obj-sum(X/2.*(1+sin(3*pi.*X)),2))];
    case 'WFG1'
        R = UniformPoint(num_n,num_obj);
        c = ones(size(R,1),num_obj);
        for i = 1 : size(R,1)
            for j = 2 : num_obj
                temp = R(i,j)/R(i,1)*prod(1-c(i,num_obj-j+2:num_obj-1));
                c(i,num_obj-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        temp = (1-sin(pi/2*x(:,2))).*R(:,num_obj)./R(:,num_obj-1);
        a = 0 : 0.0001 : 1;
        E = abs(temp*(1-cos(pi/2*a))-1+repmat(a+cos(10*pi*a+pi/2)/10/pi,size(x,1),1));
        [~,rank] = sort(E,2);
        for i = 1 : size(x,1)
            x(i,1) = a(min(rank(i,1:10)));
        end
        R      = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
        R(:,num_obj) = 1-x(:,1)-cos(10*pi*x(:,1)+pi/2)/10/pi;
        pareto_front    = repmat(2:2:2*num_obj,size(R,1),1).*R;
    case 'WFG2'
        R = UniformPoint(num_n,num_obj);
        c = ones(size(R,1),num_obj);
        for i = 1 : size(R,1)
            for j = 2 : num_obj
                temp = R(i,j)/R(i,1)*prod(1-c(i,num_obj-j+2:num_obj-1));
                c(i,num_obj-j+1) = (temp^2-temp+sqrt(2*temp))/(temp^2+1);
            end
        end
        x = acos(c)*2/pi;
        temp = (1-sin(pi/2*x(:,2))).*R(:,num_obj)./R(:,num_obj-1);
        a = 0 : 0.0001 : 1;
        E = abs(temp*(1-cos(pi/2*a))-1+repmat(a.*cos(5*pi*a).^2,size(x,1),1));
        [~,rank] = sort(E,2);
        for i = 1 : size(x,1)
            x(i,1) = a(min(rank(i,1:10)));
        end
        R      = convex(x);
        R(:,num_obj) = disc(x);
        R      = R(NDSort(R,1)==1,:);
        pareto_front    = repmat(2:2:2*num_obj,size(R,1),1).*R;
    case 'WFG3'
        X = linspace(0,1,num_n)';
        X = [X,zeros(num_n,num_obj-2)+0.5,zeros(num_n,1)];
        R = linear(X);
        pareto_front = repmat(2:2:2*num_obj,size(R,1),1).*R;
    case 'WFG4'
        R = UniformPoint(num_n,num_obj);
        R = R./repmat(sqrt(sum(R.^2,2)),1,num_obj);
        pareto_front = repmat(2:2:2*num_obj,size(R,1),1).*R;
    case 'WFG5'
        R = UniformPoint(num_n,num_obj);
        R = R./repmat(sqrt(sum(R.^2,2)),1,num_obj);
        pareto_front = repmat(2:2:2*num_obj,size(R,1),1).*R;
    case 'WFG6'
        R = UniformPoint(num_n,num_obj);
        R = R./repmat(sqrt(sum(R.^2,2)),1,num_obj);
        pareto_front = repmat(2:2:2*num_obj,size(R,1),1).*R;
    case 'WFG7'
        R = UniformPoint(num_n,num_obj);
        R = R./repmat(sqrt(sum(R.^2,2)),1,num_obj);
        pareto_front = repmat(2:2:2*num_obj,size(R,1),1).*R;
    case 'WFG8'
        R = UniformPoint(num_n,num_obj);
        R = R./repmat(sqrt(sum(R.^2,2)),1,num_obj);
        pareto_front = repmat(2:2:2*num_obj,size(R,1),1).*R;
    case 'WFG9'
        R = UniformPoint(num_n,num_obj);
        R = R./repmat(sqrt(sum(R.^2,2)),1,num_obj);
        pareto_front = repmat(2:2:2*num_obj,size(R,1),1).*R;
        
end
end








function Output = linear(x)
Output = fliplr(cumprod([ones(size(x,1),1),x(:,1:end-1)],2)).*[ones(size(x,1),1),1-x(:,end-1:-1:1)];
end

function Output = convex(x)
Output = fliplr(cumprod([ones(size(x,1),1),1-cos(x(:,1:end-1)*pi/2)],2)).*[ones(size(x,1),1),1-sin(x(:,end-1:-1:1)*pi/2)];
end

function Output = disc(x)
Output = 1-x(:,1).*(cos(5*pi*x(:,1))).^2;
end

function [FrontNo,MaxFNo] = NDSort(varargin)
PopObj = varargin{1};
[N,M]  = size(PopObj);
if nargin == 2
    nSort  = varargin{2};
else
    PopCon = varargin{2};
    nSort  = varargin{3};
    Infeasible           = any(PopCon>0,2);
    PopObj(Infeasible,:) = repmat(max(PopObj,[],1),sum(Infeasible),1) + repmat(sum(max(0,PopCon(Infeasible,:)),2),1,M);
end
if M < 3 || N < 500
    % Use efficient non-dominated sort with sequential search (ENS-SS)
    [FrontNo,MaxFNo] = ENS_SS(PopObj,nSort);
else
    % Use tree-based efficient non-dominated sort (T-ENS)
    [FrontNo,MaxFNo] = T_ENS(PopObj,nSort);
end
end

function [FrontNo,MaxFNo] = ENS_SS(PopObj,nSort)
[PopObj,~,Loc] = unique(PopObj,'rows');
Table   = hist(Loc,1:max(Loc));
[N,M]   = size(PopObj);
FrontNo = inf(1,N);
MaxFNo  = 0;
while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
    MaxFNo = MaxFNo + 1;
    for i = 1 : N
        if FrontNo(i) == inf
            Dominated = false;
            for j = i-1 : -1 : 1
                if FrontNo(j) == MaxFNo
                    m = 2;
                    while m <= M && PopObj(i,m) >= PopObj(j,m)
                        m = m + 1;
                    end
                    Dominated = m > M;
                    if Dominated || M == 2
                        break;
                    end
                end
            end
            if ~Dominated
                FrontNo(i) = MaxFNo;
            end
        end
    end
end
FrontNo = FrontNo(:,Loc);
end

function [FrontNo,MaxFNo] = T_ENS(PopObj,nSort)
[PopObj,~,Loc] = unique(PopObj,'rows');
Table     = hist(Loc,1:max(Loc));
[N,M]     = size(PopObj);
FrontNo   = inf(1,N);
MaxFNo    = 0;
Forest    = zeros(1,N);
Children  = zeros(N,M-1);
LeftChild = zeros(1,N) + M;
Father    = zeros(1,N);
Brother   = zeros(1,N) + M;
[~,ORank] = sort(PopObj(:,2:M),2,'descend');
ORank     = ORank + 1;
while sum(Table(FrontNo<inf)) < min(nSort,length(Loc))
    MaxFNo = MaxFNo + 1;
    root   = find(FrontNo==inf,1);
    Forest(MaxFNo) = root;
    FrontNo(root)  = MaxFNo;
    for p = 1 : N
        if FrontNo(p) == inf
            Pruning = zeros(1,N);
            q = Forest(MaxFNo);
            while true
                m = 1;
                while m < M && PopObj(p,ORank(q,m)) >= PopObj(q,ORank(q,m))
                    m = m + 1;
                end
                if m == M
                    break;
                else
                    Pruning(q) = m;
                    if LeftChild(q) <= Pruning(q)
                        q = Children(q,LeftChild(q));
                    else
                        while Father(q) && Brother(q) > Pruning(Father(q))
                            q = Father(q);
                        end
                        if Father(q)
                            q = Children(Father(q),Brother(q));
                        else
                            break;
                        end
                    end
                end
            end
            if m < M
                FrontNo(p) = MaxFNo;
                q = Forest(MaxFNo);
                while Children(q,Pruning(q))
                    q = Children(q,Pruning(q));
                end
                Children(q,Pruning(q)) = p;
                Father(p) = q;
                if LeftChild(q) > Pruning(q)
                    Brother(p)   = LeftChild(q);
                    LeftChild(q) = Pruning(q);
                else
                    bro = Children(q,LeftChild(q));
                    while Brother(bro) < Pruning(q)
                        bro = Children(q,Brother(bro));
                    end
                    Brother(p)   = Brother(bro);
                    Brother(bro) = Pruning(q);
                end
            end
        end
    end
end
FrontNo = FrontNo(:,Loc);
end



function [W,N] = UniformPoint(N,M)
%UniformPoint - Generate a set of uniformly distributed points on the unit
%hyperplane
%
%   [W,N] = UniformPoint(N,M) returns approximate N uniformly distributed
%   points with M objectives.
%
%   Example:
%       [W,N] = UniformPoint(275,10)

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

H1 = 1;
while nchoosek(H1+M,M-1) <= N
    H1 = H1 + 1;
end
W = nchoosek(1:H1+M-1,M-1) - repmat(0:M-2,nchoosek(H1+M-1,M-1),1) - 1;
W = ([W,zeros(size(W,1),1)+H1]-[zeros(size(W,1),1),W])/H1;
if H1 < M
    H2 = 0;
    while nchoosek(H1+M-1,M-1)+nchoosek(H2+M,M-1) <= N
        H2 = H2 + 1;
    end
    if H2 > 0
        W2 = nchoosek(1:H2+M-1,M-1) - repmat(0:M-2,nchoosek(H2+M-1,M-1),1) - 1;
        W2 = ([W2,zeros(size(W2,1),1)+H2]-[zeros(size(W2,1),1),W2])/H2;
        W  = [W;W2/2+1/(2*M)];
    end
end
W = max(W,1e-6);
N = size(W,1);
end