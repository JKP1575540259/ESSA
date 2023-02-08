function signal = ESSA_I(timeseries,index_miss,L,kk )
%% Extended singular spectrum analysis (Identity)
%  Author: Kunpu Ji 
%  Date: 2021/07/18
%  Input:  timeseries -- timeseries (Ns x 1)
%          index_miss -- gap locations ((N - Ns) x 1)
%          L -- window size
%          kk -- reconstruction order
%  Output: signal-- ectracted signals (N x 1)
%  Note: Ns is the length of available data, N is the length of complete data.
meanv = mean(timeseries);
stdv = sqrt(sum((timeseries-meanv).^2)/(length(timeseries)-1));
timeseries = (timeseries-meanv)/stdv;
% Step1 : Build trayectory matrix
N = length(timeseries) + length(index_miss);
index_total = 1:N;
index_total = index_total';
index_avail = setdiff(index_total,index_miss);
fullseries = zeros(N,1);
fullseries(index_avail) = timeseries;
summ = 0.0;
count = 0;
c = zeros(L,1);
C = zeros(L,L);
for j = 0:L-1
    summ = 0.0;
    count = 0;
    for i = 1:N-j
        if(ismember(i,index_avail)&&ismember(i+j,index_avail))
            summ = summ + fullseries(i)*fullseries(i+j);
            count = count + 1;
        end
    end
    c(j+1) = summ/count;
end
for i = 1:L
    C(i,i:end) = c(1:L-i+1)';
end

for i = 1:L
    for j = 1:i
        C(i,j) = C(j,i);
    end
end

% Step 2: Singular value decomposition
[U,autoval] = eig(C);
[~,i] = sort(-diag(autoval));
U = U(:,i);
V = U;

% Step 3: Constructing coefficient matrix 
CFMM = zeros(N,N);
for k = 1:kk
CFM = zeros(N,N);
for i = 1:L-1
    for j = 1:i-1
        h = 1:j;
        CFM(i,j) = sum(V(h,k).*V(i-j+h,k))/i;
    end
    for j = i:L
        h = 1:i;
        CFM(i,j) = sum(V(h,k).*V(j-i+h,k))/i;
    end
    for j = L+1:i-1+L
        h = 1:i+L-j;
        CFM(i,j) = sum(V(h,k).*V(j-i+h,k))/i;
    end
end

for i = L:N-L+1
    for j=i-L+1:i-1
        h = 1:L-i+j;
        CFM(i,j) = sum( V(h,k).*V(i-j+h,k))/L;
    end
    CFM(i,i) = sum(V(1:L,k).^2)/L;
    for j = i+1:i-1+L
        h = 1:i+L-j;
        CFM(i,j) = sum(V(h,k).*V(j-i+h,k))/L;
    end
end

for i = N-L+2:N
    for j = i-L+1:N-L
        h = 1:L-i+j;
        CFM(i,j) = sum(V(h,k).*V(i-j+h,k))/(N-i+1);
    end
    for j = N-L+1:i
        h = i:N;
        CFM(i,j) = sum(V(L+j-h,k).*V(i+L-h,k))/(N-i+1);
    end
    for j = i+1:N
        h = j:N;
        CFM(i,j) = sum(V(L+j-h,k).*V(i+L-h,k))/(N-i+1);
    end
end
CFMM = CFMM + CFM;
end
AA = eye(N,N) - CFMM;
A1 = AA(:,index_avail);
A2 = AA(:,index_miss);

% Step 4: Solving the equation and extracting the signals
value = -(A2'*A2)\(A2'*A1*timeseries);
signal = CFMM(:,index_avail)*timeseries + CFMM(:,index_miss)*value;
signal = signal*stdv+meanv;
end   





