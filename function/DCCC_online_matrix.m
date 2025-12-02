function [ Fc ] = DCCC_online_matrix( data, s, M )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if size(data,1) > size(data,2)
    data = data';
end

if size(s,1) > size(s,2)
    s = s';
end

% M = max(s)*4;
nch = size(data,1);
N = size(data,2);
ns = length(s);
smin = min(s);
smax = max(s);
Fc = struct('t',[],'DCCA',[],'DCCC',[]);

% sumx_tmp = zeros(ns, 4);
% initializing temporary variables
n_smax = floor(M/smax);
DCCA_tmp = zeros(nch,nch,ns,n_smax);
denom_tmp = zeros(nch,nch,ns);
for si = 1:ns
    denom_tmp(:,:,si) = repmat(1/(n_smax*smax/s(si)), [nch,nch]);
end

sx = zeros(nch,ns);
sxi = zeros(nch,ns);
sx2 = zeros(nch,nch,ns);

i = zeros(nch,ns);
wi = 0;
Mi = 0;
Fi = 0;

for t = 1:N
    if t == 1
        Xi = data(:,t);
        wi = 1;
        Mi = 1;
    else
        Xi = Xi + data(:,t); % cumulative summation
%         Xi = data(:,t);
    end
    
    i = i + 1; % increase relative index
    
    % update helper variables
    sx = sx + repmat(Xi,[1,ns]);
    sxi = sxi + repmat(Xi,[1,ns]).*i;
    sx2 = sx2 + repmat(Xi*Xi',[1,1,ns]);
    
    if i(1) < smin
        continue % if smallest window is not filled, contiune with next datapoint
    end
    
    for w = 1:ns
        if i(1,w) < s(w)
            break % only iterate through windows that are filled
        else
            % one-pass detrending
            mx = -6*(-2*sxi(:,w)+(s(w)+1)*sx(:,w))./(s(w)*(s(w)^2-1));
            bx = sx(:,w)./s(w) - mx.*(s(w)+1)./2;
            
            % one-pass calculation of detrended covariance
            
            Fw = (((s(w)+1)*(2*s(w)+1))/6).*(mx*mx') + ...
                ((s(w)+1)/2).*((mx*bx') + (mx*bx')') + (bx*bx') + ...
                (1/s(w)).*(squeeze(sx2(:,:,w)) - ...
                ((mx*sxi(:,w)') + (mx*sxi(:,w)')') - ...
                ((bx*sx(:,w)') + (bx*sx(:,w)')'));
                
            
            DCCA_tmp(:,:,w,Mi) = DCCA_tmp(:,:,w,Mi) + (Fw);
            sx(:,w) = 0;
            sxi(:,w) = 0;
            sx2(:,:,w) = 0;
            i(:,w) = 0;
        end
    end
    
    % increase/reset window offset, increase window index
    if wi < smax/smin
        wi = wi + 1;
        continue % if not all windows are filled, continue
    else
        wi = 1;
        i = zeros(nch,ns);
    end
    
    % if 4 windows are filled, calculate fluctuations and 
    if Mi == n_smax
        % save time
        Fi = Fi + 1;
        Fc(Fi).t = t;
        Fc_tmp = denom_tmp.*sum(DCCA_tmp,4);

        % save native DCCA output
        Fc(Fi).DCCA = Fc_tmp;

        % convert to DCCC
        Fc_dccc_tmp = Fc_tmp;
        for j = 1:ns
            DCCA_s = squeeze(Fc_dccc_tmp(:,:,j));
            DFA = sqrt(diag(DCCA_s))*sqrt(diag(DCCA_s)');
            DFA(eye(nch)==1) = sqrt(diag(DFA));
            Fc_dccc_tmp(:,:,j) = DCCA_s./DFA;
        end
        
        % save DCCC output
        Fc(Fi).DCCC = Fc_dccc_tmp;

        % update DCCA_tmp
        DCCA_tmp = DCCA_tmp(:,:,:,2:end); 
        DCCA_tmp = cat(4,DCCA_tmp,zeros(nch,nch,ns,1));
    end
    
    if Mi < n_smax
        Mi = Mi + 1;
    end
end

end

