function [pk,ps,spk,sps] = SLL(AF,plt,par)
% SLL       [pk,ps,spk,sps] = SLL(AF,plt,par)
%           Calculate maximum SLL and maximum shoulder for a linear antenna pattern
%           Assumes beam peak is clearly defined and 'enough' angles given
%
%           AF  = normalized pattern [th(Nth),Directivity(NthxNph)].
%           plt = [] 1 for plots (default = 0)
%           par = [No of points to use in polynomial fitting (odd),
%                  Polynomial order
%                  No of points to generate to determine maximum - this defines accuracy]
%                 (default = [7 4 1000])
%
%           pk  = [dB]  maximum (interpolated) SLL
%           ps  = [deg] angle at maximum SLL
%           spk = [dB]  maximum shoulder
%           sps = [deg] angle at maximum shoulder

%           Author: R. Lehmensiek, 08/2002.
%           Latest change: DIL de Villiers 08/2019

if nargin==1, plt = []; par = [];
elseif nargin==2, par = [];
end

if isempty(plt), plt = 0; end
if isempty(par), par = [7 4 1000]; end

pts = par(1);     % [7] No of points to use in polynomial fitting (odd)
ord = par(2);     % [4]  Polynomial order
Nop = par(3);     % [1000] No of points to generate to determine maximum - this defines accuracy
Nos = 2;          % [2] No of SLL peaks to investigate in order to determine maximum

N = size(AF,2);
[pk,ps] = deal(NaN(N-1,Nos));

warning off MATLAB:nearlySingularMatrix
if plt, hf = figure; hl = zeros(1,2+Nos*3); end

for k = 2:N
    % find the maximum sll
    pat = AF(:,k);
    dpat = gradient(pat);
    
    ind = [find(dpat>=0); length(dpat)];
    ind1 = find(diff(ind)~=1);
    
    pk_ = pat(ind(ind1));
    ps_ = ind(ind1);
    
    % If only last point found - ignore it
    if ~isempty(ps_) && ps_(end)==length(pat)
        pk_ = pk_(1:end-1); ps_ = ps_(1:end-1);
    end
    if isempty(ps_)
        % Do nothing - TODO: return the shoulder?
    else
        if (length(ps_))<Nos, Nos = length(ps_); end
        
        if plt
            clf;
            hl(1:2) = plot(AF(:,1),pat,'b-',AF(ps_,1),pk_,'ro'); hold on;
            axis([min(AF(:,1)) max(AF(:,1)) -50 0]);
        end
        
        [pk_,ind] = sort(pk_);
        pk_ind = ps_(ind(end));
        pk_ = pk_(end-Nos+1:end);
        ps_ = ps_(ind(end-Nos+1:end));
        if plt
            hl(3) = plot(AF(ps_,1),pk_,'g*');
            pause;
            delete(hl(2:3));
        end
        
        hpts = fix((pts-1)/2);
        
        for n = 1:Nos
            ind = [ps_(n)-hpts:ps_(n)+hpts];
            ind = ind(ind<=size(AF,1) & ind>0);
            ang1 = [AF(ind(1),1):(AF(ind(end),1)-AF(ind(1),1))/(Nop-1):AF(ind(end),1)];
            [P,~,MU] = polyfit(AF(ind,1),pat(ind),ord);
            Y = polyval(P,(ang1-MU(1))./MU(2));
            [pk_(n),tind] = max(Y);
            ps_(n) = ang1(tind);
            if plt
                hl(2+(n-1)*3:n*3+1) = plot(AF(ind,1),pat(ind),'b*',ang1,Y,'r-',ps_(n),pk_(n),'g*');
                pause
            end
        end
        
        pk(k-1,1:Nos) = flipud(pk_).';
        ps(k-1,1:Nos) = flipud(ps_).';
        if plt
            hl(end) = plot(ps,pk,'mo','Markersize',10);
            pause
            ind = find(hl);
            delete(hl(ind(2:end)));
        end
    end
    
    % TODO: This part not checked for compatibility
    if nargout>2
        ddpat = gradient(dpat);
        ind1 = pk_ind-1; while sign(dpat(ind1))==sign(dpat(ind1-1)), if ind1>2, ind1 = ind1-1; else ind1 = 1; end; end
        ind2 = pk_ind+1; while sign(dpat(ind2))==sign(dpat(ind2+1)), if ind2<length(dpat)-1, ind2 = ind2+1; else ind2 = length(ddpat); end; end
        ind = ind1:ind2;
        %     plot(AF(ind1:ind2,1),10000*ddpat(ind1:ind2)); pause
        
        sps_ = [];
        for k = ind1+1:pk_ind
            if (sign(ddpat(k))~=sign(ddpat(k-1))) & (ddpat(k-1)<ddpat(k)), sps_ = [sps_ k]; end
        end
        for k = pk_ind+1:ind2
            if (sign(ddpat(k))~=sign(ddpat(k-1))) & (ddpat(k-1)>ddpat(k)), sps_ = [sps_ k]; end
        end
        
        if length(sps_)==0, [spk,sps] = deal(NaN);
        else
            if length(sps_)<Nos, Nos = length(sps_); end
            
            spk_ = pat(sps_);
            if plt, hl(2) = plot(AF(sps_,1),spk_,'ro'); end
            
            [spk_,ind] = sort(spk_);
            spk_ = spk_(end-Nos+1:end);
            sps_ = sps_(ind(end-Nos+1:end));
            if plt
                hl(3) = plot(AF(sps_,1),spk_,'g*');
                pause
                delete(hl(2:3));
            end
            
            for n = 1:Nos
                ind = [sps_(n)-hpts:sps_(n)+hpts];
                ang1 = [AF(ind(1),1):(AF(ind(end),1)-AF(ind(1),1))/(Nop-1):AF(ind(end),1)];
                [t,Y] = polfit(AF(ind,1),dpat(ind),ord,ang1);
                [spk_(n),tind] = min(abs(Y));
                sps_(n) = ang1(tind);
                [t,spk_(n)] = polfit(AF(ind,1),pat(ind),ord,sps_(n));
                if plt, plot(AF(ind,1),pat(ind),'b*',[ang1(1) ang1(1)],[-50 0],'r-',[ang1(end) ang1(end)],[-50 0],'r-',sps_(n),spk_(n),'g*'); end
            end
            
            [spk,ind] = max(spk_);
            sps = sps_(ind);
            if plt
                plot(sps,spk,'mo','Markersize',10);
                pause
            end
        end
    end
end

if plt, delete(hf); end
warning on MATLAB:nearlySingularMatrix
end