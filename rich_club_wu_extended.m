function   [Rw, Rw_numerator, Rw_denominator] = rich_club_wu_extended(CIJ,varargin)
%RICH_CLUB_WU 	Rich club coefficients curve (weighted undirected graph)
%
%   Rw = rich_club_wu(CIJ,varargin) % rich club curve for weighted graph
%
%   The weighted rich club coefficient, Rw, at level k is the fraction of
%   edge weights that connect nodes of degree k or higher out of the
%   maximum edge weights that such nodes might share.
%
%   Inputs:
%       CIJ:        weighted directed connection matrix
%
%       k-level:    (optional) max level of RC(k).
%                   (by default k-level quals the maximal degree of CIJ)
%
%   Output:
%       Rw:         rich-club curve
%
%
%   References:
%       T Opsahl et al. Phys Rev Lett, 2008, 101(16)
%       M van den Heuvel, O Sporns, J Neurosci 2011 31(44)
%
%   Martijn van den Heuvel, University Medical Center Utrecht, 2011

%   Modification History:
%   2011: Original
%   2015: Expanded documentation (Mika Rubinov)

NofNodes = size(CIJ,2);     %#ok<NASGU> %number of nodes
NodeDegree = sum((CIJ~=0)); %define degree of each node

%define to which level rc should be computed
if size(varargin,2)==1
    klevel = varargin{1};
elseif isempty(varargin)
    klevel = max(NodeDegree);
else
    error('number of inputs incorrect. Should be [CIJ], or [CIJ, klevel]')
end

%wrank contains the ranked weights of the network, with strongest connections on top
wrank = sort(CIJ(:), 'descend');

%loop over all possible k-levels
%[Rw, Rw_numerator, Rw_denominator, Er] = deal(nan(klevel, 1));
%[Rw, Rw_numerator, Rw_denominator] = deal(nan(klevel, 1));
for kk = 1:klevel
    
    SmallNodes=find(NodeDegree<kk);
    
    if isempty(SmallNodes)
        [Rw(kk), Rw_numerator(kk), Rw_denominator(kk)] = deal(NaN);             %#ok<*AGROW>
        continue
    end
    
    %remove small nodes with NodeDegree<kk
    CutoutCIJ=CIJ;
    CutoutCIJ(SmallNodes,:)=[];
    CutoutCIJ(:,SmallNodes)=[];
    
    %total weight of connections in subset E>r
    Rw_numerator(kk) = sum(CutoutCIJ(:));
    
    %total number of connections in subset E>r
%     Er(kk) = length(find(CutoutCIJ~=0)); % sum((CutoutCIJ~=0),[1,2]) would be faster?
%     fprintf('%d\t', length(find(CutoutCIJ~=0)))
%     if ~mod(kk, 10)
%         fprintf('\n')
%     end
    Er = length(find(CutoutCIJ~=0));
    
    %E>r number of connections with max weight in network
%     wrank_r = wrank(1:1:Er(kk));
    wrank_r = wrank(1:1:Er);
    Rw_denominator(kk) = sum(wrank_r);
    
    %weighted rich-club coefficient
    Rw(kk) = Rw_numerator(kk) / Rw_denominator(kk);
end
