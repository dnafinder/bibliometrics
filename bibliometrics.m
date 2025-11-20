function bibliometrics(C,varargin)
%BIBLIOMETRICS Compute a broad set of bibliometric indices and plots.
%
%   Syntax
%   ------
%   bibliometrics(C)
%   bibliometrics(C, Y)
%   bibliometrics(C, Y, A)
%
%   Description
%   -----------
%   BIBLIOMETRICS computes many of the most commonly used bibliometric
%   indices for a single researcher, given citation counts for their
%   publications and (optionally) publication years and number of authors
%   per paper.
%
%   Typical use case:
%   - You export citation counts (C) from a database such as Web of Science,
%     Scopus, Google Scholar, etc.
%   - You provide publication years (Y) and number of authors (A) from your
%     own records (which avoids name ambiguity and truncated author lists).
%
%   The function:
%     - prints descriptive statistics for citations, years, and authors;
%     - computes a large set of citation-based, age-weighted, and
%       author-weighted indices;
%     - produces several diagnostic plots, including a Lorenz curve and
%       graphical interpretations of h, g, h2, hc, hI,norm, hm.
%
%   Inputs
%   ------
%   C : Row vector of nonnegative numeric values.
%       C(i) is the number of citations of the i-th paper.
%       Type: numeric, row, real, finite, nonnegative.
%       This input is mandatory.
%
%   Y : (optional) Row vector of positive integer years.
%       Y(i) is the publication year of the i-th paper.
%       Must have the same length as C if provided.
%       If omitted or empty ([]), year-weighted indices are skipped.
%
%   A : (optional) Row vector of positive integers.
%       A(i) is the total number of authors on the i-th paper.
%       Must have the same length as C if provided.
%       If omitted or empty ([]), author-weighted indices are skipped.
%
%   All three vectors C, Y, and A must have consistent lengths if given.
%
%   Outputs
%   -------
%   This function does not return variables. All results are:
%     - printed to the Command Window; and
%     - visualized in figures (Lorenz curve and several index-related plots).
%
%   The printed output includes, when possible:
%     Descriptive statistics
%       - Total number of papers and total citations
%       - Min, Max, Mode, Median, Mean citations per paper
%       - Coefficient of variation (CV, and adjusted CV)
%       - Gini's coefficient for citation inequality
%       - Years of activity, first and last year of publication
%       - Min, Max, Mode, Median, Mean papers per year
%       - Mean citations per year
%       - Min, Max, Mode, Median, Mean authors per paper
%       - Citations per author
%
%     Bibliometric indices
%       Citation indices
%         - Hirsch's h-index with a and m parameters, delta-h
%         - Fenner's chi-index
%         - Egghe's g-index, delta-g
%         - Jin's A-index
%         - Kosmulski's h2-index
%         - Zhang's e-index
%         - R-index (root of citation sum in the h-core)
%         - i10-index (number of papers with at least 10 citations)
%         - Sidiropoulos' normalized h-index
%
%       Years-weighted indices (if Y provided)
%         - Sidiropoulos' Contemporary h-index (hc-index)
%         - Jin's age-weighted citation rate (AWCR) and AR-index
%         - Harzing's age-weighted citation rate (AWCR) and AR-index
%
%       Authors-weighted indices (if A provided)
%         - Batista's Individual h-index (hI-index)
%         - Harzing's Individual h-index (hI,norm-index)
%         - Schreiber's Multi-authored h-index (hm-index)
%
%       Years and authors weighted indices (if Y and A provided)
%         - Jin's age-weighted citation rate (AWCR) and AR-index
%           normalized per author (with authors aligned to the h-core)
%         - Harzing's age-weighted citation rate (AWCR) and AR-index
%           normalized per author
%
%   Example
%   -------
%   C = [12 8 1 0 5 3 0 0];
%   Y = [2004 2007 2008 2008 2008 2009 2009 2010];
%   A = [8 9 10 7 11 11 7 5];
%   bibliometrics(C, Y, A)
%
%   Notes
%   -----
%   - The function assumes that C, Y, and A refer to the same ordered list
%     of papers. No attempt is made to disambiguate authors or to retrieve
%     data from the web.
%   - Some indices require at least one cited paper; if all citation counts
%     are zero, several indices become zero or undefined.
%   - The plotting section can be computationally heavy for very large
%     publication lists, as it builds detailed scatter plots for h-index
%     diagrams.
%
%   Citation
%   --------
%   If you use this function in academic or technical work, please cite:
%
%     Cardillo G. (2010)
%     "Bibliometrics: the art of citation indices".
%     Available from GitHub:
%     https://github.com/dnafinder/bibliometrics
%
%   Metadata
%   --------
%   Author : Giuseppe Cardillo
%   Email  : giuseppe.cardillo.75@gmail.com
%   GitHub : https://github.com/dnafinder
%   Created: 2010-01-01
%   Updated: 2025-11-20
%   Version: 2.2.0
%
%   License
%   -------
%   This function is distributed under the MIT License.
%   See the LICENSE file in the GitHub repository for details.
%

% -------------------------------------------------------------------------
% Input error handling
% -------------------------------------------------------------------------
p = inputParser;

addRequired(p,'C',@(x) validateattributes(x,{'numeric'}, ...
    {'row','real','finite','nonnan','nonnegative'},mfilename,'C',1));

validationAY = @(x) isempty(x) || ( ...
    isnumeric(x) && isrow(x) && all(isreal(x(:))) && ...
    all(isfinite(x(:))) && ~all(isnan(x(:))) && ...
    all(x(:) > 0) && all(fix(x(:)) == x(:)) );

addOptional(p,'Y',[],validationAY);
addOptional(p,'A',[],validationAY);

parse(p,C,varargin{:});
C = p.Results.C;
Y = p.Results.Y;
A = p.Results.A;
clear p validationAY

n = numel(C);
if ~isempty(Y) && numel(Y) ~= n
    error('bibliometrics:LengthMismatch', ...
        'C and Y must have the same number of elements.');
end
if ~isempty(A) && numel(A) ~= n
    error('bibliometrics:LengthMismatch', ...
        'C and A must have the same number of elements.');
end

tr = repmat('-',1,80);

disp('BIBLIOMETRICS');
disp(tr)
% -------------------------------------------------------------------------
% Descriptive statistics for citations
% -------------------------------------------------------------------------
disp('Descriptive statistics');
disp(tr)

Ctot = sum(C);
fprintf('Total number of papers: %i\n', n)
fprintf('Total number of citations: %i\n', Ctot)
fprintf('Min: %i -  Max: %i\n', min(C), max(C))
fprintf('Mode of citations per paper: %i\n', mode(C))
fprintf('Median of citations per paper: %0.1f\n', median(C))

M = mean(C);
D = std(C);
if M > 0
    CV = D/M * 100;
    fprintf('Mean number of citations per paper: %0.1f\n', M)
    fprintf('Variation coefficient (CV): %0.2f%%\n', CV)
    fprintf('Adjusted Variation coefficient (CV''''): %0.2f%%\n', CV*(1+1/(4*n)))
else
    fprintf('Mean number of citations per paper: %0.1f\n', M)
    fprintf('Variation coefficient (CV): N/A (mean = 0)\n')
    fprintf('Adjusted Variation coefficient (CV''''): N/A (mean = 0)\n')
end
clear M D CV
disp(' ')

% -------------------------------------------------------------------------
% Lorenz curve and Gini coefficient
% -------------------------------------------------------------------------
[Csorted,idx] = sort(C);
x  = 1:n;
cC = cumsum(Csorted);
F  = x ./ max(x);
L  = cC / Ctot;
Gcoeff = 1 - 2*trapz(F,L);
fprintf('Gini''s coefficient: %0.2f\n', Gcoeff)

disp(tr)

% -------------------------------------------------------------------------
% Year-based statistics (if Y is provided)
% -------------------------------------------------------------------------
if ~isempty(Y)
    currentYear = year(datetime('now'));
    Ny = currentYear - Y;

    [~,~,idxY] = unique(Y);
    cty = sort(accumarray(idxY,1)); % papers per year
    my  = max(Ny);                  % years since first publications

    fprintf('Years: %i\n', my)
    fprintf('Years of publications\t first: %i \t last: %i\n', min(Y), max(Y))
    fprintf('Papers per year\t Min: %i \t Max: %i\n', min(cty), max(cty))
    fprintf('Mode of papers per year: %i\n', mode(cty))
    fprintf('Median of papers per year: %i\n', median(cty))
    fprintf('Mean number of papers per year: %0.1f\n', mean(cty))
    fprintf('Mean number of citations per year: %0.1f\n', Ctot/my)
    disp(tr)
end

% -------------------------------------------------------------------------
% Author-based statistics (if A is provided)
% -------------------------------------------------------------------------
if ~isempty(A)
    fprintf('Authors\tMin: %i -  Max: %i\n', min(A), max(A))
    fprintf('Mode of Authors per paper: %i\n', mode(A))
    fprintf('Median of Authors per paper: %0.1f\n', median(A))
    fprintf('Mean number of Authors per paper: %0.1f\n', mean(A))
    fprintf('Citations per Author: %0.1f\n', sum(C./A))
    disp(tr)
end
disp(' ');

% -------------------------------------------------------------------------
% Bibliometric indices
% -------------------------------------------------------------------------
disp('Bibliometric indices');
disp(' ')
disp('Citations indices');
disp(tr)

% Sort citations in descending order for rank-based indices
Csorted = fliplr(Csorted);
idx     = fliplr(idx);
Hidx    = sum(Csorted >= x);
H2      = Hidx^2;

% Fenner's chi-index
[~,z]  = max(Csorted .* x);
Chiidx = sqrt(z);

fprintf('Hirsch''s h-index: %i \t a: %0.2f\t', Hidx, Ctot/H2)
fprintf('Fenner''s chi-index: %0.4f', Chiidx)

if ~isempty(Y)
    fprintf('\t m: %0.2f', Hidx/my)
else
    fprintf('\n')
end

if Hidx < n
    dH = Hidx+1 - Csorted(Hidx+1);
    fprintf('\tDelta-h: %i\n', dH);
else
    fprintf('\tThis is the max possible h-index\n');
end

x2 = x.^2;
cC = cumsum(Csorted);
Gidx = sum(cC >= x2);
fprintf('Egghe''s g-index: %i\t', Gidx)

if Gidx < n
    dG = (Gidx+1)^2 - sum(Csorted(1:Gidx+1));
    fprintf('Delta-g: %i\n', dG);
else
    fprintf('This is the max possible g-index\n');
end

Aidx = mean(Csorted(1:Hidx));
fprintf('Jin''s A-index: %0.2f\n', Aidx);

H2idx = sum(Csorted >= x2);
fprintf('Kosmulski''s h2-index: %i\n', H2idx);

Eidx = realsqrt(sum(Csorted(1:Hidx)) - H2);
fprintf('Zhang''s e-index: %0.1f\n', Eidx);

% R-index: root of the sum of citations in the h-core
Ridx = realsqrt(sum(Csorted(1:Hidx)));
fprintf('R-index: %0.2f\n', Ridx);

% i10-index: number of papers with at least 10 citations
I10 = sum(C >= 10);
fprintf('i10-index (papers with >= 10 citations): %i\n', I10);

fprintf('Sidiropoulos''es normalized h-index: %0.2f\n', Hidx/n)

disp(tr);
disp(' ');

% -------------------------------------------------------------------------
% Years weighted indices (if Y provided)
% -------------------------------------------------------------------------
if ~isempty(Y)
    disp('Years weighted indices');
    disp(tr)

    Ny  = Ny + 1;                   % avoid division by zero for current-year papers
    Sc  = sort(4 .* Ny.^-1 .* C, 'descend');
    Hcidx = sum(Sc >= x);
    Hc2   = Hcidx^2;
    fprintf('Sidiropoulos''es Contemporary h-index (hc-index): %i \t a: %0.2f\n', ...
        Hcidx, sum(Sc)/Hc2);

    Nys   = Ny(idx);
    JAWCR = sum(Csorted(1:Hidx) ./ Nys(1:Hidx));
    fprintf('Jin''s Age-weighted citation rate (AWCR): %0.2f\n', JAWCR);
    fprintf('Jin''s AR-index (AR-index): %0.2f\n', realsqrt(JAWCR));

    HAWCR = sum(Csorted ./ Nys);
    fprintf('Harzing''s Age-weighted citation rate (AWCR): %0.2f\n', HAWCR);
    fprintf('Harzing''s AR-index (AR-index): %0.2f\n', realsqrt(HAWCR));
    disp(tr);
    disp(' ');
end

% -------------------------------------------------------------------------
% Authors weighted indices (if A provided)
% -------------------------------------------------------------------------
if ~isempty(A)
    disp('Authors weighted indices');
    disp(tr)

    AS = A(idx);

    HIidx = Hidx / mean(AS(1:Hidx));
    fprintf('Batista''s Individual h-index (hI-index): %0.2f\n', HIidx);

    SAH   = sort(C ./ A, 'descend');
    HHIidx = sum(SAH >= x);
    fprintf('Harzing''s Individual h-index (hI,norm-index): %0.2f\n', HHIidx);

    xS  = cumsum(1 ./ AS)';           
    yS  = xS <= Csorted';
    zS  = xS(yS);
    Hmidx = zS(end);
    fprintf('Schreiber''s Multi-authored h-index (hm-index): %0.2f\n', Hmidx);

    disp(tr);
    disp(' ');
end

% -------------------------------------------------------------------------
% Years and Authors weighted indices (if Y and A provided)
% -------------------------------------------------------------------------
if ~isempty(Y) && ~isempty(A)
    disp('Years and Authors weighted indices');
    disp(tr)

    AS  = A(idx);
    Nys = Ny(idx);

    % Jin's per-author AWCR: authors aligned with the h-core ordering
    JAWCRN = sum(Csorted(1:Hidx) ./ Nys(1:Hidx) ./ AS(1:Hidx));
    fprintf(['Jin''s Age-weighted citation rate (AWCR) normalized ', ...
             'per authors: %0.2f\n'], JAWCRN);
    fprintf(['Jin''s AR-index (AR-index) normalized per authors: ', ...
             '%0.2f\n'], realsqrt(JAWCRN));

    HAWCRN = sum(Csorted ./ Nys ./ AS);
    fprintf(['Harzing''s Age-weighted citation rate (AWCR) normalized ', ...
             'per authors: %0.2f\n'], HAWCRN);
    fprintf(['Harzing''s AR-index (AR-index) normalized per authors: ', ...
             '%0.2f\n'], realsqrt(HAWCRN));
    disp(tr)
end

% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------
scrsz = get(groot,'ScreenSize');
hfig1 = figure;
POS   = scrsz;
POS(3)= POS(3)/2;
set(hfig1,'Position',POS)
hold on
patch([0 1 1 0],[0 1 0 0],[192 192 192]./255)
patch([0 F 1 0],[0 L 0 0],'w')
Le1 = plot([0 1],[0 0],'g','LineWidth',2);
plot([1 1],[0 1],'g','LineWidth',2)
Le2 = plot([0 1],[0 1],'b--','LineWidth',2);
Le3 = plot(F,L,'r-','LineWidth',2);
hold off
title('Lorenz curve of citations');
xlabel('% of papers');
ylabel('% of citations');
legend([Le1 Le2 Le3], ...
    'Line of perfect inequality', ...
    'Line of perfect equality', ...
    'Lorenz curve', ...
    'Location','NorthEastOutside')
axis square

hfig2   = figure;
POS(1)  = POS(1) + POS(3);
set(hfig2,'Position',POS)
subplot(2,3,2);
plot(x2,Csorted,'b.',x2,Csorted,'r-',x2,x2,'k-');
axis square
title(sprintf('Kosmulski''s\nh2-index'));
xlabel('Squared Paper Rank');
ylabel('Citations');

subplot(2,3,3);
plot(x2,cC,'b.',x2,cC,'r-',x2,x2,'k-');
axis square
title(sprintf('Egghe''s\ng-index'));
xlabel('Squared Paper Rank');
ylabel('Cumulative sum of Citations');

if ~isempty(Y)
    subplot(2,3,4);
    plot(x,Sc,'b.',x,Sc,'r-',x,x,'k-');
    axis square
    title(sprintf('Sidiropoulos''s\nhc-index'));
    xlabel('Paper Rank');
    ylabel('Age weighted citations');
end

if ~isempty(A)
    subplot(2,3,5);
    plot(x,SAH,'b.',x,SAH,'r-',x,x,'k-');
    axis square
    title(sprintf('Harzing''s\nhI,norm-index'));
    xlabel('Paper Rank');
    ylabel('Authors weighted Citations');

    subplot(2,3,6);
    plot(xS,Csorted,'b.',xS,Csorted,'r-',xS,xS,'k-');
    axis square
    title(sprintf('Schreiber''s\nhm-index'));
    xlabel('Authors weighted Paper Rank');
    ylabel('Citations');
end

% -------------------------------------------------------------------------
% Ferrers diagram for Hirsch's h-index (Durfee square)
% -------------------------------------------------------------------------
TC = sum(Csorted > 0);
X  = zeros(Ctot,1);
M  = 1:Csorted(1);
send  = cumsum(Csorted(1:TC));
sstart = [1 send(1:TC-1)+1];

for k = 1:TC
    X(sstart(k):send(k)) = M(1:Csorted(k));
end

Yc = repelem(x, Csorted)';      % rank i repeated Csorted(i) times
Xs = [X; Yc];
Ys = [Yc; X];

subplot(2,3,1);
scatter(Xs,Ys,4)
Xs2 = repmat(1:Hidx,1,Hidx);
Ys2 = repelem(1:Hidx, repmat(Hidx,1,Hidx));
hold on
scatter(Xs2,Ys2,4,'r')
axis square
hold off
xlabel('Citations');
ylabel('Citations');
title(sprintf('Hirsch''s h-index as Durfee''s\n square on a Ferrers''es diagram'));

end
