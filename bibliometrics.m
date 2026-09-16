function R = bibliometrics(C,varargin)
%BIBLIOMETRICS Bibliometric indices, longitudinal profile, diagnostics and plots.
%
%   Syntax
%   ------
%   R = bibliometrics(C)
%   R = bibliometrics(C, Y)
%   R = bibliometrics(C, Y, A)
%   R = bibliometrics(C, Y, A, 'Period', [Y1 Y2])
%   R = bibliometrics(C, Y, A, 'CareerEnd', YEAR)
%   R = bibliometrics(C, Y, A, 'Plots', false)
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
%     - builds a longitudinal publication/citation profile;
%     - returns everything in a structured output R;
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
%   Name-Value pairs
%   ----------------
%   'Period'    : [Y1 Y2] two-element vector of integer years, with Y1 <= Y2.
%                 Restricts the whole analysis to the papers published
%                 inside the closed interval [Y1 Y2]. Requires Y.
%                 Default: [] (no filter).
%
%   'CareerEnd' : scalar integer year. Splits the publication list into a
%                 pre-career-end and a post-career-end phase and reports a
%                 purely descriptive comparison of the two. Requires Y for
%                 the phase analysis. It does NOT redefine the standard
%                 m-index, which keeps using the elapsed career span.
%                 Default: [] (no career-phase analysis).
%
%   'Plots'     : logical scalar. Set to false to suppress every figure and
%                 keep only the Command Window report and the output R.
%                 Useful in batch processing or on very large publication
%                 lists, where the Ferrers diagram becomes expensive.
%                 Default: true.
%
%   Outputs
%   -------
%   R : structure collecting every computed quantity. Results are also
%       printed to the Command Window and visualized in figures.
%
%   R.version                        : version string of the function.
%   R.n                              : number of analysed papers.
%
%   R.citations.total                : total number of citations.
%   R.citations.min                  : minimum citations per paper.
%   R.citations.max                  : maximum citations per paper.
%   R.citations.mode                 : mode of citations per paper.
%   R.citations.median               : median of citations per paper.
%   R.citations.mean                 : mean citations per paper.
%   R.citations.sd                   : standard deviation of citations.
%   R.citations.CV_percent           : coefficient of variation (%).
%   R.citations.CV_adjusted_percent  : small-sample adjusted CV (%).
%   R.citations.Gini                 : Gini's coefficient of inequality.
%   R.citations.h_core_sum           : citations inside the h-core.
%   R.citations.h_core_share         : h-core citations over the total.
%   R.citations.top_g_sum            : citations inside the g-core.
%   R.citations.top_g_share          : g-core citations over the total.
%
%   R.years.first                    : first publication year.
%   R.years.last                     : last publication year.
%   R.years.span_current             : elapsed years since the first paper.
%   R.years.mean_citations_per_year  : mean citations per elapsed year.
%   R.years.unique                   : sorted vector of publication years.
%   R.years.papers_per_year          : papers published in each of them.
%   R.years.cumulative_papers        : cumulative publication count.
%   R.years.cumulative_citations     : cumulative citation count.
%
%   R.indices.h                      : Hirsch's h-index.
%   R.indices.a                      : Hirsch's a parameter (Ctot/h^2).
%   R.indices.chi                    : Fenner's chi-index.
%   R.indices.m                      : Hirsch's m parameter (h/span).
%   R.indices.delta_h                : citations needed to gain one h unit.
%   R.indices.g                      : Egghe's g-index.
%   R.indices.delta_g                : citations needed to gain one g unit.
%   R.indices.A                      : Jin's A-index.
%   R.indices.h2                     : Kosmulski's h2-index.
%   R.indices.e                      : Zhang's e-index.
%   R.indices.R                      : R-index.
%   R.indices.i10                    : i10-index.
%   R.indices.h_normalized           : Sidiropoulos' normalized h-index.
%   R.indices.hc                     : Sidiropoulos' contemporary h-index.
%   R.indices.hc_a                   : a parameter of the hc-index.
%   R.indices.Jin_AWCR               : Jin's age-weighted citation rate.
%   R.indices.Jin_AR                 : Jin's AR-index.
%   R.indices.Harzing_AWCR           : Harzing's age-weighted citation rate.
%   R.indices.Harzing_AR             : Harzing's AR-index.
%   R.indices.Batista_hI             : Batista's individual h-index.
%   R.indices.Harzing_hI_norm        : Harzing's individual h-index.
%   R.indices.Schreiber_hm           : Schreiber's multi-authored h-index.
%   R.indices.Jin_AWCR_author_normalized      : Jin's AWCR per author.
%   R.indices.Jin_AR_author_normalized        : Jin's AR-index per author.
%   R.indices.Harzing_AWCR_author_normalized  : Harzing's AWCR per author.
%   R.indices.Harzing_AR_author_normalized    : Harzing's AR-index per author.
%
%   R.authors                        : [] when A is omitted, otherwise a
%                                      structure with min, max, mode, median,
%                                      mean and citations_per_author.
%
%   R.period.enabled                 : true when 'Period' was used.
%   R.period.range                   : the requested [Y1 Y2] interval.
%
%   R.career.enabled                 : true when 'CareerEnd' was used.
%   R.career.endYear                 : the declared career-end year.
%   R.career.n_pre / n_post          : papers before/after the career end.
%   R.career.citations_pre / _post   : citations before/after it.
%   R.career.citation_share_pre/_post    : the same as fractions of the total.
%   R.career.publication_share_pre/_post : publication fractions.
%   R.career.h_pre / h_post          : h-index of each phase.
%   R.career.first_pre / last_pre    : year range of the first phase.
%   R.career.first_post / last_post  : year range of the second phase.
%   R.career.span_pre                : span of the first phase.
%   R.career.m_pre                   : descriptive h/span of the first phase.
%   R.career.note                    : set when CareerEnd is given without Y.
%
%   R.profile.year                   : publication years actually present.
%   R.profile.papers                 : papers published in each year.
%   R.profile.citations              : citations gathered by those papers.
%   R.profile.cumulativePapers       : running total of publications.
%   R.profile.cumulativeCitations    : running total of citations.
%   R.profile.careerEnd              : career-end year, NaN when unused.
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
%         - Jin's AWCR and AR-index normalized per author (with authors
%           aligned to the h-core ordering)
%         - Harzing's AWCR and AR-index normalized per author
%
%     Career-end descriptive profile (if 'CareerEnd' provided)
%
%   Figures
%   -------
%   Figure 1 : Lorenz curve of citations, with the lines of perfect equality
%              and perfect inequality.
%   Figure 2 : six panels, namely the Durfee square on the Ferrers diagram
%              for the h-index, Kosmulski's h2-index, Egghe's g-index and,
%              when the corresponding inputs are available, Sidiropoulos'
%              hc-index, Harzing's hI,norm-index and Schreiber's hm-index.
%   Figure 3 : cumulative publications and cumulative citations against the
%              publication year, with the career-end marker when requested.
%              Drawn only if Y is provided.
%
%   All figures are suppressed by 'Plots', false, and are skipped anyway
%   when the dataset carries no citations at all.
%
%   Example
%   -------
%   C = [12 8 1 0 5 3 0 0];
%   Y = [2004 2007 2008 2008 2008 2009 2009 2010];
%   A = [8 9 10 7 11 11 7 5];
%   R = bibliometrics(C, Y, A);
%
%   Restricting the analysis to a single decade:
%
%   R = bibliometrics(C, Y, A, 'Period', [2007 2009]);
%
%   Notes
%   -----
%   - The function assumes that C, Y, and A refer to the same ordered list
%     of papers. No attempt is made to disambiguate authors or to retrieve
%     data from the web.
%   - Some indices require at least one cited paper; if all citation counts
%     are zero, several indices become zero or undefined (NaN) and the plots
%     are skipped.
%   - The Gini coefficient includes the origin of the Lorenz curve in the
%     trapezoidal integration, as it must.
%   - The standard m-index is h divided by the elapsed career span. The
%     'CareerEnd' option only adds a descriptive career-phase comparison.
%   - The plotting section can be computationally heavy for very large
%     publication lists, as it builds detailed scatter plots for h-index
%     diagrams. Use 'Plots', false to skip it.
%
%   Changelog
%   ---------
%   2.3.1 - Restored the figures broken by the 2.3.0 refactoring (the
%           plotting section still referred to the old variable names x, x2
%           and cC); fixed the truncated Min/Max line; guarded the g-core
%           sum and the Lorenz curve against zero-citation datasets; added
%           the 'Plots' option; documented the output structure R.
%   2.3.0 - Corrected Gini coefficient; 'Period' filter; 'CareerEnd'
%           descriptive analysis; longitudinal profile; structured output R;
%           safer handling of zero-citation datasets and small h-cores.
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
%   Updated: 2026-09-16
%   Version: 2.3.1
%
%   License
%   -------
%   This function is distributed under the MIT License.
%   See the LICENSE file in the GitHub repository for details.
%

%% Input
p = inputParser;

addRequired(p,'C',@(x) validateattributes(x,{'numeric'}, ...
    {'row','real','finite','nonnan','nonnegative'},mfilename,'C',1));

validationAY = @(x) isempty(x) || ( ...
    isnumeric(x) && isrow(x) && all(isreal(x(:))) && ...
    all(isfinite(x(:))) && ~all(isnan(x(:))) && ...
    all(x(:) > 0) && all(fix(x(:)) == x(:)) );

addOptional(p,'Y',[],validationAY);
addOptional(p,'A',[],validationAY);
addParameter(p,'Period',[],@(x) isempty(x) || ...
    (isnumeric(x) && numel(x)==2 && all(isfinite(x)) && ...
     all(fix(x)==x) && x(1)<=x(2)));
addParameter(p,'CareerEnd',[],@(x) isempty(x) || ...
    (isnumeric(x) && isscalar(x) && isfinite(x) && fix(x)==x));
addParameter(p,'Plots',true,@(x) (islogical(x) || isnumeric(x)) && isscalar(x));

parse(p,C,varargin{:});
C = p.Results.C;
Y = p.Results.Y;
A = p.Results.A;
Period = p.Results.Period;
CareerEnd = p.Results.CareerEnd;
doPlots = logical(p.Results.Plots);
clear p validationAY

n0 = numel(C);

if ~isempty(Y) && numel(Y) ~= n0
    error('bibliometrics:LengthMismatch','C and Y must have the same number of elements.');
end
if ~isempty(A) && numel(A) ~= n0
    error('bibliometrics:LengthMismatch','C and A must have the same number of elements.');
end
if ~isempty(CareerEnd) && ~isempty(Y) && CareerEnd < min(Y)
    error('bibliometrics:CareerEndBeforeData','CareerEnd cannot precede the first publication year.');
end

%% Optional period selection
sel = true(1,n0);
if ~isempty(Period)
    if isempty(Y)
        error('bibliometrics:PeriodNeedsYears','''Period'' requires publication years Y.');
    end
    sel = Y >= Period(1) & Y <= Period(2);
    if ~any(sel)
        error('bibliometrics:EmptyPeriod','No publications fall inside the requested Period.');
    end
end

C = C(sel);
if ~isempty(Y), Y = Y(sel); end
if ~isempty(A), A = A(sel); end
n = numel(C);
currentYear = year(datetime('now'));

%% Descriptive statistics
Ctot = sum(C);
Cmean = mean(C);
Cmedian = median(C);
Cstd = std(C);

if Cmean > 0
    CV = Cstd/Cmean*100;
    CVadj = CV*(1+1/(4*n));
else
    CV = NaN;
    CVadj = NaN;
end

%% Sort citations
[Casc,idxAsc] = sort(C,'ascend');
Csorted = fliplr(Casc);
idx = fliplr(idxAsc);
rank = 1:n;
rank2 = rank.^2;

%% Correct Lorenz/Gini calculation
% The origin (0,0) must be included in the trapezoidal integration.
if Ctot > 0
    F = [0 rank./n];
    L = [0 cumsum(Casc)./Ctot];
    Gcoeff = 1 - 2*trapz(F,L);
else
    F = [0 rank./n];
    L = zeros(1,n+1);
    Gcoeff = NaN;
end

%% Year statistics
hasY = ~isempty(Y);
if hasY
    firstYear = min(Y);
    lastYear = max(Y);
    years = currentYear - firstYear;

    Ny0 = currentYear - Y;
    Ny = Ny0 + 1;

    [uniqueYears,~,idxY] = unique(Y);
    papersPerYear = accumarray(idxY,1);
    meanCitationsPerYear = Ctot/max(years,1);
    mDenominator = years;
else
    firstYear = NaN;
    lastYear = NaN;
    years = NaN;
    Ny = [];
    uniqueYears = [];
    papersPerYear = [];
    meanCitationsPerYear = NaN;
    mDenominator = NaN;
end

%% Core citation indices
Hidx = sum(Csorted >= rank);
H2 = Hidx^2;

if n > 0
    [~,z] = max(Csorted .* rank);
    Chiidx = sqrt(z);
else
    Chiidx = NaN;
end

if Hidx > 0
    HcoreCitations = sum(Csorted(1:Hidx));
    Aidx = HcoreCitations/Hidx;
    Eidx = realsqrt(max(HcoreCitations-H2,0));
    Ridx = realsqrt(HcoreCitations);
else
    HcoreCitations = 0;
    Aidx = NaN;
    Eidx = 0;
    Ridx = 0;
end

cumulativeC = cumsum(Csorted);
Gidx = sum(cumulativeC >= rank2);

if Gidx < n
    dG = (Gidx+1)^2 - cumulativeC(Gidx+1);
else
    dG = NaN;
end

if Hidx < n
    dH = Hidx+1-Csorted(Hidx+1);
else
    dH = NaN;
end

H2idx = sum(Csorted >= rank2);
I10 = sum(C >= 10);
hNorm = Hidx/n;

if hasY
    mIndex = Hidx/mDenominator;
else
    mIndex = NaN;
end

%% Year-weighted indices
Hcidx = NaN; HcA = NaN;
JAWCR = NaN; JAR = NaN;
HAWCR = NaN; HAR = NaN;
Sc = [];

if hasY
    Sc = sort(4 .* Ny.^(-1) .* C,'descend');
    Hcidx = sum(Sc >= rank);
    if Hcidx > 0
        HcA = sum(Sc)/(Hcidx^2);
    end

    Nys = Ny(idx);
    if Hidx > 0
        JAWCR = sum(Csorted(1:Hidx)./Nys(1:Hidx));
    else
        JAWCR = 0;
    end
    JAR = realsqrt(max(JAWCR,0));

    HAWCR = sum(Csorted./Nys);
    HAR = realsqrt(max(HAWCR,0));
end

%% Author-weighted indices
HIidx = NaN; HHIidx = NaN; Hmidx = NaN;
JAWCRN = NaN; JARN = NaN; HAWCRN = NaN; HARN = NaN;
SAH = []; xS = [];

if ~isempty(A)
    AS = A(idx);

    if Hidx > 0
        HIidx = Hidx/mean(AS(1:Hidx));
    end

    SAH = sort(C./A,'descend');
    HHIidx = sum(SAH >= rank);

    xS = cumsum(1./AS);
    validHm = xS <= Csorted;
    if any(validHm)
        Hmidx = xS(find(validHm,1,'last'));
    else
        Hmidx = 0;
    end

    if hasY
        Nys = Ny(idx);
        if Hidx > 0
            JAWCRN = sum(Csorted(1:Hidx)./Nys(1:Hidx)./AS(1:Hidx));
        else
            JAWCRN = 0;
        end
        JARN = realsqrt(max(JAWCRN,0));

        HAWCRN = sum(Csorted./Nys./AS);
        HARN = realsqrt(max(HAWCRN,0));
    end
end

%% CareerEnd descriptive profile
career = struct();
career.enabled = ~isempty(CareerEnd);

if career.enabled
    career.endYear = CareerEnd;

    if hasY
        preMask = Y <= CareerEnd;
        postMask = Y > CareerEnd;

        career.n_pre = sum(preMask);
        career.n_post = sum(postMask);
        career.citations_pre = sum(C(preMask));
        career.citations_post = sum(C(postMask));

        career.citation_share_pre = safeDivide(career.citations_pre,Ctot);
        career.citation_share_post = safeDivide(career.citations_post,Ctot);
        career.publication_share_pre = career.n_pre/n;
        career.publication_share_post = career.n_post/n;

        if any(preMask)
            Cpre = C(preMask);
            career.h_pre = sum(sort(Cpre,'descend') >= (1:numel(Cpre)));
            career.first_pre = min(Y(preMask));
            career.last_pre = max(Y(preMask));
            career.span_pre = career.last_pre-career.first_pre;
            if career.span_pre > 0
                career.m_pre = career.h_pre/career.span_pre;
            else
                career.m_pre = NaN;
            end
        else
            career.h_pre = 0;
            career.first_pre = NaN;
            career.last_pre = NaN;
            career.span_pre = NaN;
            career.m_pre = NaN;
        end

        if any(postMask)
            Cpost = C(postMask);
            career.h_post = sum(sort(Cpost,'descend') >= (1:numel(Cpost)));
            career.first_post = min(Y(postMask));
            career.last_post = max(Y(postMask));
        else
            career.h_post = 0;
            career.first_post = NaN;
            career.last_post = NaN;
        end
    else
        career.note = 'CareerEnd requires Y for phase analysis.';
    end
end

%% Longitudinal profile
% The variable is named prof to avoid shadowing the built-in PROFILE.
prof = struct();
if hasY
    prof.year = uniqueYears(:)';
    prof.papers = papersPerYear(:)';
    prof.citations = zeros(size(prof.year));

    for k = 1:numel(prof.year)
        prof.citations(k) = sum(C(Y == prof.year(k)));
    end

    prof.cumulativePapers = cumsum(prof.papers);
    prof.cumulativeCitations = cumsum(prof.citations);
    prof.careerEnd = CareerEnd;
else
    prof.year = [];
    prof.papers = [];
    prof.citations = [];
    prof.cumulativePapers = [];
    prof.cumulativeCitations = [];
    prof.careerEnd = NaN;
end

%% Structured output
R = struct();
R.version = '2.3.1';
R.n = n;

R.citations.total = Ctot;
R.citations.min = min(C);
R.citations.max = max(C);
R.citations.mode = mode(C);
R.citations.median = Cmedian;
R.citations.mean = Cmean;
R.citations.sd = Cstd;
R.citations.CV_percent = CV;
R.citations.CV_adjusted_percent = CVadj;
R.citations.Gini = Gcoeff;
R.citations.h_core_sum = HcoreCitations;
R.citations.h_core_share = safeDivide(HcoreCitations,Ctot);

if Gidx > 0
    R.citations.top_g_sum = cumulativeC(Gidx);
else
    R.citations.top_g_sum = 0;
end
R.citations.top_g_share = safeDivide(R.citations.top_g_sum,Ctot);

R.years.first = firstYear;
R.years.last = lastYear;
R.years.span_current = years;
R.years.mean_citations_per_year = meanCitationsPerYear;
R.years.unique = uniqueYears;
R.years.papers_per_year = papersPerYear;
R.years.cumulative_papers = prof.cumulativePapers;
R.years.cumulative_citations = prof.cumulativeCitations;

R.indices.h = Hidx;
R.indices.a = safeDivide(Ctot,H2);
R.indices.chi = Chiidx;
R.indices.m = mIndex;
R.indices.delta_h = dH;
R.indices.g = Gidx;
R.indices.delta_g = dG;
R.indices.A = Aidx;
R.indices.h2 = H2idx;
R.indices.e = Eidx;
R.indices.R = Ridx;
R.indices.i10 = I10;
R.indices.h_normalized = hNorm;
R.indices.hc = Hcidx;
R.indices.hc_a = HcA;
R.indices.Jin_AWCR = JAWCR;
R.indices.Jin_AR = JAR;
R.indices.Harzing_AWCR = HAWCR;
R.indices.Harzing_AR = HAR;
R.indices.Batista_hI = HIidx;
R.indices.Harzing_hI_norm = HHIidx;
R.indices.Schreiber_hm = Hmidx;
R.indices.Jin_AWCR_author_normalized = JAWCRN;
R.indices.Jin_AR_author_normalized = JARN;
R.indices.Harzing_AWCR_author_normalized = HAWCRN;
R.indices.Harzing_AR_author_normalized = HARN;

if ~isempty(A)
    R.authors.min = min(A);
    R.authors.max = max(A);
    R.authors.mode = mode(A);
    R.authors.median = median(A);
    R.authors.mean = mean(A);
    R.authors.citations_per_author = sum(C./A);
else
    R.authors = [];
end

R.period = struct('enabled',~isempty(Period),'range',Period);
R.career = career;
R.profile = prof;

%% Command Window output
tr = repmat('-',1,80);

disp('BIBLIOMETRICS');
disp(tr);

if ~isempty(Period)
    fprintf('Analysis period: %i-%i\n',Period(1),Period(2));
end

disp('Descriptive statistics');
disp(tr);
fprintf('Total number of papers: %i\n',n);
fprintf('Total number of citations: %i\n',Ctot);
fprintf('Min: %i - Max: %i\n',min(C),max(C));
fprintf('Mode of citations per paper: %i\n',mode(C));
fprintf('Median of citations per paper: %0.1f\n',Cmedian);
fprintf('Mean number of citations per paper: %0.1f\n',Cmean);

if isnan(CV)
    fprintf('Variation coefficient (CV): N/A (mean = 0)\n');
    fprintf('Adjusted Variation coefficient (CV''''): N/A (mean = 0)\n');
else
    fprintf('Variation coefficient (CV): %0.2f%%\n',CV);
    fprintf('Adjusted Variation coefficient (CV''''): %0.2f%%\n',CVadj);
end

if isnan(Gcoeff)
    fprintf('Gini''s coefficient: N/A (total citations = 0)\n');
else
    fprintf('Gini''s coefficient: %0.2f\n',Gcoeff);
end

if hasY
    disp(tr);
    fprintf('Years: %i\n',years);
    fprintf('Years of publications\t first: %i \t last: %i\n',firstYear,lastYear);
    fprintf('Papers per year\t Min: %i \t Max: %i\n',min(papersPerYear),max(papersPerYear));
    fprintf('Mode of papers per year: %i\n',mode(papersPerYear));
    fprintf('Median of papers per year: %0.1f\n',median(papersPerYear));
    fprintf('Mean number of papers per year: %0.1f\n',mean(papersPerYear));
    fprintf('Mean number of citations per year: %0.1f\n',meanCitationsPerYear);
end

if ~isempty(A)
    disp(tr);
    fprintf('Authors\tMin: %i - Max: %i\n',min(A),max(A));
    fprintf('Mode of Authors per paper: %i\n',mode(A));
    fprintf('Median of Authors per paper: %0.1f\n',median(A));
    fprintf('Mean number of Authors per paper: %0.1f\n',mean(A));
    fprintf('Citations per Author: %0.1f\n',sum(C./A));
end

disp(' ');
disp('Bibliometric indices');
disp(' ');
disp('Citation indices');
disp(tr);

fprintf('Hirsch''s h-index: %i \t a: %0.2f\t',Hidx,R.indices.a);
fprintf('Fenner''s chi-index: %0.4f',Chiidx);
if hasY
    fprintf('\t m: %0.2f',mIndex);
end
if ~isnan(dH)
    fprintf('\tDelta-h: %i\n',dH);
else
    fprintf('\tThis is the max possible h-index\n');
end

fprintf('Egghe''s g-index: %i\t',Gidx);
if ~isnan(dG)
    fprintf('Delta-g: %i\n',dG);
else
    fprintf('This is the max possible g-index\n');
end

fprintf('Jin''s A-index: %0.2f\n',Aidx);
fprintf('Kosmulski''s h2-index: %i\n',H2idx);
fprintf('Zhang''s e-index: %0.1f\n',Eidx);
fprintf('R-index: %0.2f\n',Ridx);
fprintf('i10-index (papers with >= 10 citations): %i\n',I10);
fprintf('Sidiropoulos'' normalized h-index: %0.2f\n',hNorm);

if hasY
    disp(tr);
    disp('Years weighted indices');
    disp(tr);
    fprintf('Sidiropoulos'' Contemporary h-index (hc-index): %i \t a: %0.2f\n',Hcidx,HcA);
    fprintf('Jin''s Age-weighted citation rate (AWCR): %0.2f\n',JAWCR);
    fprintf('Jin''s AR-index: %0.2f\n',JAR);
    fprintf('Harzing''s Age-weighted citation rate (AWCR): %0.2f\n',HAWCR);
    fprintf('Harzing''s AR-index: %0.2f\n',HAR);
end

if ~isempty(A)
    disp(tr);
    disp('Authors weighted indices');
    disp(tr);
    fprintf('Batista''s Individual h-index (hI-index): %0.2f\n',HIidx);
    fprintf('Harzing''s Individual h-index (hI,norm-index): %0.2f\n',HHIidx);
    fprintf('Schreiber''s Multi-authored h-index (hm-index): %0.2f\n',Hmidx);
end

if ~isempty(A) && hasY
    disp(tr);
    disp('Years and Authors weighted indices');
    disp(tr);
    fprintf('Jin''s AWCR normalized per authors: %0.2f\n',JAWCRN);
    fprintf('Jin''s AR-index normalized per authors: %0.2f\n',JARN);
    fprintf('Harzing''s AWCR normalized per authors: %0.2f\n',HAWCRN);
    fprintf('Harzing''s AR-index normalized per authors: %0.2f\n',HARN);
end

if career.enabled
    disp(tr);
    disp('Career-end descriptive profile');
    disp(tr);
    if hasY
        fprintf('Career end: %i\n',CareerEnd);
        fprintf('Publications through career end: %i\n',career.n_pre);
        fprintf('Publications after career end: %i\n',career.n_post);
        fprintf('Citations through career end: %i\n',career.citations_pre);
        fprintf('Citations after career end: %i\n',career.citations_post);
        fprintf('Citation share through career end: %0.1f%%\n',100*career.citation_share_pre);
        fprintf('Citation share after career end: %0.1f%%\n',100*career.citation_share_post);
        fprintf('h-index through career end: %i\n',career.h_pre);
        if ~isnan(career.m_pre)
            fprintf('Descriptive h/span through career end: %0.2f\n',career.m_pre);
        end
    else
        fprintf('%s\n',career.note);
    end
end

%% Figures
if ~doPlots
    return
end

if Ctot == 0
    disp(tr);
    disp('No citations in the dataset: plots skipped.');
    return
end

% -------------------------------------------------------------------------
% Lorenz curve
% -------------------------------------------------------------------------
scrsz = get(groot,'ScreenSize');
hfig1 = figure('Name','Bibliometrics - Lorenz curve');
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

% -------------------------------------------------------------------------
% Index diagrams
% -------------------------------------------------------------------------
hfig2   = figure('Name','Bibliometrics - Index diagrams');
POS(1)  = POS(1) + POS(3);
set(hfig2,'Position',POS)

subplot(2,3,2);
semilogy(rank2,Csorted,'b.',rank2,Csorted,'r-',rank2,rank2,'k-');
axis square
title(sprintf('Kosmulski''s\nh2-index'));
xlabel('Squared Paper Rank');
ylabel('Citations');

subplot(2,3,3);
plot(rank2,cumulativeC,'b.',rank2,cumulativeC,'r-',rank2,rank2,'k-');
axis square
title(sprintf('Egghe''s\ng-index'));
xlabel('Squared Paper Rank');
ylabel('Cumulative sum of Citations');

if hasY
    subplot(2,3,4);
    plot(rank,Sc,'b.',rank,Sc,'r-',rank,rank,'k-');
    axis square
    title(sprintf('Sidiropoulos''s\nhc-index'));
    xlabel('Paper Rank');
    ylabel('Age weighted citations');
end

if ~isempty(A)
    subplot(2,3,5);
    plot(rank,SAH,'b.',rank,SAH,'r-',rank,rank,'k-');
    axis square
    title(sprintf('Harzing''s\nhI,norm-index'));
    xlabel('Paper Rank');
    ylabel('Authors weighted Citations');

    subplot(2,3,6);
    semilogy(xS,Csorted,'b.',xS,Csorted,'r-',xS,xS,'k-');
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
send   = cumsum(Csorted(1:TC));
sstart = [1 send(1:TC-1)+1];

for k = 1:TC
    X(sstart(k):send(k)) = M(1:Csorted(k));
end

Yc = repelem(rank,Csorted)';      % rank i repeated Csorted(i) times
Xs = [X; Yc];
Ys = [Yc; X];

subplot(2,3,1);
scatter(Xs,Ys,4)
axis square
if Hidx > 0
    Xs2 = repmat(1:Hidx,1,Hidx);
    Ys2 = repelem(1:Hidx,repmat(Hidx,1,Hidx));
    hold on
    scatter(Xs2,Ys2,4,'r')
    hold off
end
xlabel('Citations');
ylabel('Citations');
title(sprintf('Hirsch''s h-index as Durfee''s\n square on a Ferrers''es diagram'));

% -------------------------------------------------------------------------
% Longitudinal profile
% -------------------------------------------------------------------------
if hasY
    figure('Name','Bibliometrics - Longitudinal profile');
    yyaxis left
    plot(prof.year,prof.cumulativePapers,'-o','LineWidth',1.2);
    ylabel('Cumulative publications');

    yyaxis right
    plot(prof.year,prof.cumulativeCitations,'-s','LineWidth',1.2);
    ylabel('Cumulative citations');

    xlabel('Publication year');
    title('Longitudinal publication and citation profile');
    grid on;

    if ~isempty(CareerEnd)
        xline(CareerEnd,'--','CareerEnd');
    end
end

end

function y = safeDivide(a,b)
%SAFEDIVIDE Divide while returning NaN for a zero denominator.
if b == 0
    y = NaN;
else
    y = a/b;
end
end