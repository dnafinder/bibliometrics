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
%   R = bibliometrics(C, Y, A, 'Interpret', false)
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
%       author-weighted indices, plus a second tier of extended indices
%       (hg, w, maxprod, o, q2, tapered h);
%     - builds a longitudinal publication/citation profile;
%     - returns everything in a structured output R;
%     - optionally prints four plain-language interpretation flags
%       (sustainedness, concentration, breadth vs depth, recent activity);
%     - produces three figures: an overview (Lorenz curve and longitudinal
%       profile), the classic index diagrams, and the extended-index
%       diagrams.
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
%   'Interpret' : logical scalar. Set to false to suppress the plain-
%                 language "Interpretation flags" block from the printed
%                 report. R.flags is still populated either way. The
%                 flags are orientation heuristics, not judgements: each
%                 compares the profile against itself or against a
%                 published rule of thumb (cited inline), never against a
%                 single external benchmark.
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
%   R.indices.Alonso_hg              : Alonso et al.'s hg-index.
%   R.indices.Woeginger_w            : Woeginger's w-index.
%   R.indices.Kosmulski_maxprod      : Kosmulski's maxprod-index.
%   R.indices.DortaGonzalez_o        : Dorta-Gonzalez's o-index.
%   R.indices.Cabrerizo_q2           : Cabrerizo et al.'s q2-index.
%   R.indices.Anderson_ht            : Anderson et al.'s tapered h-index.
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
%   R.flags.sustainedness            : level/symbol/value/message on the
%                                       m-index (pace of output).
%   R.flags.concentration            : level/symbol/value/message on the
%                                       Gini coefficient.
%   R.flags.breadth                  : level/symbol/value/message on the
%                                       g/h ratio.
%   R.flags.recent_activity          : level/symbol/value/message on the
%                                       last five years' share of papers
%                                       and citations (never coloured).
%   Each flags.* field has:
%     .level    'green'/'yellow'/'red'/'info'/'n/a'
%     .symbol   a plain-text tag ('[GREEN]'/'[YELLOW]'/'[RED]'/'[INFO]'),
%               so the report renders identically on every terminal and
%               font, with no dependency on colour-emoji glyph support
%     .value    the underlying number (or a small struct for
%               recent_activity)
%     .message  the full plain-language sentence
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
%       Extended indices (always printed; require only C)
%         - Alonso et al.'s hg-index (geometric mean of h and g)
%         - Woeginger's w-index
%         - Kosmulski's maxprod-index
%         - Dorta-Gonzalez's o-index
%         - Cabrerizo et al.'s q2-index
%         - Anderson et al.'s tapered h-index (ht)
%
%       Interpretation flags (suppressed by 'Interpret', false)
%         - Sustainedness, concentration, breadth vs depth, recent
%           activity - see R.flags above for details
%
%     Career-end descriptive profile (if 'CareerEnd' provided)
%
%   Figures
%   -------
%   Figure 1 : Overview, two panels - the Lorenz curve of citations (with
%              the lines of perfect equality and perfect inequality, and
%              the Gini coefficient in the title) and the longitudinal
%              cumulative publications/citations profile against the
%              publication year (right panel skipped if Y is not
%              provided).
%   Figure 2 : six panels, namely the Durfee square on the Ferrers diagram
%              for the h-index, Kosmulski's h2-index, Egghe's g-index and,
%              when the corresponding inputs are available, Sidiropoulos'
%              hc-index, Harzing's hI,norm-index and Schreiber's hm-index.
%   Figure 3 : six panels for the extended indices - the w-index triangle
%              threshold, the maxprod-index curve, a bar comparison for
%              hg, a bar comparison for o, the q2-index h-core with its
%              median line, and the cumulative tapered-h curve.
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
%   2.4.1 - Switched the interpretation-flag symbols from Unicode colour
%           emoji to plain-text tags ('[GREEN]'/'[YELLOW]'/'[RED]'/
%           '[INFO]'): emoji rendering in the Command Window depends on
%           the desktop/terminal font and colour-emoji glyph support,
%           which is inconsistent across platforms, while text tags are
%           font-independent.
%   2.4.0 - Added six extended indices (hg, w, maxprod, o, q2, tapered h),
%           reported both in the printed report and in R.indices under
%           the Author_index naming convention (e.g. Woeginger_w) already
%           used for the disambiguated indices; added the 'Interpret'
%           option and four plain-language interpretation flags
%           (sustainedness, concentration, breadth vs depth, recent
%           activity), reported in the printed output and in R.flags,
%           each with a Unicode traffic-light symbol; restructured the
%           figures into three: an overview (Lorenz curve + longitudinal
%           profile, merged from the former separate figures), the
%           unchanged six-panel index diagrams, and a new six-panel
%           extended-index diagram.
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
%   Updated: 2026-09-17
%   Version: 2.4.1
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
addParameter(p,'Interpret',true,@(x) (islogical(x) || isnumeric(x)) && isscalar(x));

parse(p,C,varargin{:});
C = p.Results.C;
Y = p.Results.Y;
A = p.Results.A;
Period = p.Results.Period;
CareerEnd = p.Results.CareerEnd;
doPlots = logical(p.Results.Plots);
doInterpret = logical(p.Results.Interpret);
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

%% Extended indices
% A second tier of citation indices, kept apart from the classic block
% above so that the standard report stays uncluttered. All are computed
% from C alone and therefore never require Y or A.
HGidx = realsqrt(Hidx*Gidx);
MaxProdidx = 0;
Oidx = 0;
Q2idx = 0;
HTidx = 0;
Widx = 0;

if n > 0
    MaxProdidx = max(rank .* Csorted);
    if Hidx > 0
        Oidx = realsqrt(Hidx*Csorted(1));
        Q2idx = realsqrt(Hidx*median(Csorted(1:Hidx)));
    end
    tCsorted = min(Csorted,rank);
    HTidx = sum(tCsorted./rank);

    % Woeginger's w-index: the largest w such that the i-th most-cited
    % paper has at least (w-i+1) citations for every i = 1..w. Solved in
    % closed form (no explicit loop) via a cumulative minimum: U(i) is
    % the largest w admissible by paper i alone, and w must not exceed
    % the running minimum of U up to that rank.
    U = Csorted + rank - 1;
    CM = cummin(U);
    lastValid = find(CM >= rank,1,'last');
    if ~isempty(lastValid)
        Widx = lastValid;
    end
end

%% Interpretation flags
% Four descriptive signals that translate the numeric indices into plain
% language, computed once here and reused for both the printed report
% and the R.flags output. These are orientation heuristics, not
% judgements: each compares the profile against itself, or against a
% published rule of thumb cited inline, never against a single external
% benchmark. Symbols are plain-text tags ('[GREEN]' etc.), not emoji:
% Command Window rendering of colour emoji depends on the desktop/
% terminal font, which is inconsistent across platforms, while a text
% tag always renders correctly. Suppressed entirely when 'Interpret' is
% false, but still computed for R.flags unless Plots/Interpret logic
% says otherwise - here they are always computed (cheap) and only the
% report is gated.
flags = struct();

% --- Sustainedness: pace of output, h relative to elapsed years ---
if hasY
    if mIndex >= 1
        flags.sustainedness.level = 'green';
        flags.sustainedness.symbol = '[GREEN]';
    elseif mIndex >= 0.5
        flags.sustainedness.level = 'yellow';
        flags.sustainedness.symbol = '[YELLOW]';
    else
        flags.sustainedness.level = 'red';
        flags.sustainedness.symbol = '[RED]';
    end
    flags.sustainedness.value = mIndex;
    flags.sustainedness.message = sprintf(['Sustainedness: m = %0.2f (h divided by elapsed years) - %s. ' ...
        'Hirsch (2005) suggested m~1 as typical of a continuously active scientist; ' ...
        'this is a rule of thumb, not a pass/fail threshold.'], mIndex, flags.sustainedness.level);
else
    flags.sustainedness.level = 'n/a';
    flags.sustainedness.symbol = '';
    flags.sustainedness.value = NaN;
    flags.sustainedness.message = 'Sustainedness: not available (requires Y).';
end

% --- Concentration: Gini coefficient of the citation distribution ---
if ~isnan(Gcoeff)
    if Gcoeff < 0.35
        flags.concentration.level = 'green';
        flags.concentration.symbol = '[GREEN]';
        giniWord = 'fairly even across the publication set';
    elseif Gcoeff < 0.55
        flags.concentration.level = 'yellow';
        flags.concentration.symbol = '[YELLOW]';
        giniWord = 'moderately concentrated: a handful of papers carry more weight than the rest';
    else
        flags.concentration.level = 'red';
        flags.concentration.symbol = '[RED]';
        giniWord = 'strongly concentrated: a small number of papers dominate the citation total';
    end
    flags.concentration.value = Gcoeff;
    flags.concentration.message = sprintf(['Concentration: Gini = %0.2f - citations are %s ' ...
        '(0 = perfectly even, 1 = a single paper holds everything).'], Gcoeff, giniWord);
else
    flags.concentration.level = 'n/a';
    flags.concentration.symbol = '';
    flags.concentration.value = NaN;
    flags.concentration.message = 'Concentration: not available (total citations = 0).';
end

% --- Breadth vs depth: g/h ratio ---
if Hidx > 0
    ghRatio = Gidx/Hidx;
    if ghRatio < 1.1
        flags.breadth.level = 'yellow';
        flags.breadth.symbol = '[YELLOW]';
        breadthWord = 'a tightly concentrated profile: few papers carry almost all of the measurable impact';
    elseif ghRatio <= 2
        flags.breadth.level = 'green';
        flags.breadth.symbol = '[GREEN]';
        breadthWord = 'a balanced profile between a solid citation core and a longer tail';
    elseif ghRatio <= 3
        flags.breadth.level = 'yellow';
        flags.breadth.symbol = '[YELLOW]';
        breadthWord = 'a broad, long-tailed profile: many papers sit below the h-core';
    else
        flags.breadth.level = 'red';
        flags.breadth.symbol = '[RED]';
        breadthWord = 'an extremely long-tailed profile: the h-core is a small fraction of total output';
    end
    flags.breadth.value = ghRatio;
    flags.breadth.message = sprintf('Breadth vs depth: g/h = %0.2f - %s.',ghRatio,breadthWord);
else
    flags.breadth.level = 'n/a';
    flags.breadth.symbol = '';
    flags.breadth.value = NaN;
    flags.breadth.message = 'Breadth vs depth: not available (h = 0).';
end

% --- Recent activity: purely descriptive, deliberately uncoloured ---
% A young paper's low citation count is a fact of its age, not a
% shortcoming, so this signal never gets a red/yellow/green label.
if hasY
    recentMask = Y >= (currentYear-4);
    nRecent = sum(recentMask);
    cRecent = sum(C(recentMask));
    pctPapersRecent = 100*nRecent/n;
    pctCitRecent = safeDivide(100*cRecent,Ctot);
    flags.recent_activity.level = 'info';
    flags.recent_activity.symbol = '[INFO]';
    flags.recent_activity.value = struct('n_recent',nRecent,'papers_percent',pctPapersRecent, ...
        'citations_percent',pctCitRecent);
    if isnan(pctCitRecent)
        flags.recent_activity.message = sprintf(['Recent activity: %i of %i papers (%0.0f%%) were ' ...
            'published in the last 5 years.'],nRecent,n,pctPapersRecent);
    else
        flags.recent_activity.message = sprintf(['Recent activity: %i of %i papers (%0.0f%%) were ' ...
            'published in the last 5 years, holding %0.0f%% of total citations - recent papers have had ' ...
            'less time to accumulate citations, so a lower share here is expected, not a warning sign.'], ...
            nRecent,n,pctPapersRecent,pctCitRecent);
    end
else
    flags.recent_activity.level = 'n/a';
    flags.recent_activity.symbol = '';
    flags.recent_activity.value = struct();
    flags.recent_activity.message = 'Recent activity: not available (requires Y).';
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
R.version = '2.4.1';
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

R.indices.Alonso_hg = HGidx;
R.indices.Woeginger_w = Widx;
R.indices.Kosmulski_maxprod = MaxProdidx;
R.indices.DortaGonzalez_o = Oidx;
R.indices.Cabrerizo_q2 = Q2idx;
R.indices.Anderson_ht = HTidx;

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
R.flags = flags;

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

disp(' ');
disp('Extended indices');
disp(tr);
fprintf('Alonso''s hg-index: %0.2f\n',HGidx);
fprintf('Woeginger''s w-index: %i\n',Widx);
fprintf('Kosmulski''s maxprod-index: %i\n',MaxProdidx);
fprintf('Dorta-Gonzalez''s o-index: %0.2f\n',Oidx);
fprintf('Cabrerizo''s q2-index: %0.2f\n',Q2idx);
fprintf('Anderson''s tapered h-index (ht): %0.2f\n',HTidx);

if doInterpret
    disp(' ');
    disp('Interpretation flags');
    disp(tr);
    fprintf('%s %s\n',flags.sustainedness.symbol,flags.sustainedness.message);
    fprintf('%s %s\n',flags.concentration.symbol,flags.concentration.message);
    fprintf('%s %s\n',flags.breadth.symbol,flags.breadth.message);
    fprintf('%s %s\n',flags.recent_activity.symbol,flags.recent_activity.message);
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

scrsz = get(groot,'ScreenSize');
POS   = scrsz;
POS(3)= POS(3)/3;

% -------------------------------------------------------------------------
% Figure 1: Overview - Lorenz curve and longitudinal profile
% -------------------------------------------------------------------------
hfig1 = figure('Name','Bibliometrics - Overview');
set(hfig1,'Position',POS)

subplot(1,2,1);
hold on
patch([0 1 1 0],[0 1 0 0],[192 192 192]./255)
patch([0 F 1 0],[0 L 0 0],'w')
Le1 = plot([0 1],[0 0],'g','LineWidth',2);
plot([1 1],[0 1],'g','LineWidth',2)
Le2 = plot([0 1],[0 1],'b--','LineWidth',2);
Le3 = plot(F,L,'r-','LineWidth',2);
hold off
title(sprintf('Lorenz curve of citations\n(Gini = %0.2f)',Gcoeff));
xlabel('% of papers');
ylabel('% of citations');
legend([Le1 Le2 Le3], ...
    'Line of perfect inequality', ...
    'Line of perfect equality', ...
    'Lorenz curve', ...
    'Location','SouthOutside')
axis square

subplot(1,2,2);
if hasY
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
else
    axis off
    text(0.5,0.5,'Longitudinal profile requires Y', ...
        'HorizontalAlignment','center','Units','normalized');
end

% -------------------------------------------------------------------------
% Figure 2: Index diagrams
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
% Figure 3: Extended indices
% -------------------------------------------------------------------------
hfig3   = figure('Name','Bibliometrics - Extended indices');
POS(1)  = POS(1) + POS(3);
set(hfig3,'Position',POS)

subplot(2,3,1);
stairs(rank,Csorted,'b-','LineWidth',1.2);
hold on
if Widx > 0
    triX = 1:Widx;
    triY = Widx - triX + 1;
    stairs(triX,triY,'r--','LineWidth',1.5);
end
hold off
axis square
title(sprintf('Woeginger''s\nw-index'));
xlabel('Paper Rank');
ylabel('Citations');

subplot(2,3,2);
prodVec = rank.*Csorted;
plot(rank,prodVec,'b.-');
hold on
[~,imaxprod] = max(prodVec);
plot(rank(imaxprod),prodVec(imaxprod),'ro','MarkerFaceColor','r');
hold off
axis square
title(sprintf('Kosmulski''s\nmaxprod-index'));
xlabel('Paper Rank');
ylabel('Rank x Citations');

subplot(2,3,3);
bar(categorical({'h','g','hg'},{'h','g','hg'}),[Hidx Gidx HGidx]);
axis square
title(sprintf('Alonso''s\nhg-index'));
ylabel('Value');

subplot(2,3,4);
bar(categorical({'h','C_{max}','o'},{'h','C_{max}','o'}),[Hidx Csorted(1) Oidx]);
axis square
title(sprintf('Dorta-Gonzalez''s\no-index'));
ylabel('Value');

subplot(2,3,5);
if Hidx > 0
    bar(1:Hidx,Csorted(1:Hidx));
    hold on
    yline(median(Csorted(1:Hidx)),'r--','LineWidth',1.5);
    hold off
end
axis square
title(sprintf('Cabrerizo''s\nq2-index'));
xlabel('Paper Rank (h-core)');
ylabel('Citations');

subplot(2,3,6);
cumTapered = cumsum(tCsorted./rank);
plot(rank,cumTapered,'b.-','LineWidth',1.2);
axis square
title(sprintf('Anderson''s\ntapered h-index'));
xlabel('Paper Rank');
ylabel('Cumulative tapered score');

end

function y = safeDivide(a,b)
%SAFEDIVIDE Divide while returning NaN for a zero denominator.
if b == 0
    y = NaN;
else
    y = a/b;
end
end
