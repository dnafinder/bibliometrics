[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/bibliometrics)

📌 Overview This repository provides the MATLAB function bibliometrics, which computes a broad range of bibliometric indicators for a single researcher. Given the citation counts of a set of publications (and optionally the publication years and number of authors per paper), the function reports descriptive statistics, many well-known citation indices, plain-language interpretation flags, and several diagnostic plots. The goal is to reproduce and extend the type of analysis provided by tools such as Publish or Perish or Google Scholar, but without relying on online data mining and without the ambiguities of author name matching.

✨ Features The function bibliometrics:

Computes descriptive statistics for citations, years, and authors (when provided).
Calculates a large number of bibliometric indices, including h, g, h2, e, A, R-index, i10-index, contemporary h, age-weighted citation rates, individual h-index variants, and multi-authored indices.
Calculates a second tier of extended indices: Alonso et al.'s hg-index, Woeginger's w-index, Kosmulski's maxprod-index, Dorta-Gonzalez's o-index, Cabrerizo et al.'s q2-index, and Anderson et al.'s tapered h-index.
Prints four plain-language interpretation flags (sustainedness, concentration, breadth vs depth, recent activity), each tagged [GREEN]/[YELLOW]/[RED]/[INFO] in the Command Window and mirrored in the output structure - orientation heuristics, not judgements, each cited against a published rule of thumb rather than a single external benchmark.
Provides inequality analysis via the Lorenz curve and Gini coefficient.
Generates three figures: an overview (Lorenz curve and longitudinal profile), the classic index diagrams (h, g, h2, hc, hI,norm, hm), and the extended-index diagrams (w, maxprod, hg, o, q2, tapered h). All computations are done locally on the data you supply, so there is no dependence on external web services or online databases.

🛠 Installation Download or clone this repository from GitHub: https://github.com/dnafinder/bibliometrics

Add the folder containing bibliometrics.m to your MATLAB path using the Add Folder to Path option or the addpath command.

⚙️ Requirements The function uses only core MATLAB functionality (graphics, basic statistics, and datetime) and does not require additional toolboxes. Unicode emoji symbols are not used anywhere in the code, so the Command Window output renders identically regardless of font or terminal.

▶️ Usage At minimum, you must provide a row vector of citation counts for each of your publications. Optionally, you may also provide publication years and number of authors per paper in the same order.

Examples: C = [12 8 1 0 5 3 0 0]; Y = [2004 2007 2008 2008 2008 2009 2009 2010]; A = [8 9 10 7 11 11 7 5];

bibliometrics(C) % only citation-based statistics and indices bibliometrics(C, Y) % adds year-weighted indices bibliometrics(C, Y, A) % adds author-weighted and combined indices bibliometrics(C, Y, A, 'Period', [2007 2009]) % restricts the analysis to a date range bibliometrics(C, Y, A, 'CareerEnd', 2009) % adds a descriptive pre/post career-phase comparison bibliometrics(C, Y, A, 'Plots', false) % suppresses every figure bibliometrics(C, Y, A, 'Interpret', false) % suppresses the interpretation-flags block

R = bibliometrics(C, Y, A); % also returns every computed quantity in R

The function prints the numerical results in the Command Window, opens the figures described above, and returns a structured output R (see the header of bibliometrics.m for the full field list, including R.indices, R.flags, R.profile, R.career, and R.authors).

🎛 Inputs C : Row vector of nonnegative numeric values. C(i) is the number of citations of the i-th paper. Mandatory, defines the size of the publication set.

Y : Optional row vector of publication years, same length as C. Enables years of activity, papers per year, and all age-weighted indices.

A : Optional row vector of author counts, same length as C. Enables author-weighted and per-author indices.

If Y or A are omitted or empty ([]), the corresponding groups of indices are skipped, but the rest of the analysis is still performed.

📤 Outputs Printed to the Command Window:

Descriptive statistics (total papers/citations, min/max/mode/median/mean, CV and adjusted CV, Gini coefficient, years of activity, papers per year, authors per paper and citations per author).
Citation indices: h (with a, m, delta-h), Fenner's chi, g (with delta-g), Jin's A-index, Kosmulski's h2, Zhang's e, R-index, i10-index, Sidiropoulos' normalized h.
Years-weighted indices (if Y): hc-index, Jin's and Harzing's AWCR/AR.
Authors-weighted indices (if A): Batista's hI, Harzing's hI,norm, Schreiber's hm.
Years and authors weighted indices (if Y and A): per-author AWCR/AR variants.
Extended indices (always, need only C): hg, w, maxprod, o, q2, tapered h (ht).
Interpretation flags (unless 'Interpret', false): sustainedness, concentration, breadth vs depth, recent activity.
Career-end descriptive profile (if 'CareerEnd' is given).

Returned in R: every quantity above, plus R.profile (the longitudinal series) and R.flags (level/symbol/value/message for each interpretation flag). See the function header for the complete field-by-field documentation.

🔍 Interpretation The indices implemented in bibliometrics cover multiple aspects of scientific impact:

Classical indices (h, g, h2, e, A) focus on citation counts and their distribution among the most cited papers.
The R-index summarizes the overall citation volume of the h-core as a single quantity, complementing the h-index and e-index.
The i10-index is a simple Google-style indicator that counts how many papers have accumulated at least 10 citations.
Normalized and age-weighted indices highlight the temporal evolution of impact, penalizing older work and emphasizing recent contributions.
Author-weighted indices attempt to correct for co-authorship by distributing credit among collaborators or by using fractional ranks.
Per-author age-weighted indices combine both aspects: each term in the sum refers to a specific paper characterized by its citations, age, and number of authors, all consistently aligned with the h-core ordering (for Jin's per-author AWCR) or with the full ranked list (for Harzing's per-author AWCR).
The extended indices add further nuance beyond h and g: the hg-index balances the two; the w-index and maxprod-index emphasize the very top of the citation distribution; the o-index folds in the single most-cited paper; the q2-index captures the typical (median) impact inside the h-core; the tapered h-index credits every paper proportionally instead of using a hard cutoff.
The Lorenz curve and Gini coefficient quantify inequality in the distribution of citations over the publication set.
The interpretation flags translate four of these signals (pace of output, citation concentration, breadth vs depth, and recent-paper share) into plain language, always naming the rule of thumb behind the threshold rather than presenting it as an absolute judgement.

Together, these measures provide a detailed, multi-dimensional view of a researcher's citation profile, and can help compare citation patterns beyond simple totals or a single index.

📝 Notes The function assumes that C, Y, and A refer to the same ordered list of papers, and that all data are correct and complete. No automatic disambiguation of names or merging of duplicate records is performed. For large publication lists, the plotting section may become compute-intensive, particularly the Ferrers diagram for the h-index, which draws a point for each citation; use 'Plots', false to skip all figures on very large datasets. The interpretation flags are heuristics for orientation, not pass/fail judgements, and each one states the rule of thumb it is drawing on. The internal version 2.4.1 adds the extended-index tier, the interpretation flags, and a third figure, and uses plain-text tags rather than Unicode emoji so that the report renders identically on every platform and font.

📚 Citation If you use this code in scientific, educational, or technical work, please cite it as:

Cardillo G. (2010) "Bibliometrics: the art of citation indices". Available from GitHub: https://github.com/dnafinder/bibliometrics

👤 Author Author: Giuseppe Cardillo Email: giuseppe.cardillo.75@gmail.com GitHub: https://github.com/dnafinder

⚖️ License This project is distributed under the MIT License. You are free to use, modify, and redistribute the code, provided that the original copyright notice and license text are preserved. The full license terms are provided in the LICENSE file in this GitHub repository.
This project is distributed under the MIT License. You are free to use, modify, and redistribute the code, provided that the original copyright notice and license text are preserved. The full license terms are provided in the LICENSE file in this GitHub repository.
