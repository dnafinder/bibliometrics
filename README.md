[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/bibliometrics)

📌 Overview
This repository provides the MATLAB function bibliometrics, which computes a broad range of bibliometric indicators for a single researcher. Given the citation counts of a set of publications (and optionally the publication years and number of authors per paper), the function reports descriptive statistics, many well-known citation indices, and several diagnostic plots. The goal is to reproduce and extend the type of analysis provided by tools such as Publish or Perish or Google Scholar, but without relying on online data mining and without the ambiguities of author name matching.

✨ Features
The function bibliometrics:
- Computes descriptive statistics for citations, years, and authors (when provided).
- Calculates a large number of bibliometric indices, including h, g, h2, e, A, R-index, i10-index, contemporary h, age-weighted citation rates, individual h-index variants, and multi-authored indices.
- Provides inequality analysis via the Lorenz curve and Gini coefficient.
- Generates several plots that visually illustrate the indices, such as the Lorenz curve, g-index and h2-index plots, age-weighted citation curves, and the Ferrers diagram interpretation of the h-index.
All computations are done locally on the data you supply, so there is no dependence on external web services or online databases.

🛠 Installation
Download or clone this repository from GitHub:
https://github.com/dnafinder/bibliometrics

Add the folder containing bibliometrics.m to your MATLAB path using the Add Folder to Path option or the addpath command. The function uses only core MATLAB functionality (graphics, basic statistics, and datetime) and does not require additional toolboxes.

▶️ Usage
At minimum, you must provide a row vector of citation counts for each of your publications. Optionally, you may also provide publication years and number of authors per paper in the same order.

Examples:
C = [12 8 1 0 5 3 0 0];
Y = [2004 2007 2008 2008 2008 2009 2009 2010];
A = [8 9 10 7 11 11 7 5];

bibliometrics(C)           % only citation-based statistics and indices
bibliometrics(C, Y)        % adds year-weighted indices
bibliometrics(C, Y, A)     % adds author-weighted and combined indices

The function prints the numerical results in the Command Window and opens figures with several plots. It does not return output variables.

🎛 Inputs
The function accepts up to three positional inputs:

C : Row vector of nonnegative numeric values. C(i) is the number of citations of the i-th paper. This input is mandatory and defines the size of the publication set.

Y : Optional row vector of publication years. Y(i) is the year in which the i-th paper was published. All entries must be positive integers. If provided, Y must have the same length as C. Year information is used to compute years of activity, papers per year, and all age-weighted indices.

A : Optional row vector of author counts. A(i) is the total number of authors on the i-th paper. All entries must be positive integers. If provided, A must have the same length as C. Author information is used to compute author-weighted and per-author indices.

If Y or A are omitted or empty ([]), the corresponding groups of indices (years-weighted or authors-weighted) are skipped, but the rest of the analysis is still performed.

📤 Outputs
The bibliometrics function does not return any output arguments. Instead, it:

- Prints descriptive statistics such as:
  - Total number of papers
  - Total number of citations
  - Minimum, maximum, mode, median, mean citations per paper
  - Coefficient of variation (and adjusted coefficient of variation)
  - Gini coefficient for citation distribution
  - Years of activity and papers per year (if Y is provided)
  - Authors per paper and citations per author (if A is provided)

- Prints bibliometric indices such as:
  - Hirsch's h-index, the associated proportionality constant a, and delta-h
  - Fenner's chi-index
  - Egghe's g-index and delta-g
  - Jin's A-index
  - Kosmulski's h2-index
  - Zhang's e-index
  - R-index, defined as the square root of the total citations in the h-core
  - i10-index, defined as the number of papers with at least 10 citations
  - Sidiropoulos' normalized h-index
  - Sidiropoulos' contemporary h-index (hc-index) and its a parameter (if Y)
  - Jin's and Harzing's age-weighted citation rates (AWCR) and AR-indices (if Y)
  - Batista's individual h-index, Harzing's hI,norm, and Schreiber's hm-index (if A)
  - Age-weighted and per-author variations of Jin's and Harzing's indices (if Y and A)

🔍 Interpretation
The indices implemented in bibliometrics cover multiple aspects of scientific impact:
- Classical indices (h, g, h2, e, A) focus on citation counts and their distribution among the most cited papers.
- The R-index summarizes the overall citation volume of the h-core as a single quantity, complementing the h-index and e-index.
- The i10-index is a simple Google-style indicator that counts how many papers have accumulated at least 10 citations.
- Normalized and age-weighted indices highlight the temporal evolution of impact, penalizing older work and emphasizing recent contributions.
- Author-weighted indices attempt to correct for co-authorship by distributing credit among collaborators or by using fractional ranks.
- Per-author age-weighted indices combine both aspects: each term in the sum refers to a specific paper characterized by its citations, age, and number of authors, all consistently aligned with the h-core ordering (for Jin's per-author AWCR) or with the full ranked list (for Harzing's per-author AWCR).
- The Lorenz curve and Gini coefficient quantify inequality in the distribution of citations over the publication set.

Together, these measures provide a detailed, multi-dimensional view of a researcher's citation profile, and can help compare citation patterns beyond simple totals or a single index.

📝 Notes
The function assumes that C, Y, and A refer to the same ordered list of papers, and that all data are correct and complete. No automatic disambiguation of names or merging of duplicate records is performed. For large publication lists, the plotting section may become compute-intensive, particularly the Ferrers diagram for the h-index, which draws a point for each citation. If necessary, you can comment out or adapt the plotting section for very large datasets. The internal version 2.2.0 adds the R-index and i10-index while keeping Jin's per-author AWCR aligned with the h-core ordering.

📚 Citation
If you use this code in scientific, educational, or technical work, please cite it as:

Cardillo G. (2010)
"Bibliometrics: the art of citation indices".
Available from GitHub:
https://github.com/dnafinder/bibliometrics

👤 Author
Author: Giuseppe Cardillo
Email: giuseppe.cardillo.75@gmail.com
GitHub: https://github.com/dnafinder

⚖️ License
This project is distributed under the MIT License. You are free to use, modify, and redistribute the code, provided that the original copyright notice and license text are preserved. The full license terms are provided in the LICENSE file in this GitHub repository.
