# bibliometrics
Compute pratically all bibliometric indices actually known.
Actually, every researcher come to terms with his/her bibliometric
performance. There are many software to compute these indices, like
Publish or Perish or Google Scholar. Anyway, these softwares require an
additional work like to identify the proper papers to exclude papers by
authors with the same name and surname. More over, Google Scholar is not
efficient to list all the authors of a papers (it cuts the list), so some
metrics are inaccurate.
So I decided to set-up this algorithm to compute these metrics without web
data-mining.

Syntax: 	bibliometrics(C,Y,A)
     
    Inputs:
          C - an array with the number of Citations of published papers
          (mandatory). These data can be obtained using several databases
          as ISI, Researcher ID, Google Scholar, Scopus, Citeseeker,...
          Y - an array with the Years of publication of the papers
          (optional). Everybody know when an own paper was published...
          A - an array with the number of Authors of published papers.
          Everybody know the number of co-authors of an own paper...
    Outputs (if all vectors were given):
          - Descriptive statistics:
              ° Total number of papers
              ° Total number of citations
              ° Min, Max, Mode, Median and Mean citations per papers (with 95% confidence intervals)
              ° Variation Coefficient (normal and adjusted)
		      ° Gini's Coefficient
              ° Years of activity, first and last year of publication
              ° Min, Max, Mode, Median and Mean papers per year (with 95% confidence intervals)
              ° Mean number of citations per year
              ° Min, Max, Mode, Median and Mean Authors per paper (with 95% confidence intervals)
              ° Citations per author
          - Bibliometrics indices
              ° Citations indices
                  * Hirsch's h-index with a and m parameters, delta-h
                  * Egghe's g-index, delta-g
                  * Jin's A-index
                  * Kosmulski's h2-index
                  * Zhang's e-index
                  * Sidiropoulos'es normalized h-Index 
              ° Years weighted indices
                  * Sidiropoulos'es Contemporary h-Index (hc-index)
                  * Jin's Age weighted citation rate (AWCR) and AR-index
                  * Harzing's Age weighted citation rate (AWCR) and AR-index
              ° Authors weighted indices
                  * Batista's Individual h-index (hI-index)
                  * Harzing's Individual h-index (hI,norm-index)
                  * Schreiber's Multi-authored h-index (hm-index)
              ° Years and Authors weighted indices
                  * Jin's Age weighted citation rate (AWCR) and AR-index
                  * Harzing's Age weighted citation rate (AWCR) and AR-index

              Several plots

  Example:
              C=[12 8 1 0 5 3 0 0];
              Y=[2004 2007 2008 2008 2008 2009 2009 2010];
              A=[8 9 10 7 11 11 7 5];
              bibliometrics(C,Y,A)

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2010) Bibliometrics: the art of citations indices
http://www.mathworks.com/matlabcentral/fileexchange/28161