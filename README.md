# hlsgr
Hiroshima Liver Study Group R Utilities

git clone https://github.com/hcnh174/hlsgr.git

devtools::install_github("hcnh174/hlsgr")

#install devtools
install.packages('devtools')
install.packages("roxygen2")
install.packages("testthat")

# install Rtools
https://cran.rstudio.com/bin/windows/Rtools/

.Renviron
PATH="${RTOOLS40_HOME}\usr\bin;${PATH}"

#install packages
install.packages('openxlsx')
install.packages('pheatmap')
install.packages('R.oo')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#BiocManager::install("KEGGprofile", update=FALSE)
BiocManager::install("biovizBase", update=FALSE)
BiocManager::install("ggbio", update=FALSE)
BiocManager::install("tximport", update=FALSE)
BiocManager::install("edgeR", update=FALSE)
BiocManager::install("DESeq2", update=FALSE)

#install.packages("WGCNA")
#BiocManager::install('GO.db', update=FALSE)
#BiocManager::install('org.Hs.eg.db', update=FALSE)
#BiocManager::install('impute', update=FALSE)
#BiocManager::install('preprocessCore', update=FALSE)

# links and tutorials
https://r4ds.had.co.nz/tibbles.html
http://r-pkgs.had.co.nz/tests.html
http://r-pkgs.had.co.nz/git.html
https://rawgit.com/rstudio/cheatsheets/master/package-development.pdf
https://tinyheero.github.io/jekyll/update/2015/07/26/making-your-first-R-package.html
https://stackoverflow.com/questions/45044269/how-to-use-data-within-a-function-in-an-r-package
https://mastering-shiny.org/
https://community.rstudio.com/t/magrittr-inside-a-package/2033/9
https://community.rstudio.com/t/import-from-magrittr-package/4313
https://stackoverflow.com/questions/27947344/r-use-magrittr-pipe-operator-in-self-written-package
https://towardsdatascience.com/creating-r-packages-what-you-need-to-know-2a20233b328a
https://travis-ci.org/
https://daattali.gitbooks.io/stat545-ubc-github-io/bit003_api-key-env-var.html

