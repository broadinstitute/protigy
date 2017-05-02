pck <- c('shiny', 'scales', 'gtable', 'limma', 'RColorBrewer',
         'hexbin', 'Hmisc', 'grid', 'scatterplot3d', 'plotly',
         'WriteXLS', 'reshape', 'nlme', 'BlandAltmanLeh',
         'preprocessCore', 'mice', 'mixtools', 'mclust',
         'shinydashboard', 'ChemometricsWithR', 'devtools',
         'maptools', 'DT', 'ggrepel', 'dplyr')

setRepositories(graphics=F)
install.packages(pck, dependencies=T)


pck.bioc <- c("org.Hs.eg.db")

source("https://bioconductor.org/biocLite.R")
biocLite(pck.bioc)

