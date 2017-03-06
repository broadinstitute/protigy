pck <- c('shiny', 'scales', 'gtable', 'limma', 'RColorBrewer',
         'hexbin', 'Hmisc', 'grid', 'scatterplot3d', 'plotly',
         'WriteXLS', 'reshape', 'nlme', 'BlandAltmanLeh',
         'preprocessCore', 'mice', 'mixtools', 'mclust',
				 'shinydashboard', 'ChemometricsWithR')

setRepositories(graphics=F)

install.packages(pck, dependencies=T)
