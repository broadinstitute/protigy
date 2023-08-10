# Last updated May 4, 2022 by Natalie Clark (nclark@broadinstitute.org) - v 1.0.1

############################
## UI part
printHTMLUI <- function(id) {
    ns <- NS(id)

    htmlOutput(ns("html"))
}

#############################
## server part
printHTML <- function(input, output, session, what, error=NULL, global.input=NULL, global.param=NULL){

    txt=''

    ## ##@#############################
    ## ## getting started
    if(what == 'gs'){
      
      ## render HTML
      output$html <- renderText({
        if(!is.null(global.input$file)) return()
        if(!is.null(error$msg)) return()
        includeMarkdown('README.md')
        
      })
    }
    
    ##@###############################
    ## change log
    if(what == 'cl'){
      txt <- '<h4><font color="red">What\'s new?</font></h4>
<font size=\"3\">
<b>v1.1.7 August 10, 2023</b>\
<ul>
<li>Fixed an issue where statistical analysis results would not display properly for files with only one comparison after a missing value filter was applied.
</ul>
<b>v1.1.6 June 27, 2023</b>\
<ul>
<li>Fixed an issue where, with certain input files, the raw average expression/raw fold change columns were not sorted properly in the output data files.
</ul>
<b>v1.1.5 May 17, 2023</b>\
<ul>
<li>Fixed an issue where columns were erroneously named if one column name was a substring of another column name. This fix required re-naming some of the raw expression value columns.
<li> Replaced the zero-centered model coefficients in the moderated F test table with the normalized average expression values.
<li> Removed pairwise information from moderated F test table. 
<li> Generalized profile plot axis labels.
</ul>
<b>v1.1.4 March 7, 2023</b>\
<ul>
<li>Fixed an issue with the one-sample and two-sample T-test. Previously, in the rare case where a protein was not detected in all samples in a certain group, any test comparisons involving that group would be erroneously discarded, which would cause an error when attempting to export the signed log-transformed p-values. Now, all test results are reported for all proteins: NAs are reported when a protein was excluded from a particular test. 
<li> Fixed an issue where exporting the RMarkdown report would fail when attempting to export the heatmap.
<li> Fixed an issue where a GCT file with only one cdesc column could not be imported correctly.
<li> Fixed some column descriptors in output GCT files.
<li> Updated various help text within the app.
</ul>
<b>v1.1.3 February 8, 2023</b>\
<ul>
<li>Occasionally, eBayes(trend=TRUE) fails for intensity-based data, particularly when the distribution of quantified features is not uniform across samples. In these cases, eBayes(trend=FALSE) is run instead, and a warning message is printed. We highly encourage users who encounter this warning to carefully examine their data, and re-perform statistical analysis as needed. Typically, setting a stricter missing value filter will fix the issue.
<li> The one-sample T-test is now fixed. It was not working due to an error in passing a parameter to the function.
<li> An issue with heatmap color annotations was fixed. The issue occured during conversion of certain annotation group names to R-compliant names.
</ul>
<b>v1.1.2 January 18, 2023</b>\
<ul>
<li>Previously, raw log-fold changes were median-centered by default. Raw log-fold changes are no longer median-centered so that they accurately reflect the fold change before any normalization is performed.
<li>When a one-sample T-test or no test is performed, raw average abundances per group are reported in the place of fold change. Column headers and descriptors in the exported excel file have been updated to reflect this for these cases.
<li>Help text has been updated to clarify that data must be log-transformed prior to statistical analysis. The statistical testing performed in Protigy is not valid for non-log-transformed data. 
</ul>
<b>v1.1.1 January 5, 2023</b>\
<ul>
<li>A checkbox has been added for group-wise normalization on the group assignment screen. When this box is checked, group normalization will be performed based on the selected column. When the box is left unchecked, group normalization will not be performed. This replaces the previous checkbox on the statistical analysis screen.
<li>There is now a help button that explains how the missing value filter works, as well as help text within the application.
<li>The correlation boxplot is now always scaled from 0 to 1 for positive values. For negative values, the boxplot will scale based on the minimum value. This should produce more comparable boxplots across datasets.
<li>More information has been added to the workflow summary including if intensity data were used, if group-wise normalization was performed, and if QC.fail samples were filtered.
</ul>
<b>v1.1 November 28, 2022</b>\
<ul>
<li>A checkbox has been added to filter out samples where QC.status=QC.fail. This option is for GCT files ONLY. If the QC.status column does not exist and the box is checked, a warning message is printed and analysis will proceed using all samples.
<li>A warning message appears when attempting to log-transform a dataset that contains negative values (i.e., a dataset that has already been log-transformed). In this case, log-transformation will not occur to prevent downstream analysis errors.
<li>An error message appears when a column label only has one sample assigned to it, and analysis will not proceed. ProTIGY requires all column labels have more than one sample assigned to them. When this error appears, the user must either remove that sample from the file, or use another column for statistical testing and/or group-wise normalization.
<li>A checkbox has been added to indicate the data are intensity data rather than ratios. When the box is checked, normalization methods, filtering methods, and statistical testing will be automatically filtered to only methods appropriate for intensity data. Further, the missing value rate will be capped at 99%, as the statistical analysis can fail when the missing value rate is 100%.
<li> For the moderated two-sample t-test and the moderated F test, eBayes(trend=true) is used for intensity data. This option is more stable for intensity data, but does require the capped missing value rate of 99%. eBayes(trend=FALSE) continues to be used for ratios.
<li> Help text has been edited for non-GCT files to clarify that NA must be used for missing sample annotations. Further, a sample can be excluded in a non-GCT file by leaving the Experiment and Group annotations as NA.
</ul>
<b>v1.0.4 September 19, 2022</b>\
<ul>
<li>Fixes multiple plotting issues.
<li>Volcano plot labels may now be chosen from the following: ID_Symbol (default), ID, or Symbol.
<li>Gene symbol support for Ensembl protein IDs is now included.
<li>For gct files, if a geneSymbol column is included in the row descriptors (rdesc), it is used to determine the symbol rather than using the available database.
</ul>
<b>v1.0.3 July 21, 2022</b>\
<ul>
<li>Fixes issues with heatmap visualization and export.
</ul>
<b>v1.0.2 May 10, 2022</b>\
<ul>
<li>Fixes missing sample annotation columns in output .gct files.
</ul>
<b>v1.0.1 May 4, 2022</b>\
<ul>
<li>Blanks ("") are now read in as blanks rather than missing values. This is important to retain sample annotation information from .gct files.
<li>When not performing statistics (statistical test set to "none"), the expression values in the .gct file are no longer repeated twice.
<li>The export template file is now fixed to contain NA (missing values) rather than blanks.
<li>Multiple characters such as "na" (and all capitalization variations of NA) are classified as missing values.
</ul>
<b>v1.0.0 May 2, 2022</b>
<ul>
<li>Heatmap export is now fixed. 
<li>App no longer crashes when attempting to export all the files without performing statistics (statistical test set to "none").
<li>For non-gct files, the experimental design file has been streamlined. Users may still download the template experimental design file for use with Protigy. Alternatively, a experimental design file with any number of sample annotations and any type of text delimter (tab-delimited, comma-separated, etc) may be uploaded by the user as long as the first column of the design file contains the sample (column) names. These must match the column names of the table exactly.
<li>Non-gct file users may now select the columns they wish to use for statistical analysis and/or group-wise normalization using drop-down menus as they would for a .gct file.
<li>Help text now says to refresh the page to start a new session. F5 is a Windows-only shortcut and not applicable for all users.
<li>Blanks ("") are now read in as missing values (NA) rather than characters.
</ul>
<b>v0.9.1.5 Apr 13, 2022</b>
<ul>
<li>Now accepts a third annotation column to denote group-wise normalization separately from groups for statistical testing. 
<li>Fixes an issue with missing sample annotations. Previously, samples with missing annotations were sometimes included in the analysis. Now, samples with missing annotations are not included in any part of the analysis.
<li>Other quality of life changes, such as automatically adjusting plot axes to specify the type of p-value used (non-adjusted vs adjusted)
<li>Improved help text in various regions of the application.
</ul>
<b>v0.9.1.4 Oct 21, 2021</b>
<ul>
<li>mod F test: use row centered data that makes the F test more interpretable (to identify groups that are different, as opposed to groups with non-zero average values)
<li><code>logFC.raw</code></li>-columns are now median-centered
</ul>
<b>v0.9.1.3 Aug 30, 2021</b>
<ul>
<li>Result-tables: Fixed a bug in reporting <code>logFC.raw</code> that triggered a crash when the missing data filter was applied.</li>
<li>Misc: Updated parsing of user database for shared sessions (RSC/SSP only) to be compatible with the user authentication that RSC at Broad is using. Shared sessions are now visible again.</li>
<li>Misc: Fixed a few typos.</li>
</ul>
<b>v0.9.1.2 Aug 26, 2021</b>
<ul>
<li>Result tables: Added column for the log fold change before normalization (<code>logFC.raw</code>).</li>
<li>Table-tab: Set number of decimals to 3.</li>
<li>Excel-export: Made the error message more explicit.</li>
<li>Misc: Some edits to the README file.</li>
</ul>
<b>v0.9.1.1 July 30, 2021</b>
<ul>
<li>Profile plots: Switched to interactive plots using <code>plotly</code>. Profile plots in the zip-file remain static.</li>
<li>Bugfix: Distributions of unnormalized data in the .pdf verions of the box and profile plots when downloaded as .zip file are no shown correctly (pdf version would show normlaized distributions). This bug was probably introduced with the missing value filter in version v0.8.9.6.</li>
</ul>
<b>v0.9.1 June 26, 2021</b>
<ul>
<li>Misc: Fixed a bug that prevented the upload of the experimental design file, if an alternate id column was chosen (e.g. <code>accession_number</code> or <code>geneSymbol</code> in Spectrum Mill reports) AND a column <code>id</code> was already present in the uploaded text file.</li>
</ul>
<b>v0.9.0 May 27, 2021</b>
<ul>
<li>PC plots: Fixed colors if annotation data contains NA.</li>
<li>SD filter: Fixed crash triggered by selecting SD filter after the "run-test" button has been pressed.</li>
<li>Missing value filter: Replace slider input by numeric input widget.</li>
<li>Misc: Updated URL to new conf-app (RStudio Connect only).</li>
<li>Misc: Reorganized R-package managment for deployment to RStudio Connect. See parameter <code>PACMAN</code> in <code>global.R</code>.</li>
</ul>
<b>v0.8.9.7 May 11, 2021</b>
<ul>
<li>Normalization: Added checkbox for group-level normalization. If enabled the normalization will be performed within a particular group (Median, Median-MAD, Quantile, VSN). For Median and Median-MAD normalization, the group-level median of sample medians is added to each normaized data value.</li>
</ul>
<b>v0.8.9.6 April 30, 2021</b>
<ul>
<li>Filter: Added filter for missing values.</li>
<li>Scatterplots: Added correlation coefficients.</li>
<li>Scatterplots: Fixed legend for filtered values if StdDev was selected.</li>
<li>Misc: User selections when setting up the analysis worflow are now remembered when changing selecting a different value for <i>Filter data</i></li>
<li>Misc: Removed dependency with deprecated package <code>prada</code>.</li>
</ul>
<b>v0.8.9.5 April 29, 2021</b>
<ul>
<li>PCA: Added the possibility to use different annotations to color the PC plots. Only available for GCT v1.3 input files.</li>
</ul>
<b>v0.8.9.4 April 28, 2021</b>
<ul>
<li>Normalization: Added normalization methods for intensity data (e.g. label-free quantification in proteomics): <b>Median (log-intensity)</b>, <b>Median-MAD (log-intensity)</b>, <b>VSN (intensity)</b></li>
<li>Profile plots: Users can now control the scale of the x-axis (<code>symmetric</code> or <code>as-is</code>)</li>
<li>Misc: Changed color scheme to <i>pink</i> in preparation to transition from SSP to RSC.</li>
</ul>
<b>v0.8.9.3 March 26, 2021</b>
<ul>
<li>Misc: Switched to GitHub version of package <code>ChemometricsWithR</code> as it got removed from CRAN.</li>
<li>Misc: Robustified category names of the class vector if selected from a GCT 1.3 file.</li>
</ul>
<b>v0.8.9.2 February 11, 2021</b>
<ul>
<li>Misc: SSP/RSC only - path to configuration file for sharing sessions between users can be specified in <code>global.R</code>.</li>
</ul>
<b>v0.8.9.1 February 05, 2021</b>
<ul>
<li>Misc: Fixed the error message "Experimental design" file does not match the table you have uploaded (different number of rows/columns)!" that was falsely triggered whenever the selected id column was anything other than <code>id</code>. The bug was intrduced in v0.8.9</li>
<li>Misc: Fixed a typo in the column description of the Excel result sheet ("Nomical P-value" -> "Nominal P-Value")
</ul>

<b>v0.8.8.1 December 3, 2020</b>
<ul>
<li>PCA: Added the number of total features to the title of the plot.</li>
</ul>

<b>v0.8.8 October 23, 2020</b>
<ul>
<li>GCT export: Fixed bug that caused the GCT export to fail if the column name for ids was not <code>id</code>.</li>
<li>Normalization: Improved error handling if 2-component normalization fails to converge.</li>
<li>Normalization: Added upper quartile normalization (subtract 75th percentile).</li>
<li>Misc: Some code cleanup.</li>
</ul>

<b>v0.8.7 August 05, 2020</b>
<ul>
<li>GCT v1.3 import: Error message when name of group variable contains special characters.</li>
<li>GCT v1.3 import: Robustified column names of cdesc object.</li>
<li>Misc: fixed inconsistent capitalization of file extentions (<code>.r</code> vs. <code>.R</code>)</li>
<li>Misc: data folder set to <code>tempdir()</code> when running locally</li>
</ul>

<b>v0.8.6.3 April 01, 2020</b>
<ul>
<li>Heatmap: Disabled interactive heatmap due to incompatibility with newer verions of the heatmaply package.</li>
<li>Excel Sheet: Fixed broken column desciptions that occured when testing for prefix/suffix adn the column the respective column had been annotated already.</li>
<li>Scatterplots: Columns from other experiments can now be selected.</li>
<li>Misc: Improved error handling and error messages when uploading the experimental desing file.</li>
</ul>
<b>v0.8.6.2 March 27, 2020</b>
<ul>
<li>Misc: Improved handling of redundant ids.</li>
</ul>
<b>v0.8.6.1 March 18, 2020</b>
<ul>
<li>Heatmap: Fixed number of significant features.</li>
<li>Heatmap: Interactive heatmap working again.</li>
<li>Correlation boxplots: Included in exports (.zip) and Rmarkdown reports.</li>
</ul>
<b>v0.8.6 March 9, 2020</b>
<ul>
<li>Excel Sheet: Added column descriptions for Protigy and Spectrum Mill-specific columns.</li>
<li>Experimental design: Robustified handling of special characters in experiment names.</li>
<li>Volcano plots: Added box around legend for highlighted proteins/PTM-sites.</li>
<li>PCA plots: Number of features used for PCA shown in  title.</li>
<li>PCA plots: Central function in <code>src/plots.r</code></li>
</ul>
<b>v0.8.5.5 September 3, 2019</b>
<ul>
<li>Summary tab: Fixed the bug that would cause an error in the \'Dataset:\'-widget.</li>
</ul>
<b>v0.8.5.4 August 13, 2019</b>
<ul>
<li>Excel sheet/GCT export: Fixed a bug that that occasionally would mess up the header in the Excel/GCT file. Only happened in <b>Two-sample test</b> and if very similar names for experiments were used (e.g. one experiment name is a substring of another experiment name).</li>
<li>Summary tab: Number of features w/o any quant is now correctly reported in previously saved sessions.</li>
</ul>
<b>v0.8.5.3 June 11, 2019</b>
<ul>
<li>Heatmap: Accession and gene name are shown as row annotation.</li>
<li>SD filter: Number of remaining features is reported correctly on the summary page.</li>
</ul>
<b>v0.8.5.2 June 10, 2019</b>
<ul>
<li>Misc: Added <i>toggle all</i> buttons in group selection dialog.</li>
<li>Misc: Draggable modal windows using the <code>shinyjqui</code> R package.</li>
<li>Misc: Improved error reporting using the <code>shinyalert</code> R package.</li>
</ul>
<b>v0.8.5.1 June 6, 2019</b>
<ul>
<li>GCT: Added GCT file with singed, log-transformed p-values to output zip-archive which can be used as ranking input for ssGSEA/PTM-SEA.</li>
<li>PPI: Fixed bug that would show protein-protein interactions of the first protein in the list as default in volcano and scatterplots in R vesion >=3.5.</li>
</ul>
<b>v0.8.5 Jan 21, 2019</b>
<ul>
<li>Misc: added BSD-3 license.</li>
<li>Misc: updated Readme file.</li>
</ul>
<b>v0.8.4.3 Jan 17, 2019</b>
<ul>
<li>Misc (SSP only): Started to implement a shiny module for session management which will enable users to share saved sessions with collaborators/team members. Not used in this version though.</li>
</ul>
<b>v0.8.4.2 Jan 10, 2019</b>
<ul>
<li>Misc: compatible with both, R>=3.5 AND R<3.5</li>
</ul>
<b>v0.8.4.1 Dec 5, 2018</b>
<ul>
<li>Misc: compatible with R >=3.5</li>
</ul>
<b>v0.8.4 Nov 16, 2018</b>
<ul>
<li>Multiscatter: Fixed bug that would show straight lines for each pairwise plot. Occured when column ids where longer than 20 characters.</li>
<li>GCT 1.3: Error message if GCT file does not contain any column meta data tracks.</li>
<li>Misc: Fixed a bug causing the app to crash if a GCT 1.3 file with single <b>column meta data track</b> was uploaded.</li>
<li>Misc: Renamed  <i>Modify selected groups</i> to  <i>Select groups</i>.</li>
</ul>
<b>v0.8.3.1 July 24, 2018</b>
<ul>
<li>Misc: Robustified filtering of significant features for plotting purposes.</li>
<li>Misc: Disabled Javascript code in the <code>datatable()</code>.</li>
<li>Correlation boxplots: Changed some more aesthetics.</li>
</ul>
<b>v0.8.3 July 23, 2018</b>
<ul>
<li><mark>BUG:</mark> fixed a bug resulting in an inaccurate number of significant features reported in the heatmap.</li>
<li>UpSet-plot: small bugfix causing a crash under certein circumstances.</li>
<li>Correlation boxplots: changed some aesthetics of the plot.</li>
</ul>
<b>v0.8.2.8 July 2, 2018</b>
<ul>
<li>Misc: fixed a bug causing the app to crash if a GCT 1.3 file <b>without row meta data</b> was uploaded.</li>
</ul>
<b>v0.8.2.7 June 28, 2018</b>
<ul>
<li>Misc: fixed a bug causing the app to crash under certain combinations of <i>Modify selected groups</i> and test selections.</li>
</ul>
<b>v0.8.2.6 June 28, 2018</b>
<ul>
<li>Misc: disabled the cmapR-package because of installation problems of the required package <code>rhdf5</code> on a <i>Red Hat Enterprise Linux 6.9</i> machine. The io.R file from the cmapR GitHub repository is used instead. </li>
</ul>
<b>v0.8.2.5 June 28, 2018</b>
<ul>
<li>Excel sheet: in case of <b>Two sample moderated T-test</b>, the table header will now report <b>KO.over.WT</b> instead of <b>WT.vs.KO</b>.</li>
<li>Export: result files will also be epxorted in GCT v1.3 format.</li>
<li>Barplot: fixed a bug causing the barplot of identified features mislabel the colors if the <i>Modify selected groups</i>-feature was used.</li>
<li>Boxplots: fixed a bug resulting in slightly different numbers reported in the exported pdf file compared to the numbers shown in the app. This only happened in boxplots depicting values after normalization.</li>
</ul>
<b>v0.8.2.4 June 27, 2018</b>
<ul>
<li><mark>BUG:</mark> duplicated session ids: included a timestamp-based seed before generating the seesion id. It happened that the sample function returned the same string/number combination.</li>
<li>Misc: session id is doubled checked whether it exists as folder on the server.</li>
<li>Export: all files except RData-session files and zip-archives are removed from the server.</li>
</ul>
<b>v0.8.2.3 April 20, 2018</b>
<ul>
<li>Misc: updated code for 2-component normalization (by D. R. Mani).</li>
</ul>
<b>v0.8.2.2 April 17, 2018</b>
<ul>
<li>Heatmap: row and column labels can be disabled now.</li>
<li>UpSet plots: inter-group comparison of significantly regulated features.</li>
</ul>
<b>v0.8.2.1 March 14, 2018</b>
<ul>
<li>Correlation: correlation matrix calculated centrally in a separate function and shared with plots using correlations: multiscatter, correlation heatmap, correlation boxplot.</li>
<li>Correlation: novel QC-tab depicting pairwise intra-group correlations as boxplots.</li>
</ul>
<b>v0.8.2 February 27, 2018</b>
<ul>
<li>Misc: Installable on Mac OS.</li>
<li>Misc: group selection now correctly updated in saved sessions.</li>
</ul>
<b>v0.8.1 February 26, 2018</b>
<ul>
<li>Misc: simplified installation under Windows OS.</li>
<li>Misc: if Perl and/or Pandoc are not availbale the app will show a corresponding messsage.</li>
<li>PCA: "Run me first"-tab became obsolete.</li>
<li>Scatterplots: added trace of filtered values for reprodicibility filter.</li>
<li>Scatterplots: separated data tracks and added legend.</li>
</ul>
<b>v0.8.0.9 February 24, 2018</b>
<ul>
<li>GCT 1.3: robustified import of GCT 1.3 files. If not unique, row and column identifiers are made unique.</li>
<li>Volcano: Labeled points can be removed individually from the table.</li>
<li>Table-tab: page overhaul</li>
</ul>
<b>v0.8.0.7 February 22, 2018</b>
<ul>
<li>Misc: PPI queries now work after export of results.</li>
</ul>
<b>v0.8.0.6 February 21, 2018</b>
<ul>
<li>UI: updated help text</li>
<li>UI: Only first 20 characters of column names are shown when prompted to select ID column.</li>
</ul>
<b>v0.8.0.5 February 20, 2018</b>
<ul>
<li>Export: fixed a bug preventing the export of results as zip-archive.</li>
</ul>
<b>v0.8.0.4 February 15, 2018</b>
<ul>
<li>Misc: Robustified import of gct 1.3 files (row and column names are made unique).</li>
</ul>
<b>v0.8.0.3 February 14, 2018</b>
<ul>
<li>Misc: new session import/export features</li>
</ul>
<b>v0.8.0.2 February 14, 2018</b>
<ul>
<li>Export: page overhaul</li>
<li>Export: generation of Rmarkdown-reports (still under developement).</li>
<li>Export: added option to download rmarkdown, xls, zip, separately.</li>
<li>Volcanos: color overhaul.</li>
<li>PPI: fixed a bug in which multiple occurences of selected bait proteins were not shown in zoomed view.</li>
<li>PPI: fixed the <i>all-turns-green</i> bug.</li>
<li>PPI: added ID mapping support for <i>mus musculus</i>, <i>rattus norvegicus</i> and <i>danio rerio</i>.</li>
<li>Multiscatter: robust determination of plotting limits.</li>
<li>Multiscatter: re-drawing only after button was pressed.</li>
</ul>
<b>v0.8.0.1 February 06, 2018</b>
<ul>
<li>Misc: Piwik integration.</li>
<li>Fanplot: colors are synchronized with current group selection.</li>
<li>Fanplot: added legend and possibility to modify labels.</li>
</ul>
<b>v0.8.0 January 25, 2018</b>
<ul>
<li>Release version for SSP (dev),</li>
<li>Session import: improved backwards compatibility.</li>
<li>Export: data directory is cleaned up now. Only .RData session files and the latest zip archive remain in the user/session data directory.</li>
</ul>
<b>v0.7.8.4 January 24, 2018</b>
<ul>
<li>Heatmap: interactive heatmap using "heatmaply".</li>
<li>Heatmap: annotation tracks (GCT 1.3) can be selected/deselected.</li>
<li>Clustering: default distance metric switched from <b>euclidean</b> to <b>1-Pearson</b>.</li>
<li>Clustering: Fanplot v0.1 - circular dendrogram to visualize sample clustering.</li>
<li>PCA: added legend to plots.</li>
<li>Misc: links to "Genecards" if ids are not UniProt.</li>
<li>Multiscatter: BA-filtered values shown in blue.</li>
</ul>
<b>v0.7.8.3 January 22, 2018</b>
<ul>
<li>PPI scatterplots: reduced opacity for non-interactors.</li>
<li>PPI: robustified extraction of gene symbols in function"get.interactors()".</li>
<li>Summary: new plot for missing values.</li>
<li>SSP import: backwards compatibility, sessions saved from older versions can be imported.</li>
<li>Gene name mapping: fixed bug that would cause a crash if neither UniProt nor RefSeq ids were found.</li>
<li>Gene name mapping: gene names that could not mapped are indicated by "NotFound".</li>
<li>Misc: improved start up time of the app using function "import.ppi.db()"</li>
<li>Misc: working button in the "Select Groups" modal window.</li>
</ul>
<b>v0.7.8.2 December 29, 2017</b>
<ul>
<li>Volcano: fixed overlaping legends.</li>
<li>Volcano: fixed fdr line bug.</li>
<li>Volcano: IDs are site-specific. Also effects PPI panel, i.e. a query always returns a single site rather than all sites mapping to a gene symbol.</li>
<li>Heatmap: GCT v1.3 annotation columns shown as tracks.</li>
<li>Misc: Groups defined in the experimental design or in GCT v1.3 annotation tracks can be enabled/disabled for testing.</li>
<li>Gene name mapping: finally works with RefSeq ids.</li>
</ul>
<b>v0.7.8.1 December 25, 2017</b>
<ul>
<li>Misc: added support for GCT v1.3 files. Class vector can be selected from column meta data.</li>
</ul>
<b>v0.7.8 December 4, 2017</b>
<ul>
<li>Heatmap: had to disable Morpheus widget since it would interfere with interactivity of volcono plots.</li>
<li>Misc: switched to "selectizeInput" to select saved sessions.</li>
<li>Misc: re-organization of navbarPage creation to fix an error thrown after Shiny R-packge update (v1.0.5)</li>
<li>Misc: integrated Readme.html into entry page.</li>
</ul>
<b>v0.7.7 October 30, 2017</b>
<ul>
<li>Heatmap: Morpheus integration (ALPHA)</li>
<li>Gene name mapping: robustified mapping if no RefSeq or UniProt ids were used.</li>
</ul>
<b>v0.7.6 September 02, 2017</b>
<ul>
<li>Misc: switched to <i>pacman</i> R package managment system.</li>
<li>Misc: added Readme on GitHub.</li>
<li>Normalization: Turned off automatic centering of Quantile-normalized data.</li>
</ul>
<b>v0.7.5 August 18, 2017</b>
<ul>
<li>Volcano: re-organization of PPI legends.</li>
<li>Scatterplots: PPI analysis is now fully integrated.</li>
<li>Heatmap: row annotations are shown correctly again.</li>
<li>Misc: gene mapping doesn\'t crash if no test was selected.</li>
<li>Misc: fixed a couple of other smaller bugs, mostly related to data exploration without performing a test.</li>
</ul>
<b>v0.7.4 August 10, 2017</b>
<ul>
<li>Scatterplots: new tab that provides interactive scatterplots between replicate measurements. For One-sample moderated T-test and F-test the significant proteins are marked in red.</li>
<li>Volcano: PPI - search mask keeps working after multiple rounds of analysis.</li>
<li>Misc: parameter file: fixed NA for data filter.</li>
</ul>
<b>v0.7.3.1 August 4, 2017</b>
<ul>
<li>Volcano plots: Fixed a bug causing volcano plots to crash when points were selected, but no protein-protein interactors were found.</li>
</ul>
<b>v0.7.3 June 20, 2017</b>
<ul>
<li>Misc: unified the naming of the id-column throughout the code.</li>
<li>Volcano plots: Integration of Reactome (human) protein-protein interactions.</li>
</ul>
<b>v0.7.2 June 6, 2017</b>
<ul>
<li>Volcano plots: added hyperbolic curves based on a minimal fold change cut-off and adjusted p-values.</li>
</ul>
<b>v0.7.1 June 1, 2017</b>
<ul>
<li>Misc: fixed a bug preventing the filter to be triggered after the user re-runs an analysis.</li>
<li>Volcano plots: integration of BioGRID (human) protein-protein interactions.</li>
<li>Volcano plots: <i>selectizeInput</i> now rendered on the server. Significantly speeded up page respond times.</li>
<li>Correlation matrix: updated color scheme to better visualize subtle differences.</li>
<li>Gene name mapping: fixed some bugs causing the app to crash if no gene names could be mapped or if other accessions than UniProt or RefSeq were used..</li>
<li>Misc: loading animation.</li>
</ul>
<b>v0.7.0 May 5, 2017</b>
<ul>
<li>Misc: automatic mapping to gene names if RefSeq or UniProt accession numbers were found in "id" column.</li>
<li>Volcano plots: integration of InWeb database.</li>
<li>Volcano plots: paramater "max. Log10(p-value)" works in all volcano plots. Before, changing the parameter only worked in the parameter panel of the first volcano.</li>
<li>Volcano plots: completely zoomable</li>
<li>Volcano plots: button to reset volcano annotations</li>
</ul>
<b>v0.6.6.1 Mar 17, 2017</b>
<ul>
<li>Export-tab: PCA loadings can be exported as Excel sheet (by Ozan Aygun).</li>
<li>PCA-tab: New PCA loadings plot (by Ozan Aygun).</li>
<li>Export-tab: included button for PCA loadings in \"toggle all\".</li>
<li>Heatmap-tab: Default vaules for row/column fon size read from \"plotparams\", if defined.</li>
</ul>
<b>v0.6.6 Mar 17, 2017</b>
<ul>
<li>Fixed the \"incorrect number of dimensions\"-error in the table preview tab, if only a single annotation column is present.</li>
<li>Prevented the automatic switch to the \"Summary\"-tab after changing the filter.</li>
<li>Related to the previous point, the result filter is now implemented as observer rather than a reactive function.</li>
<li>Summary-tab: fixed the workflow box showing NA when selecting filter \"none\" or \"top.n\".</li>
<li>Dynamic UI elements will not switch back to \"One-sample modT\" after running an analysis.</li>
<li>Table-tab: switched to DT package.</li>
</ul>
<b>v0.6.5 Mar 7, 2017</b>
<ul>
<li>Fixed a bug that resulted in not listing all saved session for a user.</li>
<li>Worked on the filenames of exported RData and Excel files.</li>
<li>modF: In case of too many missing values the test would not return a p-value which resulted in NA for the enumber of significant hits on the summary page.</li>
</ul>
<b>v0.6.4 Mar 6, 2017</b>
<ul>
<li>Summary tab: number of significant hits are now reported correctly.</li>
<li>Summary tab: Missing value distribution after log-transformation shown correctly.</li>
<li>Changed cluster method from \'complete\' to \'ward\'.</li>
<li>Fixed a bug that was triggered if a project is defined and shared in \'user-roles.txt\' but has been deleted from the server.</li>
</ul>
<b>v0.6.3 Feb 2, 2017</b>
<ul>
<li>Commited to GitHub for debugging purposes. Do not use this version!</li>
<li>Re-organization of UI elements when setting up the analysis.</li>
<li>Implementation of SD filter across all samples.</li>
</ul>
<b>v0.6.2 Jan 31, 2017</b>
<ul>
<li>UI elements for setting up an anlysis workflow are now dynamically generated, e.g. if reproducibility filter is chosen, only "One-sample modT" or "none" can be selected.</li>
<li>Reproducibility filter: users can choose bewteen (predefined) alpha-values.</li>
<li>Increased number of colors by 60 (85 total).</li>
<li>Correlation matrix: increased the size of exported heatmap to 12x12 inches.</li>
<li>Multiscatter: increased number of digits to three.</li>
<li>Some more error handling when exporting analysis results.</li>
<li>Previously saved sessions are not deleted anymore, if checkbox "save session" is not enabled.</li>
</ul>
<b>v0.6.1 Jan 12, 2017</b>
<ul>
<li>Session managment: Added possibility to delete saved sessions and to choose whether to save a session on the server in the first place.</li>
<li>User role managment (alpha status): A project saved on the server has an owner and (optional) collaborators. Collaborators can \"see\" projects they are assigned to in the dropdown menu \"Saved sessions\".</li>
</ul>
<b>v0.6.0 Jan 4, 2017</b>
<ul>
<li>Switched to <a href="https://rstudio.github.io/shinydashboard/">Shiny Dashboards</a>.</li>
<li>Extented PCA analysis using the <i>ChemometricsWithR</i> R package.</li>
</ul>
<b>v0.5.4  Dec 27, 2016</b>
<ul>
<li>Filter values and plotting parameters are now restored after session import (except for volcano plot...).</li>
<li>Changed visual style of volcano plots.</li>
</ul>
<b>v0.5.3 Dec 21, 2016</b>
<ul>
<li>Minor fixes due to shiny update.</li>
<li>User can now specify label names used to create file and session names when  exporting results. Initial values are taken from filenames of the input and experimental design file.</li>
<li>Experimental design file is now part of the results.</li>
</ul>
<b>v0.5.2 Dec 1, 2016</b>
<ul>
<li>Rudimentary support of \'gct\' files, i.e. files can be imported by ignoring the first two lines (gct header). </li>
<li>Figured out the issue with  the 2-sample T-test volcanos. The functions in \'limma\' always report fold changes group factor variable \'0\'. The original \'moderated.t\' alphabetically orders the class names and then converts class names to factors. First class name will become zero. I make sure that class names are alphabeticaly sorted before calling \'moderated.t\'.</li>
</ul>
<b>v0.5.1 Nov 26, 2016</b>
<ul>
<li><mark>BUG: </mark>Reverted the indication of direction in volcano plots for <b>2-sample tests</b>. The direction was inferred from the sign of \'logFC\' returned by function \'topTable\' (limma) which cannot be used to do that.</li>
<li>Updated shiny R package from 0.12/0.13.2 to 0.14.2 resulting in some minor changes in the <i>look and feel</i> of the app. Code needed some adaptions (navbarPage, navbarMenu) to run poperly with 0.14.2 version.</li>
<li>Outsourced HTML instructions to a separate file using Shiny-module framework.</li>
<li>Changed how heatmap dimensions are determined to better show very large and very small heatmaps.</li>
<li>Scaling of heatmap done after clustering.</li>
</ul>
<b>v0.5.0 Nov 7, 2016</b>
<ul>
<li>Exported sessions are saved on the server and can be re-imported. Each user has its own folder on ther server in which an R-sessions file is stored.</li>
<li>Non-unique entries in the id column are made unique, e.g. \'Abl\', \'Abl\' -> \'Abl\', \'Abl_1\'. Empty entries will be replaced by \'X\', e.g. \'Abl\', \'\', \'\' -> \'Abl\', \'X\', \'X_1\'.</li>
</ul>
<b>v0.4.5 Sep 1, 2016</b>
<ul>
<li>Multiscatter: log-transformed values wil be used if log-transformation has been applied.</li>
<li>For each user a new folder on the server is created. Every session that gets exported will be saved there.</li>
<li>A copy of the original data file will be part of the results (zip-file).</li>
</ul>
<b>v0.4.4 Aug 19, 2016</b>
<ul>
<li>New \'Export\'-tab to download a zip-file containing:
 <ul>
   <li>all figures (pdf).</li>
   <li>result table (xlsx).</li>
   <li>session file (Rdata) which can be imported back into the app.</li>
   <li>parameter file (txt)</li>
 </ul>
<li>Directionality of two-sample test is now indicated in the volcano plots.</li>
<li>Error handling for two-component normalization.</li>
<li>Profile plots under \'QC\'-tab</li>
</ul>
<b>v0.4.3 Aug 16, 2016</b>
<ul>
<li>Session export/import.</li>
<li>"#VALUE!"-entries from Excel can be handeled now.</li>
<li>Fixed bug causing PDF export of heatmap with user defined max. values to crash.</li>
</ul>
<b>v0.4.2 Jul 21, 2016</b>
<ul>
<li><mark>BUG:</mark> Bugfix in 2-sample test that occured whenever the names of different groups defined the experimental design file started with the same series of characters, e.g. \'ABC\' and \'ABCD\'.</li>
</ul>
<b>v0.4.1 Jul 1, 2016</b>
<ul>
<li>Novel tab summarizing the analysis.</i>
<li>Data can now be log-transformed, e.g. for MaxQuant LFQ results.</li>
<li>Added option to skip testing, e.g. for PCA analysis.</li>
<li>User can specify principle components in the PCA scatterplot.</li>
</ul>
<b>v0.4 Jun 29, 2016</b>
<ul>
<li>Integration of moderated F statistics</li>
<li>Disabled column-based clustering one-sample and two-sample tests if multiple groups are being compared.</li>
</ul>
<b>v0.3 Mar 11, 2016</b>
<ul>
<li>Data normalization.</li>
<li>Reproducibility filter.</li>
<li>Upload/download of experimental design files.</li>
<li>Download of native Excel files.</li>
<li>Integration of the Javascript D3-based plotly library.</li>
</ul>
<b>v0.2 Feb 23, 2016</b>
<ul>
<li>Working version on server.</li>
</ul>
<b>v0.1 Dec 20, 2015</b>
<ul>
<li>First prototype.</li>
</ul>
</font>'
#,sep='')

        ## render HTML
        output$html <- renderText({
            if(!is.null(global.input$file)) return()

            HTML(txt)
        })
    }

    ##@#########################################
    ## id column / exp design template
    if(what == 'id'){

        txt <- paste('<br><br><p><font size=\"4\"><b>Group assigment</b></br>
Download a template of an experimental design file. You can open this file in Excel and define the groups you want to compare. Replicate measurements have to be grouped under a single name in the \'Experiment\'-column. You can use the \'Group\' column to denote group-wise normalization. If you are not performing group-wise normalization, you may leave this as NA. If you would like to exclude a sample from analysis, you may leave both Experiment and Group as NA. <mark>Please don\'t use special characters, like blanks or any punctuation, when defining these names!</mark> </font></p>
<br><p><font size=\"4\"><b>Select ID column</b></br>
Choose a column from the list on the left that contains <b>unique</b> identifiers for the features in the data table. If the enntries are not unique, uniqueness will enforces by appending \"_1\". Preferably, IDs should be unique protein accession numbers (e.g. <font face=\"Courier\">NP_073737</font>) or a combination of protein accession and residue number in case of PTM analysis (e.g. <font face=\"Courier\">NP_073737_S544s _1_1_544_544</font>).</p>  
<br><p><font size=\"4\"><b>Automatic retrieval of gene symbols</b></br>
If the ID column contains <a href=\"http://www.uniprot.org/\" target=\"_blank_\">UniProt</a> or <a href=\"https://www.ncbi.nlm.nih.gov/refseq/\" target=\"_blank_\">RefSeq</a> accession numbers, the software will try to map those ids to gene symbols. Currently, mapping of following organisms is supported:
<ul>
<li>human (<i>Homo sapiens</i>)</li>
<li>mouse (<i>Mus musculus</i>)</li>
<li>rat (<i>Rattus norvegicus</i>)</li>
<li>zebrafish (<i>Danio rerio</i>)</li>
</ul>
</font></p>')

        ## render HTML
        output$html <- renderText({

            if(global.param$analysis.run) return()
            if(global.param$id.done) return()
            if(!global.param$file.done) return()
            
            #if(is.null(global.input$id.col)) return() ## start page

            #if(global.input$id.col > 0 && !is.null(global.param$id.col.value)) return() ## after id column is choosen

            HTML(txt)
        })
    }
    #####################################################################
    ## upload of experimental design file
    if(what == 'ed'){

        txt <- paste('<br><br><p><font size=\"4\">Please upload the experimental design file that you have created using the upload button on the left. <mark>If using your own experimental design file (not created from the template on the previous screen), the first column MUST contain the sample (column) names and must match the names in your table exactly!</mark> It is optional for your experimental design file to include the ID column. <mark> All other columns must be present in your experimental design file </mark> </p></font></p>')

        ## render HTML
        output$html <- renderText({

            if(global.param$analysis.run) return()
            if(!global.param$file.done) return()

            if( is.null(global.param$id.col.value) ) return()
          
            #if(global.input$id.col ==0) return()
            if(global.param$file.gct3) return()
            #if(global.param$id.done) return()  
            if(global.param$grp.done) return()

            HTML(txt)
        })
    }
    #####################################################################
    ## gct v3
    if(what == 'gct3'){
      
      txt <- paste('<br><p><font size=\"4\">Choose the annotation column to use as class vector for marker selection.</p></font></p>')
      
      ## render HTML
      output$html <- renderText({
        
        if(global.param$analysis.run) return()
        if(!global.param$file.gct3) return()
        
        #if(global.param$id.done) return()  
        if(global.param$grp.done) return()
        
        HTML(txt)
      })
    }
    
    
    
    
    ## ####################################################################
    ## analysis
    if(what == 'ana'){

        txt <- paste('<font size=\"4\">
        <p><h3>Intensity data</h3>Check the box if you are using raw or log-transformed intensity data. Only the relevant normalization methods and statistical tests for intensity data will now show under \'Data normalization\' and \'Select test\'.</p>
<p><h3>Log-transformation</h3>Apply log transformation to the data. <b>Data MUST be log-transformed for statistical testing!</b> If you select "none," we assume you have already log-transformd your data.</p>
<p><h3>Normalize per group</h3>If enabled the normalization will be performed within a particular group (Median, Median-MAD, Quantile, VSN). For Median and Median-MAD normalization, the group-level median of sample medians is added to each normaized data value.</p>
<p><h3>Data normalization</h3>You can apply different normalization methods to the data prior to testing. The methods are applied for each sample (column) separately, except for \'Quantile\' and \'VSN\' normalization which take the entire matrix into account.</p>
<p>
<ul>
<li><b>Median</b>: Subtract the sample median from each value (centering). After normalization, all samples have a median of zero. Intended to be used with <b>log-transformed ratios</b>.</li>
<li><b>Median (non-zero)</b>: Subtract the sample median from each value and add back the median of all sample medians (which will be the common sample median after normalization). Intended to be used with <b>log-transformed intensities</b>.</i> 
<li><b>Median-MAD</b>: Subtract the sample median and divide by sample MAD (centering plus scaling). After normalization, all samples have a median of zero. Intended to be used with <b>log-transformed ratios</b>.</li>
<li><b>Median-MAD (non-zero)</b>: Subtract the sample median and divide by sample MAD, and add back the median of all sample medians (which will be the common sample median after normalization). Intended to be used with <b>log-transformed intensities</b>.</li>
<li><b>Upper quartile</b>: Subtract the sample\'s 75th percentile from each value. Intended to be used with <b>log-transformed intensities</b>.</li>
<li><b>2-component</b>: Use a mixture-model approach to separate non-changing from changing features. Z-score all features using the mean and standard deviation of the non-changing features. Intended to be used with <b>log-transformed ratios</b>.</li>
<li><b>Quantile</b>: Transform the data such that the quantiles of all sample distributions are the equal. <b>Use with caution as this type of normalization can remove potentially meaningful outliers from the data</b>.</li>
<li><b>VSN</b>: Variance stabilizing normalization. Intended to be used with <b>raw intensity values</b>.</li>
<li><b>none</b>: The data will be taken as is. Use this option if the data has already been normalized.</li>
</ul>



<h3>Filter data</h3>

<b>You can filter the data by p-value (non-adjusted or adjusted) and change the p-value cutoff once running the analysis (on the next screen).</b>

<br><br><b>Missing data:</b><br>
This represents the maximum number of missing values allowed. For example: A value of 70 means a feature may have up to 70% missing values (must be quantified in at least 30% of samples). If you do not wish to perform missing value filtering, please leave the value at 100 (99 for intensity-based data). For intensity data, the missing data rate is capped at 99% to prevent downstream statistical testing from throwing an error.

<br><br><b>Reproducibility:</b><br>
Remove features that were not reproducibly quantified across replicate measurements of a group. For duplicate measurements a Bland-Altman Filter of 99.9% (+/-3.29 sigma) will be applied. For more than two replicate measurements per group a generalized reproducibility filter is applied which is based on a linear mixed effects model to model the within-group variance and between-group variance (See \'MethComp book (pp 58-61). <i>Comparing Clinical Measurement Methods</i> by Bendix Carstensen\' for more details). You can inspect the results of the filtering step in the multiscatter plot under the \'QC\'-tab as well as in the interactive scatterplots. Data points removed prior to testing will be depicted in blue. <b>This type of filter is applied separately to each group.</b>  <b>This filter is only appropriate for ratio data.</b>

<br><br><b>StdDev:</b><br>
Remove features with low standard deviation across all samples. Only useful if applied to sample cohorts that were quantified against a common reference. The percentile <b><i>P</i></b> you specify 
in the slider refers to the <b><i>P</i></b> percent of features having the <b>lowest standard deviation</b> across sample columns which will be <b>excluded prior to analyis</b>.
Using this type of filter is useful to explore result of unsupervised clustering of the data without running a statistical test.

<br><h3>Select test</h3>You can choose between a one-sample, two-sample moderate T-tests, moderated F-test or no testing.
<ul>
<li><b>One-sample mod T</b>: For each group test whether the group mean is significantly different from zero. Only meaningful to <b>ratio data</b>!</li>
<li><b>Two-sample mod T</b>: For each possible pairwise comparison of groups test whether the group means are significantly different from each other.</li>
<li><b>mod F</b>: Test whether there is a significant difference between any of the defined groups. Should be used if more than 2 groups are being compared. </li>
<li><b>none</b>: Don\'t do any test. Useful for exploratory data analysis such as PCA.</li>
</ul>
<br></font></p>')

        ## render HTML
        output$html <- renderText({
             if(global.param$analysis.run) return()
##             if( !is.null(error$msg) ) return()

             if(global.param$grp.done == F) return()
             if(!is.null(global.input$run.test))
                 if(global.input$run.test > 0) return()

             HTML(txt)
        })
    }

    ## ####################################################################
    ## analysis
    if(what == 'res'){

        txt <- paste('<p><font size=\"4\">This page allows you to interactively explore the results of your analyis. On the left you can choose between different filters, the results will be updated immediately. <mark> Note that by scrolling down you can see more options, such as setting the type of p-value and the cutoff. </mark> The filter that you specify applies to all tabs (\'Heatmap\', \'Volcanos\', ...), except the \'QC\' which shows the entire dataset. You can change the appearance of the heatmap by modifying the parameters below, you can select points shown in the Volcano plots and browse through the result table.</font></p><br>')

        ## render HTML
        output$html <- renderText({
        ##     if( !is.null(error$msg) ) return()

            if(global.param$grp.done == F) return()
            if(!is.null(global.input$run.test))
                if(global.input$run.test == 0) return()

            HTML(txt)
        })
    }


} ## end printHTML

