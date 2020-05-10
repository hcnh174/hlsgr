#' Title
#'
#' @param pathway
#' @param genes
#' @param kegg_dir
#' @param bg_col
#'
#' @return
#' @export
#'
#' @examples
#' hlsgr::plotKeggPathway(pathway='hsa05200', genes=c(672, 811, 2065, 2260, 201163, 2321, 6416, 5727, 5888, 5925, 7157), kegg_dir='c:/workspace/hlsg/data/kegg')
plotKeggPathway <- function(pathway, genes, kegg_dir, bg_col='firebrick1')
{
  expr <- data.frame(expr=rep(1, length(genes)))
  rownames(expr) <- genes
  species <- substr(pathway, 1, 3)
  pathway_id <- substr(pathway, 4, 100)
  outfile <- tempfile(pattern=paste0(pathway,'_profile_bg'), fileext='.png')
  xml <- KEGGprofile::parse_XMLfile(pathway_id=pathway_id, species=species, database_dir=kegg_dir)
  KEGGprofile::plot_profile(expr, pathway_name=pathway, KEGG_database=xml, database_dir=kegg_dir,
                            bg_col=bg_col, type='bg', result_name=outfile)
  return(outfile)
}

