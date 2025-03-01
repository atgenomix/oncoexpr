#' @title KEGG Enrichment Analysis
#' @description a set of genes for enrichment analysis with KEGG database
#' @param gene_list_ a set of genes with specific expression level
#' @param save_path_ what path you want to save the files
#' @param save_filename_ what name you call the files
#' @param showCategory_ the number of GO term you want to show in doplot, default: 10
#' @return the dotplot for KEGG enrichment
#' @export
#' 
kegg_enrich_dotplot <- function(gene_list_, save_path_=NULL, save_filename_=NULL, showCategory_ = 10){
  save_parameter <- paste0(save_path_, save_filename_)
  gene_ids <- bitr(gene_list_, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  kegg_enrich_result <- enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
  kegg_dotplot_ <- dotplot(kegg_enrich_result, showCategory = showCategory_ , font.size = 15)
  if(length(save_path_)!=0){
    ggsave(   save_parameter, 
              plot = kegg_dotplot_, 
              width = 32, 
              height = 28, 
              units="cm", 
              dpi = 300)
  }else{
    print("not saving the figure")
  }
  
  return(kegg_dotplot_)
}
