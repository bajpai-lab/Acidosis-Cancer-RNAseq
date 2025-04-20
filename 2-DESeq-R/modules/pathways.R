toppathwaysIPA <- function (
    input,
    tag,
    limit=30,
    exclude_from_top_by_index = c(),
    height=9
){
  print(paste0("Reading ",input))
  # Reading in CSV
  IPA_ALL <- read.csv(input)
  
  # Removing errored / incomplete rows
  IPA_S <- IPA_ALL %>%
    filter(complete.cases(.)) %>%
    filter(combined != "#N/A")
  
  # Ensuring relevant columns are numbers
  IPA_S$combined = as.numeric(as.character(IPA_S$combined)) 
  IPA_S$nlogp = as.numeric(as.character(IPA_S$nlogp)) 
  
  # Sorting now that data is valid type
  IPA_S <- IPA_S %>%
    arrange(combined)
  
  # First removing entries requested to be removed
  if (length(exclude_from_top_by_index) > 0){
    IPA_S <- IPA_S[-exclude_from_top_by_index,]
  }
  
  # Then cap off to limit
  if (dim(IPA_S)[1] - limit > 0){
    # Above checks that there are enough entries to bother cutting the set
    IPA_S <- IPA_S[c(1:limit),]
  }
  
  IPA_GENES <- unlist(strsplit(str_replace_all(paste0(IPA_S$Genes, collapse=""),c("(\\[|, )([A-Z0-9a-z]+)( +\\-).+?" = "PseqL_\\2_PseqR", "_PseqR.*?PseqL_" = ",", "^.*PseqL_"="","_PseqR.*$"="")),","))
  IPA_GENES <- IPA_GENES[!duplicated(IPA_GENES)]
  
  IPA_S <- IPA_S %>% mutate(combined=combined)
  
  IPA_S <- IPA_S[complete.cases(IPA_S),]#remove any rows with NA
  
  IPA_P <- ggplot(IPA_S, aes(x = reorder(Ingenuity.Canonical.Pathways, -combined), y = -combined, main="IPA pathways inhibited by pHe=6.4 vs. 7.4")) +
    geom_bar(stat = "identity", fill="#13aeef",width = 0.75) +
    geom_text(aes(label=paste0("p = ",format(10^nlogp, digits=3, scientific=TRUE))), hjust=1.05) +
    coord_flip(clip = "off") +
    scale_x_discrete(position = "top") +
    theme(axis.text.x = element_text(color="#000",
                                     size=11, angle=0),
          axis.text.y = element_text(color="#000",
                                     size=11, angle=0)) +
    labs(title=paste0("IPA: ",tag),x="")+
    scale_y_reverse(name="Combined Score")+
    theme_classic()+
    theme(
      plot.margin = margin(t = 1,  # Top margin
                           r = 1,  # Right margin
                           b = 1,  # Bottom margin
                           l = 3,  # Left margin
                           unit = "cm"))
  
  print(paste0("Exporting ","results/IPATopPathways_",tag,".pdf"))
  ggsave(filename=paste0("results/IPATopPathways_",tag,".pdf"), IPA_P, width=9, height=height)
  return (IPA_P)
}
