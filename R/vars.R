
####################################
###Variables for figure formating###
####################################


#set the global ggplot theme (full frame)
#' @export
theme_full <- ggplot2::theme_bw() + ggplot2::theme(axis.text = ggplot2::element_text(size=15),
                             axis.title = ggplot2::element_text(size=16),
                             axis.line = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(size=1.5),
                             axis.ticks = ggplot2::element_line(size=1.5),
                             plot.title = ggplot2::element_text(size = 16, hjust =0.5, face="bold"),
                             panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank())
#set the global ggplot theme (half frame)
#' @export
theme_half <- ggplot2::theme_bw() + ggplot2::theme(axis.text = ggplot2::element_text(size=15),
                                 axis.title = ggplot2::element_text(size=16),
                                 axis.line = ggplot2::element_line(size=0.8),
                                 panel.border = ggplot2::element_blank(),
                                 axis.ticks = ggplot2::element_line(size=1.5),
                                 plot.title = ggplot2::element_text(size = 16, hjust =0.5, face="bold"),
                                 panel.grid.major = ggplot2::element_blank(),
                                 panel.grid.minor = ggplot2::element_blank())

# Defien a color scheme, based on ggsci_NEJM panel, for the paper
#' @export
colList <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")

# Color scheme for CyTOF
#' @export
colorMap <- c(high = colList[3], low = colList[6], M = colList[1], U=colList[2],
              `high CLL-PD` = colList[3], `low CLL-PD` = colList[6])

# Format conditions, one line version
#' @export
condiMap <- c(DMSO = "DMSO", CpG = "CpG ODN", untreated = "Day0",
              CpG_Ever = "CpG ODN + Everolimus", Ever = "Everolimus")

# Format conditions, two line version
#' @export
condiMapS <- c(DMSO = "DMSO", CpG = "CpG ODN", untreated = "Day0",
              CpG_Ever = "CpG ODN\n+ Everolimus", Ever = "Everolimus", day0 = "untreated\n(day0)")

# Format cell type clusters
#' @export
clustMap <- c(CLL_proliferating = "proliferating CLL cells", CLL_resting = "non-proliferating CLL cells", Apoptotic = "dead/apoptotic cells",
                Tcell = "T cells", Neutrophil = "myeloid cells")

