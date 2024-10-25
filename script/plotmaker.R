#!/home/pgs8070/miniconda3/bin/Rscript
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
main                <- function(){
    Peak_scatter('Early', 'ChIP_CRE')
	Peak_scatter('Full', 'ChIP_CRE')
	Peak_scatter('Late', 'ChIP_CRE')
	Peak_scatter('Mid', 'ChIP_CRE')
}

Peak_scatter        <- function(stage, peak_type){
    data            <- read.table(paste0('../output/pvalue_maker/results_',peak_type,'_',stage,'.txt'), header = T, comment.char= '')
    iteration       <- as.integer(sub(".*_", "", colnames(data)[3]))
    colnames(data)  <- c("ID", "real_peak", "Pval", "sum_peak")
    data$sum_peak   <- data$sum_peak / iteration
    data$Pval       <- (iteration - data$Pval) / iteration
    data$Pval[data$Pval == 0]   <- 1/iteration
    data$adjPval    <- p.adjust(data$Pval, method='BH')
    data$fcVal      <- data$real_peak / data$sum_peak
    data$size       <- rank(data$real_peak) / 4
    p_threshold     <- 5e-3
    fc_threshold    <- 2

    insig           <- data[-which(data$adjPval < p_threshold & data$fcVal > fc_threshold), ]
    significantpval <- data[data$adjPval < p_threshold,]
    significant     <- significantpval[significantpval$fcVal > fc_threshold, ]
    write.table(significant, file = paste0('/home/pgs8070/o828_GWAS_practice/output/plotmaker/Peak_scatter_',peak_type,'_',stage,'_',iteration,'_significant.txt'))
    write.table(data, file = paste0('/home/pgs8070/o828_GWAS_practice/output/plotmaker/Peak_scatter_',peak_type,'_',stage,'_',iteration,'_all.txt'))

    g1  <- ggplot() +
        geom_point(aes(data$fcVal, -log10(data$adjPval), size=data$real_peak), shape=1, col='darkgrey') +
        geom_point(aes(significant$fcVal, -log10(significant$adjPval), size=significant$real_peak, col='red')) +
        geom_hline(yintercept=-log10(p_threshold), linetype="dashed", color = "darkgrey") +
        geom_vline(xintercept=fc_threshold, linetype="dashed", color = "darkgrey") +
        labs(x='EC-specific cRE enrichment(FC over expectation)', y='-log10(adjusted P)') +
        theme_classic()

    # insig2          <- data[-which(data$adjPval < p_threshold & data$fcVal > fc_threshold), ]
    # significantpval <- data[data$adjPval < p_threshold,]
    # significant     <- significantpval[significantpval$fcVal > fc_threshold, ]
    # write.table(significant, file = paste0('/home/pgs8070/o828_GWAS_practice/output/plotmaker/Peak_scatter_',peak_type,'_',stage,'_',iteration,'_significant.txt'))
    # write.table(data, file = paste0('/home/pgs8070/o828_GWAS_practice/output/plotmaker/Peak_scatter_',peak_type,'_',stage,'_',iteration,'_all.txt'))

    # others_list     <- c("A_body_shape_index", "Abdominal_aortic_aneurysm", "Aspartate_aminotransferase_to_alanine_aminotransferase_ratio", "Balding_type_1", "Birth_weight", "Blond_vs._brown_black_hair_color", "Lung_function_(FEV1_FVC)", "Offspring_birth_weight", "Serum_alkaline_phosphatase_levels", "Testosterone_levels")
    # hema_list       <- c("Basophil_count", "Coronary_artery_disease_or_fibrinogen_levels_(pleiotropy)", "Neutrophil_percentage_of_white_cells", "Plateletcrit", "Platelet-to-lymphocyte_ratio")
    # disease_list    <- c("Colorectal_cancer", "Testicular_germ_cell_tumor", "Vertical_cup-disc_ratio_(multi-trait_analysis)", "Rheumatoid_arthritis_or_type_1_diabetes")

    # others_data     <- data[data$ID %in% others_list,]
    # hema_data       <- data[data$ID %in% hema_list,]
    # disease_data    <- data[data$ID %in% disease_list,]
    

    # g1  <- ggplot() +
    #     geom_point(aes(insig2$fcVal, -log10(insig2$adjPval), size=insig2$real_peak), shape=1, col='darkgrey') +
    #     geom_point(aes(hema_data$fcVal, -log10(hema_data$adjPval), size=hema_data$real_peak, col='red')) +
    #     geom_point(aes(others_data$fcVal, -log10(others_data$adjPval), size=others_data$real_peak, col='green')) +
    #     geom_point(aes(disease_data$fcVal, -log10(disease_data$adjPval), size=disease_data$real_peak, col='blue')) +
    #     geom_text(aes(hema_data$fcVal, -log10(hema_data$adjPval), label=hema_data$ID, col='red'), vjust=1, size=3, hjust=0, nudge_x=0.05) +
    #     geom_text(aes(others_data$fcVal, -log10(others_data$adjPval), label=others_data$ID, col='green'), vjust=1, size=3, hjust=0, nudge_x=0.05) +
    #     geom_text(aes(disease_data$fcVal, -log10(disease_data$adjPval), label=disease_data$ID, col='blue'), vjust=1, size=3, hjust=0, nudge_x=0.05) +
    #     geom_hline(yintercept=-log10(p_threshold), linetype="dashed", color = "darkgrey") +
    #     geom_vline(xintercept=fc_threshold, linetype="dashed", color = "darkgrey") +
    #     labs(x='EC-specific cRE enrichment(FC over expectation)', y='-log10(adjusted P)') +
    #     theme_classic()

    # g2  <- ggplot() +
    #     geom_point(aes(data$fcVal, -log10(data$Pval), size=data$real_peak), shape=1, col='darkgrey') +
    #     geom_point(aes(label_data$fcVal, -log10(label_data$Pval), size=label_data$real_peak, col='red')) +
    #     geom_text(aes(label_data$fcVal, -log10(label_data$Pval), label=label_data$ID, col='red'), vjust=1, size=3, hjust=0, nudge_x=0.05) +
    #     geom_hline(yintercept=-log10(p_threshold), linetype="dashed", color = "darkgrey") +
    #     geom_vline(xintercept=fc_threshold, linetype="dashed", color = "darkgrey") +
    #     labs(x='EC-specific cRE enrichment(FC over expectation)', y='-log10(P)') +
    #     theme_classic()

    output_directory    <- paste0('../plot/Peak_scatter_',peak_type,'_',stage,'_',iteration,'.pdf')
    pdf(output_directory, width=8, height=8, pointsize=4)
    plot(g1)
    # plot(g2)
    dev.off()
}

main()