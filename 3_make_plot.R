#!/home/ajl1213/anaconda2/envs/r_plot/bin/R



library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(viridis)


dir.create('Plots')

resolution <- '20kb'

sample_list <- c(
'hESC',
'ME',
'EC4',
'EC12',
'EC24',
'EC48',
'nonEC12',
'nonEC48'
)


##================================================================================================
## siginter type
data1 <- read.table('res_fithic/siginter.20kb.merged.txt', header=T)

##
pdf('Plots/siginter.nsample.bar.pdf')
barplot(table(data1$nsample), xlab='number of samples', ylab='number of sig inter')
dev.off()

## 3d piechart
library(plotrix)

#
table <- data.frame(table(data1$inter_type))
table$frac <- round(table$Freq / sum(table$Freq) * 100, digits=2)
table$inter_type <- paste0(table$Var1, ' (', table$frac, '%)')

#
col_list <- c(brewer.pal(8, 'Accent')[c(5, 3)])
pdf('Plots/siginter.inter_type.3d_pie.pdf')
pie3D(rev(table$frac), labels=rev(table$inter_type), radius=1,  explode=0.08, col=rev(col_list))
dev.off()


##================================================================================================
## inter class barplot

my_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(panel.background = element_rect(fill='white', color='black', linetype='solid')) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(legend.position='bottom')

## static inter
merged_df <- NULL
for (sample in sample_list){
	print(sample)

	data1 <- read.table(paste0('res_fithic/',sample,'.',resolution,'.sig_inter.filtered.anno.txt'), header=T)
	tmp_df <- data.frame(inter_id=data1$inter_id, inter_dir=data1$inter_dir, inter_class=data1$inter_class, inter_type=data1$inter_type, dist=data1$dist, sample=sample)

	merged_df <- rbind(merged_df, tmp_df)
}

##
count_table <- as.matrix(table(merged_df$sample, merged_df$inter_class))

get_frac <- function(x){x / sum(x) * 100}
frac_mat <- t(apply(count_table, 1, get_frac))

df <- as.data.frame(as.table(frac_mat))
colnames(df) <- c('sample','inter_class','frac')

df$sample <- factor(df$sample, levels=sample_list)
df$inter_class <- factor(df$inter_class, levels=c('peak_prom', 'peak_none', 'none_none'))

print(df)

#
colors <- rev(c('lightgrey','darkgrey','black'))

g1 <- ggplot(df) +
    geom_bar(aes(x=sample, y=frac, fill=inter_class), stat='identity') +
    scale_fill_manual(values=colors) +
    labs(x='', y='Composition of siginter') +
    my_theme

pdf('Plots/siginter.inter_class.bar.pdf')
plot(g1)
dev.off()


##================================================================================================
## peak composition of sig inters
stage_list <- c('hESC','Mesoderm','EarlyEC','MidEC','LateEC','FullEC','MidNonEC','FullNonEC')

for (datatype in c('static','variable')){
	print(datatype)

	anno_df <- read.table('anno_res/sig_inter.peak_label.txt', header=T)
	anno_df <- anno_df[anno_df$peak_label %in% stage_list, ]
	anno_df <- anno_df[which(anno_df$inter_type==datatype), ]

	count_table <- as.matrix(table(anno_df$sample, anno_df$peak_label))

	get_frac <- function(x){x / sum(x) * 100}
	frac_mat <- t(apply(count_table, 1, get_frac))

	df <- as.data.frame(as.table(frac_mat))
	colnames(df) <- c('sample','peak_label','frac')

	df$sample <- factor(df$sample, levels=sample_list)
	df$peak_label <- factor(df$peak_label, levels=stage_list)

	print(df)

	#
	colors <- c(
	brewer.pal(8, 'Dark2')[4],
	brewer.pal(9, 'YlOrRd')[5],
	brewer.pal(9, 'PuBu')[5],
	brewer.pal(9, 'PuBu')[6],
	brewer.pal(9, 'PuBu')[7],
	brewer.pal(9, 'PuBu')[8],
	brewer.pal(9, 'Greys')[6],
	brewer.pal(9, 'Greys')[8]
	)

	g1 <- ggplot(df) +
	    geom_bar(aes(x=sample, y=frac, fill=peak_label), stat='identity') +
	    scale_fill_manual(values=colors) +
	    labs(x='', y='Composition of stage-specific cREs', main=datatype) +
	    my_theme

	pdf(paste0('Plots/siginter.', datatype, '.peak_composition.bar.pdf'))
	plot(g1)
	dev.off()
}


## log2 enrichment
# background
data1 <- read.table('anno_res/sig_inter.merged.peak_label.txt', header=T)
data1 <- data1[data1$peak_label %in% stage_list, ]
bg_table <- table(data1$peak_label) 
bg_total <- sum(bg_table)
bg_pct <- bg_table/bg_total

#
for (datatype in c('static','variable')){
	print(datatype)

	anno_df <- read.table('anno_res/sig_inter.peak_label.txt', header=T)
	anno_df <- anno_df[anno_df$peak_label %in% stage_list, ]
	anno_df <- anno_df[which(anno_df$inter_type==datatype), ]

	fc_mat <- NULL
	pval_mat <- NULL
	for (sample in sample_list){

		tmp_df <- anno_df[which(anno_df$sample==sample), ]

		test_table <- table(tmp_df$peak_label)
		test_total <- sum(test_table)
		test_pct <- test_table/test_total

		fc_list <- test_pct/bg_pct
		fc_mat <- rbind(fc_mat, fc_list)

		pval_vec <- c()

		#
		for (peak_label in stage_list){
			obs <- test_table[[peak_label]]
			bg <- bg_table[[peak_label]]

			p <- phyper(q = obs - 1,
			m = bg,
			n = bg_total - bg,
			k = test_total,
			lower.tail = FALSE
			)

			pval_vec <- c(pval_vec, p)
		}
		pval_mat <- rbind(pval_mat, pval_vec)
	}

	rownames(fc_mat) <- sample_list
	fc_mat <- fc_mat[, stage_list]

	rownames(pval_mat) <- sample_list
	colnames(pval_mat) <- stage_list

	pval_mat[pval_mat==0] <- .Machine$double.xmin
	logP_mat <- -log10(pval_mat)

	print(fc_mat)
	print(logP_mat)


	##
	ncolors <- 100 
	heatmap_colors <- colorRampPalette(c(brewer.pal(9, 'BuPu'), brewer.pal(9, 'BuPu')[9], 'black'))(ncolors)

	min_val <- 0
	max_val <- 300
	logP_mat[logP_mat < min_val] <- min_val
	logP_mat[logP_mat > max_val] <- max_val

	pdf(paste0('Plots/siginter.', datatype, '.peak_enrichment.heatmap.pdf'))
	pheatmap(
		logP_mat, color= heatmap_colors,
		cluster_rows=F, cluster_cols=F,
		show_colnames=T, show_rownames=T,
		breaks=seq(min_val, max_val, length.out=ncolors),
		border_color = NA
	)
	dev.off()

}


##================================================================================================
## sig inter dist boxplot

## by sample
df <- NULL
for (sample in sample_list){
        print(sample)

        data <- read.table(paste0('res_fithic/', sample, '.', resolution, '.sig_inter.filtered.txt'), header=T)
        tmp_df <- data.frame(sample=sample, dist=data$dist)

        df <- rbind(df, tmp_df)
}
df$sample <- factor(df$sample, levels=sample_list)

#
my_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(panel.background = element_rect(fill='white', color='black', linetype='solid')) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme(legend.position='bottom')

p1 <- ggplot(df) +
        geom_boxplot(aes(x=sample, y=dist), notch=T, outlier.shape=NA, width=0.9) +
        ylim(0, 1200000) +
        labs(x=NULL, y='Distance for significant interactions') +
        my_theme

pdf('Plots/siginter.dist.by_sample.box.pdf')
plot(p1)
dev.off()


## merged siginter
data1 <- read.table('res_fithic/siginter.20kb.merged.txt', header=T)

pdf('Plots/siginter.dist.merged.hist.pdf')
hist(data1$dist, breaks=50, freq=F, col='lightblue', xlab='Distance (bp)', main=NULL )
lines(density(data1$dist, adjust=10), col='navy', lwd=2)
abline(v=mean(data1$dist), col='darkgrey', lwd=2, lty='dashed')
dev.off()

print(mean(data1$dist))

# static vs variable
point_df <- data.frame(group=c('static','variable'), vals=c(mean(data1$dist[which(data1$inter_type=='static')]), mean(data1$dist[which(data1$inter_type=='variable')]))) 

p1 <- ggplot() +
	geom_violin(aes(x=data1$inter_type, y=data1$dist), trim=F, scale='width', width=0.9) +
	geom_point(aes(x=point_df$group, y=point_df$vals), shape=13, col='black', size=5) +
        ylim(0, 1200000) +
        labs(x=NULL, y='Distance for significant interactions') +
        my_theme

pdf('Plots/siginter.dist.inter_type.box.pdf')
plot(p1)
dev.off()


# inter class
data1 <- read.table('res_fithic/siginter.20kb.merged.anno.txt', header=T)
data1$inter_class <- factor(data1$inter_class, levels=c('peak_prom', 'peak_none', 'none_none'))

point_df <- data.frame(group=c('peak_prom', 'peak_none', 'none_none'), vals=c(mean(data1$dist[which(data1$inter_class=='peak_prom')]), mean(data1$dist[which(data1$inter_class=='peak_none')]), mean(data1$dist[which(data1$inter_class=='none_none')]))) 

p1 <- ggplot() +
	geom_violin(aes(x=data1$inter_class, y=data1$dist), trim=F, scale='width', width=0.9) +
	geom_point(aes(x=point_df$group, y=point_df$vals), shape=13, col='black', size=5) +
        ylim(0, 1200000) +
        labs(x=NULL, y='Distance for significant interactions') +
        my_theme

pdf('Plots/siginter.dist.inter_class.box.pdf')
plot(p1)
dev.off()


##================================================================================================
## rna correlation by number of promoter-centered siginter

data1 <- read.table('anno_res/sig_inter_count.per_gene.txt', header=T)

##
hesc_0 <- data1$log2tpm[which(data1$sample=='hESC' & data1$prom_siginter==0)]
hesc_1 <- data1$log2tpm[which(data1$sample=='hESC' & data1$prom_siginter>=1 & data1$prom_siginter<4)]
hesc_4 <- data1$log2tpm[which(data1$sample=='hESC' & data1$prom_siginter>=4 & data1$prom_siginter<7)]
hesc_7 <- data1$log2tpm[which(data1$sample=='hESC' & data1$prom_siginter>=7 & data1$prom_siginter<10)]
hesc_10 <- data1$log2tpm[which(data1$sample=='hESC' & data1$prom_siginter>=10)]

me_0 <- data1$log2tpm[which(data1$sample=='ME' & data1$prom_siginter==0)]
me_1 <- data1$log2tpm[which(data1$sample=='ME' & data1$prom_siginter>=1 & data1$prom_siginter<4)]
me_4 <- data1$log2tpm[which(data1$sample=='ME' & data1$prom_siginter>=4 & data1$prom_siginter<7)]
me_7 <- data1$log2tpm[which(data1$sample=='ME' & data1$prom_siginter>=7 & data1$prom_siginter<10)]
me_10 <- data1$log2tpm[which(data1$sample=='ME' & data1$prom_siginter>=10)]

ec4_0 <- data1$log2tpm[which(data1$sample=='EC4' & data1$prom_siginter==0)]
ec4_1 <- data1$log2tpm[which(data1$sample=='EC4' & data1$prom_siginter>=1 & data1$prom_siginter<4)]
ec4_4 <- data1$log2tpm[which(data1$sample=='EC4' & data1$prom_siginter>=4 & data1$prom_siginter<7)]
ec4_7 <- data1$log2tpm[which(data1$sample=='EC4' & data1$prom_siginter>=7 & data1$prom_siginter<10)]
ec4_10 <- data1$log2tpm[which(data1$sample=='EC4' & data1$prom_siginter>=10)]

ec12_0 <- data1$log2tpm[which(data1$sample=='EC12' & data1$prom_siginter==0)]
ec12_1 <- data1$log2tpm[which(data1$sample=='EC12' & data1$prom_siginter>=1 & data1$prom_siginter<4)]
ec12_4 <- data1$log2tpm[which(data1$sample=='EC12' & data1$prom_siginter>=4 & data1$prom_siginter<7)]
ec12_7 <- data1$log2tpm[which(data1$sample=='EC12' & data1$prom_siginter>=7 & data1$prom_siginter<10)]
ec12_10 <- data1$log2tpm[which(data1$sample=='EC12' & data1$prom_siginter>=10)]

ec24_0 <- data1$log2tpm[which(data1$sample=='EC24' & data1$prom_siginter==0)]
ec24_1 <- data1$log2tpm[which(data1$sample=='EC24' & data1$prom_siginter>=1 & data1$prom_siginter<4)]
ec24_4 <- data1$log2tpm[which(data1$sample=='EC24' & data1$prom_siginter>=4 & data1$prom_siginter<7)]
ec24_7 <- data1$log2tpm[which(data1$sample=='EC24' & data1$prom_siginter>=7 & data1$prom_siginter<10)]
ec24_10 <- data1$log2tpm[which(data1$sample=='EC24' & data1$prom_siginter>=10)]

ec48_0 <- data1$log2tpm[which(data1$sample=='EC48' & data1$prom_siginter==0)]
ec48_1 <- data1$log2tpm[which(data1$sample=='EC48' & data1$prom_siginter>=1 & data1$prom_siginter<4)]
ec48_4 <- data1$log2tpm[which(data1$sample=='EC48' & data1$prom_siginter>=4 & data1$prom_siginter<7)]
ec48_7 <- data1$log2tpm[which(data1$sample=='EC48' & data1$prom_siginter>=7 & data1$prom_siginter<10)]
ec48_10 <- data1$log2tpm[which(data1$sample=='EC48' & data1$prom_siginter>=10)]

nonec12_0 <- data1$log2tpm[which(data1$sample=='nonEC12' & data1$prom_siginter==0)]
nonec12_1 <- data1$log2tpm[which(data1$sample=='nonEC12' & data1$prom_siginter>=1 & data1$prom_siginter<4)]
nonec12_4 <- data1$log2tpm[which(data1$sample=='nonEC12' & data1$prom_siginter>=4 & data1$prom_siginter<7)]
nonec12_7 <- data1$log2tpm[which(data1$sample=='nonEC12' & data1$prom_siginter>=7 & data1$prom_siginter<10)]
nonec12_10 <- data1$log2tpm[which(data1$sample=='nonEC12' & data1$prom_siginter>=10)]

nonec48_0 <- data1$log2tpm[which(data1$sample=='nonEC48' & data1$prom_siginter==0)]
nonec48_1 <- data1$log2tpm[which(data1$sample=='nonEC48' & data1$prom_siginter>=1 & data1$prom_siginter<4)]
nonec48_4 <- data1$log2tpm[which(data1$sample=='nonEC48' & data1$prom_siginter>=4 & data1$prom_siginter<7)]
nonec48_7 <- data1$log2tpm[which(data1$sample=='nonEC48' & data1$prom_siginter>=7 & data1$prom_siginter<10)]
nonec48_10 <- data1$log2tpm[which(data1$sample=='nonEC48' & data1$prom_siginter>=10)]

##
pdf('Plots/siginter.log2tpm_rna.box.pdf', height=4, width=12, pointsize=3)
boxplot(
hesc_0,hesc_1,hesc_4,hesc_7,hesc_10,
me_0,me_1,me_4,me_7,me_10,
ec4_0,ec4_1,ec4_4,ec4_7,ec4_10,
ec12_0,ec12_1,ec12_4,ec12_7,ec12_10,
ec24_0,ec24_1,ec24_4,ec24_7,ec24_10,
ec48_0,ec48_1,ec48_4,ec48_7,ec48_10,
nonec12_0,nonec12_1,nonec12_4,nonec12_7,nonec12_10,
nonec48_0,nonec48_1,nonec48_4,nonec48_7,nonec48_10,
names=c(
'hesc_0','hesc_1','hesc_4','hesc_7','hesc_10',
'me_0','me_1','me_4','me_7','me_10',
'ec4_0','ec4_1','ec4_4','ec4_7','ec4_10',
'ec12_0','ec12_1','ec12_4','ec12_7','ec12_10',
'ec24_0','ec24_1','ec24_4','ec24_7','ec24_10',
'ec48_0','ec48_1','ec48_4','ec48_7','ec48_10',
'nonec12_0','nonec12_1','nonec12_4','nonec12_7','nonec12_10',
'nonec48_0','nonec48_1','nonec48_4','nonec48_7','nonec48_10'
),
outline=F,notch=T, las=2,
xlab='Number of P-cRE interactions', ylab='RNA log2TPM'
)
dev.off()






