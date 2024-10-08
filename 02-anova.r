

library(stringr)
input="results/rpkm_table.tsv";
output="results/l2fc_table";

tt=read.table(input,header=T);
meta=as.data.frame(do.call("rbind", strsplit(colnames(tt)[-1],"_")))
colnames(meta)=c("s","g","t")

dd=NULL;
for( i in unique(meta$s)){
	d0=tt[,-1][,meta$s %in% i & meta$t == "0h" ] + 0.1;
	d1=tt[,-1][,meta$s %in% i & meta$t != "0h" ] + 0.1;
	x=log2(d1/d0)
	if(is.null(dd)){ dd=x;
	}else{ dd=cbind(dd,x); }
}
row.names(dd)=tt[,1]
dd.meta = meta[ meta$t != "0h",]
dd.meta$g = factor(dd.meta$g,levels=c("StableCOPD","AECOPD"))



library(lme4)
library(lmerTest)  # For p-values
library(emmeans)
perform_anova <- function(x) {
  anova_data <- data.frame(
    expression = unlist(x),
    subject = factor(dd.meta$s),
    group = factor(dd.meta$g,levels=c("StableCOPD","AECOPD")),
    time = factor(dd.meta$t)
  )
 
  fit <- tryCatch({
	lmer(expression ~ group * time + (1 | subject), data = anova_data)
  }, warning = function(w) return(NULL),
     error = function(e) return(NULL))
  ret= data.frame(gt.pval=1,g.pval=1,t.pval=1);
  if(!is.null(fit)){
 	ar=anova(fit)
  	ret$gt.pval= ar["group:time","Pr(>F)"]
  	ret$g.pval= ar["group","Pr(>F)"]
  	ret$t.pval= ar["time","Pr(>F)"]
  }
  return(ret)
  #emm <- emmeans(fit, ~ group | time)
  #pairwise_comparisons <- pairs(emm, adjust = "tukey")  # Adjust method can be changed as needed
}
r=apply(dd, 1, perform_anova)
d <- do.call(rbind,x) lapply(r, function(x) data.frame(pval = x$p_val)))
write.table(cbind(dd,d), file=paste0(output,"_l2fc_over_0h_anova_group_time_pval.tsv"),col.names=T,row.names=F,quote=F,sep="\t"))

### group-diff
i=d$g.pval < 0.01
m=as.matrix(dd[i,])
o=order(dd.meta$t)
library(ComplexHeatmap)
pdf(file=paste0(output,"_group_pval_001_full_heatmap.pdf"),height=14);
Heatmap(t(scale(t(m[,o]))), column_split=dd.meta$g[o], cluster_columns=F,show_row_names=T,row_names_gp = gpar(fontsize =6))
dev.off();

i=d$t.pval < 0.01
m=as.matrix(dd[i,])
o=order(dd.meta$t)
library(ComplexHeatmap)
#pdf(file=paste0(output,"_group_time_pval_001_full_heatmap.pdf"),height=14);
Heatmap(t(scale(t(m[,o]))), column_split=dd.meta$t[o], cluster_columns=F,show_row_names=T,row_names_gp = gpar(fontsize =6));
#dev.off();


## average patient expression
library(dplyr)
x=paste(dd.meta$g,dd.meta$t)
m1=data.frame(row.names=row.names(m))
for( i in x ){
	m1[[i]] = apply(m[,x %in% i],1,mean);
}

m1= as.matrix(m1)
g1= str_extract(colnames(m1),"\\D+COPD")
row_clusters <- cutree(hclust(dist(m1)), k = 8)  # Define 5 clusters for rows
top_annotation <- HeatmapAnnotation(bar = anno_barplot(column_means))

pdf(file=paste0(output,"_pval_001_short_heatmap.pdf"),height=21,width=5);
Heatmap(m1, column_split=g1, cluster_columns=F,show_row_names=F,top_annotation=top_annotation, row_names_gp=gpar(fontsize=6))
dev.off();

library(ComplexHeatmap)

# Assuming m1 is your matrix and g1 defines 5 clusters for rows
# Calculate row means for each cluster
m1=t(scale(t(m1)))
row_means <- rowMeans(m1)

# Create a barplot annotation for the row clusters
left_annotation <- rowAnnotation(bar = anno_barplot(row_means[row_clusters]))
top_annotation <- HeatmapAnnotation(bar = anno_barplot(colMeans(m1)))

# Create heatmap with row split and the barplot annotation on the left
Heatmap(m1, 
        row_split = row_clusters, 
        cluster_rows = TRUE, 
        cluster_columns = FALSE, 
        show_row_names = FALSE, 
        top_annotation = top_annotation,
        left_annotation = left_annotation)




m1_df <- as.data.frame(m1)
m1_df$Gene <- rownames(m1_df)  # Add gene names to the data frame
m1_df$Cluster <- factor(row_clusters)  # Add the cluster information

# Convert the matrix to long format for ggplot
m1_long <- m1_df %>%
  tidyr::pivot_longer(-c(Gene, Cluster), names_to = "Sample", values_to = "Expression")

# Extract the group and time information from the sample names (assuming they are embedded)
m1_long$Group <- str_extract(m1_long$Sample, "\\D+COPD")  # Extract group information (e.g., StableCOPD, AECOPD)
m1_long$Time <- str_extract(m1_long$Sample, "\\d+h")  # Extract time points (e.g., 0h, 2h, 4h, 6h)
m1_long$Time <- factor(m1_long$Time, levels = c("0h", "2h", "4h", "6h"))  # Ensure time is treated as a factor

# Create line and dot plots for each gene cluster
ggplot(m1_long, aes(x = Time, y = Expression, color = Group, group = interaction(Gene, Group))) +
  geom_line() +  # Connect the points with lines for each gene and group
  geom_point(size = 2) +  # Draw dots for each time point
  facet_wrap(~ Cluster, scales = "free_y") +  # Create a separate plot for each gene cluster
  theme_minimal() +  # Use a minimal theme
  labs(title = "Gene Expression Patterns per Cluster",
       x = "Time Points", y = "Expression Level") +
  scale_color_manual(values = c("blue", "red"))  # Adjust color palette for groups


