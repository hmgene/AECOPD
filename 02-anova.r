

library(stringr)
library(data.table)
input="results/rpkm_table.tsv.gz";

tt <- fread(input, header = TRUE)  # Read the file
tt <- as.data.frame(tt)            # Convert to a data frame
rownames(tt) <- tt$gene            # Set row names to "gene"
tt$gene <- NULL                    # Remove the "gene" column if not needed
tt$V1 = NULL

j=colnames(tt)
meta=as.data.frame(do.call("rbind", strsplit(j,"_")))
colnames(meta)=c("s","g","t")


i=names(sort(apply(tt,1,mean,na.rm=T),decreasing=T)[1:10000])
x=tt[i,]

i = names(sort(apply(x,1,var,na.rm=T),decreasing=T)[1:1000])

#i = names(sort(apply(x,1,var,na.rm=T),decreasing=F)[1:1000])

x=t(scale(t(x[i,])))
Heatmap(x, column_split=meta$g,cluster_columns=F)



dd=NULL;
for( i in unique(meta$s)){
	d0=tt[,j][,meta$s %in% i & meta$t == "0h" ] + 0.1;
	d1=tt[,j][,meta$s %in% i & meta$t != "0h" ] + 0.1;
	x=log2(d1/d0)
	if(is.null(dd)){ dd=x;
	}else{ dd=cbind(dd,x); }
}
row.names(dd)=tt$gene
meta = meta[ meta$t != "0h",]
meta$g = factor(meta$g,levels=c("StableCOPD","AECOPD"))


library(lme4)
library(lmerTest)  # For p-values
library(emmeans)
perform_anova <- function(x) {
  anova_data <- data.frame(
    expression = unlist(x),
    subject = factor(meta$s),
    group = factor(meta$g,levels=c("StableCOPD","AECOPD")),
    time = factor(meta$t)
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
d =cbind(dd, do.call(rbind,r))

output="results/l2fc_test";
write.table(d, file= gzfile(paste0(output,".tsv.gz")),col.names=T,row.names=T,quote=F,sep="\t")
write.table(meta, file=gzfile(paste0(output,"_meta.tsv.gz")),col.names=T,row.names=T,quote=F,sep="\t")

### group-diff
m=as.matrix(d[,!grep("pval",colnames(d))]);
m = as.matrix(d[,colnames(d)[ -grep("pval",colnames(d)) ]])
row.names(row.names(d))

j=order(meta$t)
i=d$gt.pval < 0.001
library(ComplexHeatmap)
Heatmap(t(scale(t(m[i,j]))), column_split=meta$t[j], cluster_columns=F,show_row_names=T,row_names_gp = gpar(fontsize =6))

## average patient expression
library(dplyr)
x=paste(meta$g,meta$t)
m1=data.frame(row.names=row.names(m))
for( i in x ){
	m1[[i]] = apply(m[,x %in% i],1,mean);
}
p=d$gt.pval;
p=p[p<1];
q=p.adjust(p,method="fdr")
pthre=max(p[ q < 0.1 ])

i=d$gt.pval < 0.01

m2= t(scale(t(as.matrix(m1[i,]))))
g1 = str_extract(colnames(m2),"\\D+COPD")
g1 = factor(g1,levels=c("StableCOPD","AECOPD"))
row_clusters <- cutree(hclust(dist(m2)), k = 4)  # Define 5 clusters for rows

cluster_means <- sapply(unique(row_clusters), function(clust) {
  colMeans(m2[row_clusters == clust, , drop = FALSE])
})

colnames(cluster_means)=c("1","2","3","4")
c_colors <- c("1" = "blue", "2" = "green", "3" = "orange", "4" = "purple")
top_annotation <- HeatmapAnnotation(
 	cluster = anno_lines(cluster_means, gp = gpar(col = c_colors[colnames(cluster_means)]))
)
right_annotation <- rowAnnotation(
  cluster = anno_block(gp = gpar(fill = c_colors[as.character(unique(row_clusters))]),
                       labels = as.character(unique(row_clusters)), 
                       labels_gp = gpar(col = "black", fontsize = 10))
)
output="average_l2fc_gt_lt_0.01"

cat(paste(row.names(m2)[ row_clusters == "1" ],collapse="\n"))

# Plot the heatmap with both top annotation and right annotation
pdf(paste0(output,"_heatmap.pdf"),width=4,height=12)
Heatmap(m2,
        column_split = g1, cluster_columns = FALSE, show_row_names = TRUE, row_split = row_clusters, row_names_gp = gpar(fontsize = 2),
        top_annotation = top_annotation, right_annotation = right_annotation)
dev.off();



