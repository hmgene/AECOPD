
out="results/out"

library(stringr)

tt=read.table("out_table.tsv",header=T);
tt_sig=read.table("out.tukey_results.tsv",header=T);
g = str_extract(colnames(tt),"(\\D+COPD)",group=T)[-1]
t = str_extract(colnames(tt),"_(\\d+)",group=T)[-1]
s= str_extract(colnames(tt),"^(\\w+\\d+)\\D+_",group=T)[-1]

dd=NULL;
for( i in unique(s)){
	d0=tt[,-1][,s %in% i & t == "0" ] + 0.1;
	d1=tt[,-1][,s %in% i & t != "0" ] + 0.1;
	x=log2(d1/d0)
	if(is.null(dd)){	
		dd=x;
	}else{
		dd=cbind(dd,x);
	}
}

g = str_extract(colnames(dd),"(\\D+COPD)",group=T)
t = str_extract(colnames(dd),"_(\\d+)",group=T)
s= str_extract(colnames(dd),"^(\\w+\\d+)\\D+_",group=T)

# Function to perform ANOVA for each gene
perform_anova <- function(x) {
  anova_data <- data.frame(
    expression = x,
    group = factor(g),
    time = factor(t)
  )
  aov_results <- aov(expression ~ time + group, data = anova_data)
  aov_summary <- summary(aov_results)
  time_pvalue <- aov_summary[[1]]["time", "Pr(>F)"]
  group_pvalue <- aov_summary[[1]]["group", "Pr(>F)"]
  return( list(time.pval=time_pvalue,group.pval=group_pvalue))
}
r=apply(dd, 1, perform_anova)
d <- do.call(rbind, lapply(r, function(x) data.frame(g.pval = x$group.pval, t.pval=x$time.pval)))
x=data.frame(gene=tt$gene, cbind(dd,d))
write.table(x, file=paste0(out,"_l2fc_over_0h_anova_group_time_pval.tsv"),col.names=T,row.names=F,quote=F,sep="\t"))

i=!is.na(d$t.pval) & d$t.pval < 0.01
#i=!is.na(d$g.pval) & d$g.pval < 0.01
#i=d$t.pval < 0.01 & d$g.pval < 0.01

m=as.matrix(dd[i,])
o=order(t)
library(ComplexHeatmap)
pdf(file=paste0(out,"_pval_001_full_heatmap.pdf"),height=21);
Heatmap(m[,o], column_split=g[o], cluster_columns=F,show_row_names=F)
dev.off();

library(dplyr)

x=paste(g,t)
m1=data.frame(row.names=row.names(m))
for( i in x ){
	m1[[i]] = apply(m[,x %in% i],1,mean);
}

m1= as.matrix(m1)
g1= str_extract(colnames(m1),"\\D+COPD")
pdf(file=paste0(out,"_pval_001_short_heatmap.pdf"),height=21,width=5);
Heatmap(m1, column_split=g1, cluster_columns=F,show_row_names=F)
dev.off();



