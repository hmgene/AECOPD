
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

i=!is.na(d$t.pval) & d$t.pval < 0.001

m=as.matrix(dd[i,])
o=order(t)
Heatmap(m[,o], column_split=g[o], cluster_columns=F)

m1= split(m, paste(c(g,t)))

.data.frame(

split(m,paste(c(g,t)))




 )))

g1=str_extract(colnames(m1),"(\\D+COPD)",group=T)

Heatmap(split




aecopd_cols <- grep("AECOPD", colnames(m), value = TRUE)
stablecopd_cols <- grep("StableCOPD", colnames(m), value = TRUE)

# Function to normalize each group by the 0h time point
calculate_l2fc <- function(data, group_cols) {
  # Extract 0h time points for the group
  group_0h_cols <- grep("_0h", group_cols, value = TRUE)
  
  # Loop over each 0h sample to calculate L2FC for subsequent time points
  for (col_0h in group_0h_cols) {
    # Extract the prefix to match with other time points (remove "_0h")
    sample_prefix <- sub("_0h", "", col_0h)
    time_cols <- grep(paste0(sample_prefix, "_[246]h"), colnames(data), value = TRUE)
    data[time_cols] <- log2((data[time_cols]+0.1) / (data[[col_0h]]+0.1))
  }
  
  return(data)
}

# Apply normalization for AECOPD and StableCOPD groups
m_l2fc <- m  # Create a copy of the data
m_l2fc[aecopd_cols] <- calculate_l2fc(m[, aecopd_cols], aecopd_cols)
m_l2fc[stablecopd_cols] <- calculate_l2fc(m[, stablecopd_cols], stablecopd_cols)
m_l2fc=m_l2fc[,- grep("_0h",colnames(m_l2fc))];

o=order(str_extract(colnames(m_l2fc),"\\d+h"))
Heatmap(m_l2fc[,o],cluster_columns=F, column_split = str_extract(colnames(m_l2fc)[o],"(\\D+COPD)",group=T))








library(ComplexHeatmap)
Heatmap(t(scale(t(m))),column_split=g)


