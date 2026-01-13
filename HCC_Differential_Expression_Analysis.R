################################################################################
# Project: Differential Gene Expression Analysis (HCC vs Normal)
# Author: Farah
# Dataset: GSE112790
################################################################################

##--- 1. تحميل المكتبات الأساسية ---
library(GEOquery)
library(limma)
library(hgu133plus2.db)
library(ggplot2)
library(ggrepel)
library(pheatmap)

#--- 2. تحميل البيانات وتصنيف المجموعات ---
gse = getGEO( filename = "GSE112790_series_matrix.txt.gz" ,getGPL = T)
metadata <- pData(gse)
expr <- exprs(gse)

# تعريف المجموعات بناءً على وصف العينات
group  <- ifelse(grepl("normal" , metadata$source_name_ch1 , ignore.case = T  ), "Normal" , "Tumor")
group <- factor(group)

#--- 3. التحليل الإحصائي (Differential Expression Analysis) ---

design = model.matrix(~0 + group )
colnames(design) = levels(group)
# تحديد المقارنة المطلوبة (Tumor vs Normal)
contrast = makeContrasts(Tumor_vs_Normal = Tumor - Normal ,levels = design)
fit = lmFit(expr,design)
fit2 =contrasts.fit(fit , contrast)
fit2 = eBayes(fit2)

# استخراج الجينات ذات الدلالة الإحصائية (Fold Change > 2 and P-value < 0.05)
sigGenes = topTable(fit2 , number = Inf , p.value = 0.05 , lfc = 1)

#--- 4. ربط النتائج بأسماء الجينات (Annotation) ---

my_probes= rownames(sigGenes)
geneMap = select(hgu133plus2.db , keys  =my_probes , columns = "SYMBOL" , keytype = "PROBEID"  )
# دمج الأرقام مع الأسماء في جدول واحد
final_results = merge(sigGenes ,geneMap  , by.x ="row.names" , by.y = "PROBEID" )
colnames(final_results)[1]= "ProbeID"
final_results=final_results[order(final_results$adj.P.Val) ,]

#--- 5. رسم الـ Volcano Plot ---
# تحديد الألوان للرسمة

final_results$diffexpressed ="NO"
final_results$diffexpressed[final_results$logFC > 1 & final_results$adj.P.Val< .05 ] = "UP"
final_results$diffexpressed[final_results$logFC < -1 & final_results$adj.P.Val< .05 ] = "DOWN"

ggplot( data= final_results , aes(x = logFC , y = -log10(adj.P.Val) , col =diffexpressed ))+
  geom_point(alpha =.4 , size = 1.5)+
  theme_minimal()+
  scale_color_manual( values =c("DOWN" = "blue" ,"UP" = "red" ,"NO" = "grey" ))+
  geom_vline(xintercept = c(-1 , +1)  , color = "black", linetype = "dashed" )+
  geom_hline(yintercept = -log10(.05), color = "black", linetype = "dashed" )+
  geom_text_repel(data =head(final_results ,10) , aes (label = SYMBOL) ,size = 4 )+
  labs(title="Volcano Plot: Liver Cancer Analysis", x="log2 Fold Change", y="-log10 Adj P-value")

#--- 6. رسم الـ Heatmap (أهم 50 جين) ---
# تنظيف البيانات من الـ NA والتكرار
final_filtred = final_results[!is.na(final_results$SYMBOL) , ]
final_filtred = final_filtred[!duplicated(final_filtred$SYMBOL) , ]

top_50_ids = final_filtred$ProbeID[1:50]
plot_matrix_clean  = expr[top_50_ids ,]
rownames(plot_matrix_clean ) = final_filtred$SYMBOL[1:50]
# عمل الـ Scaling لجعل البيانات قابلة للمقارنة بصرياً
plot_matrix_clean = t(scale(t(plot_matrix_clean)))

# تجهيز معلومات العينات والألوان للرسمة
sample_info = data.frame(Group = group)
rownames(sample_info) = colnames(expr)
ann_color = list(Group= c(Normal = "blue" , Tumor = "red"))

pheatmap(plot_matrix_clean, 
         annotation_col = sample_info,
         annotation_colors = ann_color,
         main = "Heatmap: Top 50 Differentially Expressed Genes",
         color = colorRampPalette(c("navy" ,"white" ,"firebrick3"))(50),
         show_colnames = F,
         fontsize_row = 8)

# حفظ النتائج النهائية في ملف CSV
rownames(final_results) <- NULL
write.csv(final_results ,"HCC_Significant_Genes.csv")
