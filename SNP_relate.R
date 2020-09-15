#install.packages("SNPrelate")
library(gdsfmt)
library(SNPRelate)

# both plink and VCF files can be used for PCA analysis 

# import plink files 

bed.fn<-"My_data.bed"
fam.fn<-"My_data.fam"
bim.fn<-"My_data.bim"

#-------------------------------------------------------------------#
# import vcf file 

#vcf.fn <- "My_data.vcf"

#snpgdsVCF2GDS(vcf.fn, "My_data.gds", method="biallelic.only")

#-------------------------------------------------------------------#

# Convert to GDS 

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "My_data.gds")

genofile_My_data <- snpgdsOpen("My_data.gds")

#run pca 

pca_My_data  <- snpgdsPCA(genofile_My_data,num.thread=50)

# extracting principal components

pca_My_data_eigvect <- data.frame(PC1=pca_My_data$eigenvect[,1],  # the first eigenvector
                            PC2=pca_My_data$eigenvect[,2], # the second eigenvector
                            PC3=pca_My_data$eigenvect[,3],
                            PC4=pca_My_data$eigenvect[,4],
                            PC5=pca_My_data$eigenvect[,5],
                            PC6=pca_My_data$eigenvect[,6],
                            PC7=pca_My_data$eigenvect[,7],
                            PC8=pca_My_data$eigenvect[,8],
                            PC9=pca_My_data$eigenvect[,9],
                            PC10=pca_My_data$eigenvect[,10],
                            stringsAsFactors = FALSE)


# plot PCs
plot(pca_My_data_eigvect $PC1,pca_My_data_eigvect $PC2)

# write to a table
write.csv(pca_My_data_eigvect ,"pca_My_data_eigvect.csv",row.names = T,quote = F)


# variant components
pca_My_data$varprop

PC_variance_my_data<-data.frame(pca_My_data$varprop,stringsAsFactors = FALSE)

write.csv(PC_variance,"PC_variance_my_data.csv")
