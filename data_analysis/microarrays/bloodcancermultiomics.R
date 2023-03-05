library("BloodCancerMultiOmics2017")
library("Biobase")

data("exprTreat", "drugs")
e <- exprTreat

keep_cols <- which( pData(e)$Drug %in% c("none","D_002"))
e <- e[,keep_cols]
colnames( pData(e) ) <- sub( "PatientID", "Patient", colnames( pData(e) ) )
colnames( pData(e) ) <- sub( "DrugID", "Drug", colnames( pData(e) ) )
pData(e)$Drug[ is.na(pData(e)$Drug) ] <- "none"
pData(e)$Drug <- relevel( factor( pData(e)$Drug ), "none" )
pData(e)$SampleID <- colnames(e)
colnames(e) <- paste( pData(e)$Patient, pData(e)$Drug, sep=":" )

head( pData(e) )

fData(e) <- fData(e)[ , c( "ProbeID", "Entrez_Gene_ID", "Symbol",
                           "Cytoband", "Definition" ) ]

mm <- model.matrix( ~ 0 + Patient + Drug, pData(e) )
colnames(mm) <- sub( "Patient", "", colnames(mm) )
colnames(mm) <- sub( "Drug", "", colnames(mm) )
head(mm)

fit <- lmFit( e, mm )
ebfit <- eBayes( fit)

cont_interest <- "D_002"
hist(ebfit$p.value[,cont_interest])


sum(p.adjust(ebfit$p.value[,cont_interest],method="BH") <= 0.05)


ebfit$df.prior # 10.42875
ebfit$s2.prior * fit$stdev.unscaled[1,cont_interest]^2 # 0.007079249

library(tidyverse)
tbl_df <- tibble(id = rownames(fit),
                 limma_pvalue = fit$p.value[,cont_interest],
                 mu_hat = fit$coefficient[,cont_interest],
                 se_hat = fit$sigma * fit$stdev.unscaled[,cont_interest],
                 residual_dof = fit$df.residual)
write_csv(tbl_df, "microarray_ibrutinib.csv")
