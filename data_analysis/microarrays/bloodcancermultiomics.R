library("BloodCancerMultiOmics2017")
library("Biobase")
library("limma")
library("tidyverse")

data("exprTreat", "drugs")
e <- exprTreat
pd <- pData(e)

colnames( pData(e) ) <- sub( "PatientID", "Patient", colnames( pData(e) ) )
colnames( pData(e) ) <- sub( "DrugID", "Drug", colnames( pData(e) ) )
pData(e)$Drug[ is.na(pData(e)$Drug) ] <- "none"
keep_cols <- which( pData(e)$Drug %in% c("none","D_002"))
e <- e[,keep_cols]
pData(e)$Drug <- relevel( factor( pData(e)$Drug ), "none" )
pData(e)$SampleID <- colnames(e)
colnames(e) <- paste( pData(e)$Patient, pData(e)$Drug, sep=":" )
na_rows <- apply(exprs(e), 1, function(x) any(is.na(x))) 
e <- e[!na_rows,]

#e_df <- bind_cols(
#  data.frame(ProbeID = fData(e)$ProbeID, Symbol = fData(e)$Symbol),
#  as.data.frame(exprs(e))
#  )
# write_csv(e_df, "microarray_ibrutinib_individuals.csv")


fData(e) <- fData(e)[ , c( "ProbeID", "Entrez_Gene_ID", "Symbol",
                           "Cytoband", "Definition" ) ]

mm <- model.matrix( ~ 0 + Patient + Drug, pData(e) )
colnames(mm) <- sub( "Patient", "", colnames(mm) )
colnames(mm) <- sub( "Drug", "", colnames(mm) )
head(mm)

dim(mm)

fit <- lmFit( e, mm )
ebfit <- eBayes( fit)

cont_interest <- "D_002"
hist(ebfit$p.value[,cont_interest])


sum(p.adjust(ebfit$p.value[,cont_interest],method="BH") <= 0.05)


ebfit$df.prior # 10.43969
ebfit$s2.prior * fit$stdev.unscaled[1,cont_interest]^2 # 0.007078982

tbl_df <- tibble(id = rownames(fit),
                 limma_pvalue = fit$p.value[,cont_interest],
                 mu_hat = fit$coefficient[,cont_interest],
                 se_hat = fit$sigma * fit$stdev.unscaled[,cont_interest],
                 residual_dof = fit$df.residual)
                 
write_csv(tbl_df, "microarray_ibrutinib.csv")
