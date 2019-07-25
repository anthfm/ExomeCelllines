library(seqminer)
library(dplyr)
library(tidyverse)

census <- read.csv(file="cancer_gene_census.csv", stringsAsFactors=FALSE)
	  census <- census[!grepl(":-", census$Genome.Location),] #get rid of ranges with no values (":-")
	  coding <- read.csv(file = "ccdsGene.txt",header = FALSE, sep = "\t",stringsAsFactors = FALSE)
	  colnames(coding)[which(names(coding) == "V7")] <- "cdsStart"
	  colnames(coding)[which(names(coding) == "V8")] <- "cdsEnd"
	  coding$V3 = gsub('chr', '', coding$V3)
	  allele <- read.csv(file = "MinorAlleleFreq.txt",header = TRUE, sep = "\t",stringsAsFactors = FALSE)
	  allele <- allele[allele$MAF>0.01,]
	  
	  variants <- list()
	  a<-0
	  
	  files <- list.files(path="~/Desktop/vcf", pattern=".vcf.gz$", full.names=TRUE, recursive=FALSE)
	  
	  for (file in files) {
	    df2 <- data.frame(CHROM=integer(0), POS=integer(0), ID=character(0), REF=character(0), ALT=character(0), QUAL=character(0), FILTER=character(0), GENE=character(0))
	    a <- a + 1
	    tt <-  gsub("\\..*$", "", basename(file))
	    print(file)
	    
	    #CANCER GENE CENUS MATCHING + PASS
	    for (i in census$Genome.Location) {
	      sample <- tabix.read.table(file, i)
	      sample <- sample[sample$FILTER == "PASS",]
	      if(nrow(sample) != 0){
	        
	        position <- i
	        genefound <- census$Gene.Symbol[census$Genome.Location == position]
	        sample$GENE <- genefound
	        df2 <- rbind(df2, sample)
	        
	      }
	    }
	    
	    df2 <- select(df2, CHROM, POS, ID, REF, ALT, FILTER, GENE)
	    df2$CODING <- NA
	    
	    # CHECK FOR CODING GENE's
	    #df2 <- df2 %>%
	    #mutate(CODING = map_chr(
	    #.x = POS,
	    #.f = ~ if_else(
	    #condition = any(.x >= coding$cdsStart & .x <= coding$cdsEnd),
	    #true = "YES",
	    #false = NA_character_
	    #)
	    #))
	    
	    df2$CODING <- ifelse(sapply(seq_along(df2$POS), function(i) {
	      inds <- coding$cdsStart <= df2$POS[i] & coding$cdsEnd >= df2$POS[i]
	      any(inds) & (df2$CHROM[i] == coding$V3[which.max(inds)])
	    }), "YES", NA)
	    
	    df2 <- df2[!is.na(df2$CODING),]
	    
	    
	    
	    #CHECK ALLELIC FREQUENCY (COMMON VARIANTS)
	    
	    df3 <- df2 %>% mutate(
	      temp = paste(CHROM, POS, sep = ":")
	    )
	    
	    df3$COMMON <- NA
	    
	    df3 <- df3 %>%
	      mutate(COMMON= map_chr(
	        .x = temp,
	        .f = ~ if_else(
	          condition = any(.x == allele$Location),
	          true = "YES",
	          false = NA_character_
	        )
	      ))
	    
	    df3$temp <- NULL 
	    df3 <- df3[!is.na(df3$COMMON),]
	    
	    a <- a - 1 
	    variants <- append(variants, list(df3), a)
	    names(variants)[(a+1)] <- tt
	    
	    
	    
	  } #end of analysis
	  
