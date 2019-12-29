import pandas as pd
import os
import allel as allel
import tabix as tabix
import collections
import pickle

# import cancer gene census
census = pd.read_csv("/Users/imanthonny/Desktop/mutation/cancer_gene_census.csv",delimiter=",")
census = census[~census["Genome Location"].str.contains(":-")]
census["Genome Location"] = "chr" + census["Genome Location"]

# import coding regions
coding = pd.read_csv("/Users/imanthonny/Desktop/mutation/ccdsGene.txt", delimiter="\t", header=None)
coding.rename(columns={coding.columns[6]:'cdsStart', coding.columns[7]:'cdsEnd'}, inplace=True)
coding[2] = coding[2].replace("chr","", regex=True)

# MAF
MAF = pd.read_csv("/Users/imanthonny/Desktop/mutation/MinorAlleleFreq.txt", delimiter="\t")
MAF = MAF[MAF["MAF"]>0.01]

# loop through each VCF to identify if census genome location present
filtered_vcf = pd.DataFrame()
for root, dirs, files in os.walk("/Users/imanthonny/Desktop/mutation/CCLE_VCFs"):
    for file in files:
        if file.endswith("vcf"):
            print(file)
            path = root + "/" + file
            for i in census["Genome Location"]:
                print(i)
                vcf = allel.vcf_to_dataframe(path, region=i)
                if isinstance(vcf,collections.Sized):  # checks if vcf is NOT blank when imported (has census[Genome Location]region)
                        vcf["Sample"] = file.replace(".vcf", "")
                        vcf["Gene"] = census["Gene Symbol"][census["Genome Location"].str.contains(i)].values[0]
                        # if vcf["FILTER_PASS"].any():  # if FILTER_PASS has any True values. 'if not vcf["FILTER_PASS"].any()' would be used for if FILTER pass has any False values.
                        filtered_vcf = filtered_vcf.append(vcf, ignore_index=True)

filtered_vcf = filtered_vcf.drop(columns=["ALT_2","ALT_3", "QUAL"])

# check for coding genes in filtered vcf's
filtered_vcf_2 = pd.DataFrame()
for v in range(0, len(filtered_vcf["POS"])):
    df3 = coding[(coding["cdsStart"] <= filtered_vcf["POS"][v]) & (coding["cdsEnd"] >= filtered_vcf["POS"][v])]
    if isinstance(df3,collections.Sized):
        chr = filtered_vcf["CHROM"][v].replace("chr","")
        df3.rename(columns={df3.columns[2]:'CHR'}, inplace=True)
        # df3["CHR"] = df3["CHR"].astype(str) # isinstance(df3["CHR"], str), checks if object is a string
        matched_chr = df3["CHR"][df3["CHR"].str.contains(chr)]
        # vcf_final = filtered_vcf.iloc[[v]]
        if len(matched_chr) !=  0:
            filtered_vcf.loc[v, "Coding"] = "YES"
        if len(matched_chr) ==  0:
            filtered_vcf.loc[v,"Coding"] = "NO"
        filtered_vcf_2 = filtered_vcf_2.append(filtered_vcf.iloc[[v]], ignore_index=True)


# isolate coding genes only
vcf_coding_only = filtered_vcf_2[filtered_vcf_2["Coding"]=="YES"]

# check for alleleic frequency (common variants only)

vcf_coding_only["CHROM"] = vcf_coding_only["CHROM"].replace("chr","", regex=True)
combined_coordinate = vcf_coding_only["CHROM"] + ":" + vcf_coding_only["POS"].astype(str)
vcf_coding_only["coordinate"] = combined_coordinate
vcf_coding_only = vcf_coding_only[vcf_coding_only["coordinate"].isin(MAF["Location"])]

final_variants = vcf_coding_only
final_variants["Common"] = "YES"

with open('/Users/imanthonny/Desktop/mutation/final_variants.obj', 'wb') as f:
    pickle.dump(final_variants, f)
