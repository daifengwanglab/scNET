#source code : https://davetang.org/muse/2010/11/10/gene-ontology-enrichment-analysis/



library("biomaRt")
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
 
my_chr <- c(1:22, 'M', 'X', 'Y')
my_ensembl_gene <- getBM(attributes='ensembl_gene_id',
                    filters = 'chromosome_name',
                    values = my_chr,
                    mart = ensembl)
 

my_entrez_gene <- getBM(attributes='entrezgene_id',
                    filters = 'chromosome_name',
                    values = my_chr,
                    mart = ensembl)

# get some more info on the entrez_gene
my_attribute <- c('entrezgene_id',
                  'hgnc_symbol',
                  'chromosome_name',
                  'start_position',
                  'end_position',
                  'strand')
 
my_entrez_gene_info  <- getBM(attributes=my_attribute,
                        filters = c('entrezgene_id', 'chromosome_name'),
                        values = list(entrezgene=my_entrez_gene$entrezgene_id, chromosome_name=my_chr),
                        mart = ensembl)
 

my_entrez_gene_info=my_entrez_gene_info[,1:2]
