################################################################################
# Title: Known mutation fasta extraction
# Author: Maria del Carmen Martos Contreras mcar.martos@ipb.csic.es
# Date: 27.10.2026
################################################################################
library(Biostrings)
library(argparse)
library(stringr)

parser <- ArgumentParser()

parser$add_argument(
  "-r", "--reference",
  help = "Reference proteins fasta",
  required = TRUE
)

parser$add_argument(
  "-a", "--assembly",
  help = "Assembly proteins fasta",
  required = TRUE
)

parser$add_argument(
  "-m", "--mutations",
  help = "Known mutations table",
  required = TRUE
)

parser$add_argument(
  "-o", "--output",
  help = "Output file",
  required = TRUE
)

args <- parser$parse_args()

print(args$reference)
print(args$assembly)
print(args$mutations)
print(args$output)

#Process known mutation file
library(dplyr)
library(tidyr)
library(stringr)

mut_wide <- read.csv(args$mutations, header=T)


known_mutations <- mut_wide %>%
  separate_longer_delim(mutation, delim="; ") %>%
  mutate(
    RefAA = str_extract(mutation, "^[A-Z]"),
    Position = as.integer(str_extract(mutation, "[0-9]+")),
    MutAA = str_extract(mutation, "[A-Z/]+$")
  )

known_mutations$Gene<-paste0(known_mutations$Gene, ".1-p1")

ref <- readAAStringSet(args$reference)
asm <- readAAStringSet(args$assembly)

ref_seq <- setNames(as.character(ref), names(ref))
asm_seq <- setNames(as.character(asm), names(asm))
names(ref_seq) <- sub("\\s*\\|.*$", "", names(ref_seq))

#Obtain all mutations
mutations <- data.frame()

for(gene in intersect(names(ref_seq), names(asm_seq))){
  
  r <- strsplit(ref_seq[gene], "")[[1]]
  a <- strsplit(asm_seq[gene], "")[[1]]
  
  idx <- which(r != a)
  print(idx)
  if(length(idx)>0){
    
    mutations <- rbind(
      mutations,
      data.frame(
        Gene = gene,
        Mutation = paste0(r[idx], idx, a[idx])
      )
    )
  }
}
if (nrow(mutations) > 0){
  print(mutations)
  write.table(mutations,paste0(args$output,".all_mutations.tsv"), sep="\t", row.names=FALSE)
  
}

#Obtain know mutation
known_mutations$ObservedAA <- mapply(
  function(gene, pos){
    
    seq <- asm_seq[gene]
    
    if(length(seq)==0) return(NA)
    
    substr(seq, pos, pos)
    
  },
  known_mutations$Gene,
  known_mutations$Position
)

known_mutations$Status <- ifelse(
  known_mutations$ObservedAA == known_mutations$MutAA,
  "Resistant_mutation",
  ifelse(
    known_mutations$ObservedAA == known_mutations$RefAA,
    "Wild_type",
    "Other_variant"
  )
)

print(known_mutations)

#Microsatellite information

microsatellite_results <- data.frame(
  Gene = character(),
  Motif = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

gene <- "PF3D7_1303500.1-p1"

if(gene %in% names(asm_seq)){
  
  seq <- asm_seq[[gene]]
  
  microsatellite_results <- rbind(
    microsatellite_results,
    data.frame(
      Gene = gene,
      Motif = "DNNND",
      Count = str_count(seq, "DNNND")
    )
  )
  
  microsatellite_results <- rbind(
    microsatellite_results,
    data.frame(
      Gene = gene,
      Motif = "DDNHNDNHNND",
      Count = str_count(seq, "DDNHNDNHNND")
    )
  )
  print(microsatellite_results)
  write.table(microsatellite_results,paste0(args$output, ".microsatellite.tsv"),sep="\t", row.names=FALSE, quote=FALSE)
}

# result <- dplyr::left_join(
#   mutations,
#   known_mutations,
#   by = "Gene"
# )
write.table(known_mutations,paste0(args$output,".known_mutations.tsv"), sep="\t", row.names=FALSE)
