### FASTA processing -----

## Load packages ----

library(tidyverse)
library(seqinr)
library(here)

## helper function ----

remove_duplicate_string <- function(x){
   
   teo <- str_split(x, pattern = " ") %>%
      unlist()
   
   teo2 <- paste(unique(teo), collapse = "--")
   
   return(teo2)
}

# >sp|P02489|CRYAA_HUMAN Alpha-crystallin A chain OS=Homo sapiens OX=9606 GN=CRYAA PE=1 SV=2
# this is how the headers should look like

## Load the fasta ----

fasta_file_location <- "data/sample.fasta"

fasta_file_name <- basename(fasta_file_location) %>% 
   str_remove(".fasta$")

# read fasta file
fasta <- read.fasta(file = fasta_file_location, 
                    seqtype = "AA", 
                    as.string = TRUE)

## Checking ----

#head(fasta)
#
#head(names(fasta))
#
#getAnnot(head(fasta))
#
#getName(fasta) %>% head()
#
#tnames <- getName(fasta) 

## Get annotation of headers ----

headers <- getAnnot(fasta) %>% 
      unlist() %>% 
      data.frame(header = .) %>%
      mutate(header = str_remove(header, ">"))

# separate regular headers from variant headers ----
regular_headers <- headers %>% 
   dplyr::filter(str_detect(string = header,
                            pattern = "variant:", 
                            negate = TRUE))

variant_headers <- headers %>% 
   dplyr::filter(str_detect(string = header,
                            pattern = "variant:", 
                            negate = FALSE))

# get organism annotation from the header information  
os_obj <- regular_headers$header[1] %>% 
   str_extract(pattern = "OS=.*") %>%
   str_remove(pattern = " GN=.*") 

## prepare tabular annotation of variant identifiers -----

sep_var_head1 <- variant_headers %>% 
   mutate(gene = str_extract(header, "(?<=GN=).*")) %>%
   mutate(header_mod = str_remove(header, "OS=.*")) %>%
   separate(col = header_mod, 
            into = c("first_part", "description"),
            sep = " ",
            remove = TRUE, 
            extra = "merge")

sep_var_head2 <- sep_var_head1 %>%
   mutate(description_simple = map_chr(description, remove_duplicate_string)) %>%
   #mutate(nvariants = str_count(description_simple, pattern = "variant:")) %>%
   separate_rows(description_simple, sep = "--") %>% # use the pattern used to collapse to separate the different variants into rows
   filter(description_simple != "") %>% 
   separate(col = description_simple,
            sep = "\\|",
            into = c("var_code", "type", "recurrence", "gene",
                     "ensembl_gene", "is_transcript", "ensembl_transcript",
                     "orf_type", "fraction_x", "substitution_x", "substitution_aa",
                     "fraction_y", "fraction_z","fraction_zz", "coord"),
            extra = "merge") %>%
   separate(col = first_part, 
            into = c("prefix","protein_seq_id","transcript_id"),
            sep = "\\|") %>%
   mutate(protein_seq_id = str_remove(protein_seq_id, "^transcript:"),
          transcript_id = str_remove(transcript_id, "^transcript:"),
          prefix = str_replace(prefix, "mz", "sp")) %>%
   mutate(protein_seq_id = str_trunc(protein_seq_id, 70, "right", ellipsis = "X"))

variant_annotation <- sep_var_head2 %>%
   dplyr::select(-prefix)

nasthere <- unique(variant_annotation$type) %>% # check if there are NAs in the type of variant
   is.na() %>% 
   any()

if(nasthere == TRUE){
   inter_vars <- c(unique(na.omit(variant_annotation$type)), "NA")
} else {
   inter_vars <- c(unique(na.omit(variant_annotation$type)))
}



## new headers for variable sequences ----

new_headers_var <- sep_var_head2 %>%
   distinct(header, .keep_all = TRUE) %>%
   dplyr::select(header, prefix, protein_seq_id, transcript_id, type, description, var_code, gene) %>%
   dplyr::mutate(init_desc = paste(transcript_id, "variant ")) %>%
   pivot_wider(id_cols = c(header, prefix, protein_seq_id, transcript_id, description, init_desc, var_code, gene), 
               names_from = type,
               values_from = type,
               values_fn = function(x) x[1]) %>%
   unite(col = "variant_desc1",
         all_of(inter_vars),
         na.rm = TRUE) %>%
   unite(col = "var_description",
         init_desc,variant_desc1,
         sep = " ",
         na.rm = TRUE) %>%
   mutate(identifier = paste(prefix,protein_seq_id,transcript_id, sep = "|"),
          description = var_description,
          rest = paste(
             os_obj,
             "OX=0000",
             paste0("GN=",gene),
             "PE=0 SV=0"
          )) %>%
   mutate(new_header = paste(identifier, description, rest)) %>%
   dplyr::select(header, new_header)

## new headers for regular/canonical sequences -----

new_headers_reg <- regular_headers %>% 
   distinct(header) %>%
   mutate(new_header = str_replace(header, "mz", "sp")) %>%
   mutate(new_header = str_remove_all(new_header, "transcript:"))


# merge tables of new headers (variant and regular)

tab_new_headers <- bind_rows(new_headers_reg,
                             new_headers_var) %>%
   left_join(headers,.) # merge with original headers so they have the same order  

# rename fasta -----

rename_fasta <- fasta

names(rename_fasta) <- tab_new_headers$new_header

## Save new fasta ----  

write.fasta(rename_fasta, names = names(rename_fasta), 
            file.out = here(paste0("results/",fasta_file_name,"_phil_compatib.fasta")))

## Save variant annotation tabular file

write_tsv(variant_annotation,
          file = here(paste0("results/",fasta_file_name,"_variant_annotation.tsv")))


################# TESTS AND EVALUATION #############

#nrow(tab_new_headers)
#nrow(headers)
#
#
#unique(sep_var_head4$type)
#
#table(sep_var_head3$type)
#head(sep_var_head3$first_part)
#sep_var_head4$len_identifier %>% hist()
#
#ggplot(sep_var_head4,
#       aes(x = type, y = len_identifier, fill = type)) + 
#   geom_boxplot() + 
#   theme(axis.text.x = element_text(angle = 90))
#
#sep_headers <- headers %>% 
#   separate(col = header,
#            into = c("initiator", "descrition", "rest"),
#            sep = " ") #%>% 
#


#head(headers$header, 3)

#   separate(col = rest, 
#            into = c("ident_2", "description_1", "description_2", "description_3", "desc_4", "desc_5"),
#            sep = " ")


