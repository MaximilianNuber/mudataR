# Format constants (MuData on-disk RFC, https://mudata.readthedocs.io/stable/io/spec.html)

.mudata_encoding_type <- "MuData"
.mudata_encoding_version <- "0.1.0"
.mudata_format_version <- "0.1.0"
.h5mu_userblock_size <- 512L

# Global slots stored with AnnData element encodings at the store root
.mudata_df_slots <- c("obs", "var")
.mudata_mapping_slots <- c("obsm", "varm", "obsp", "varp", "uns")
.mudata_map_slots <- c("obsmap", "varmap")

.mudata_axes <- c(-1L, 0L, 1L)

mudatar_version <- function() {
  as.character(utils::packageVersion("mudataR"))
}
