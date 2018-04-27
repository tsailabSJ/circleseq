######################################################################################################quote
### test_regions_BED.R: make bed file with regions including 
###                     on-target site, 2 off-target sites without variants,
###                     2 off-target sites with variants, and 1 region without off-targets.
######################################################################################################
bed = data.frame(chr=c('2', '8', '1', '12', '4'), start=c(73160981, 120587494, 234492858, 73504668, 48639390), end=c(73161159, 120587517, 234492881, 73504691, 48639413), name=c('2', '8', '1', '12', '4'))

write.table(bed, 'CIRCLEseq_test.bed', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')

### Information about the sites
on_target="2:73160981-73161004"
off_target01="8:120587494-120587517"
off_target02="1:234492858-234492881"
off_target_with_variantWindowOnly="12:73504668-73504691"
off_target_with_variants="4:48639390-48639413"
hotspots="2:73161104-73161159"
######################################################################################################
######################################################################################################
