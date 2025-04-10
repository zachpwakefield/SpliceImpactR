# Data within /extData/rawData was generated through manually generating a 
# limited number of alternative RNA processing events across AFE and SE that
# span events which will be significantly different across condition and events
# which will not be significantly different across conditions. Along with this,
# we generated .exon files to go with the AFE data.
#   
# Having both significant and non-significant evens ensures the rest of the
# function of the package will be able to be minimally used.
# 
# Data within extData/ was extracted by limiting the files generated through
# the getTranscripts getTranslations setupBiomart and setupAnnotation to a 
# defined set of genes -- partially extracted from the genes contained in the 
# rawData files and partially randomly selected.
# 
# As is done in normal usage of SpliceImpactR, these are generated through 
# biomaRt and from GENCODE.
# 
# This is done because the normal versions of these can be large to store and 
# so that the examples run quickly.
# 
# The PPIDM_GoldDDIs.csv.gz is also generated through a limited run of the 
# getTTI setup function for speed and memory efficiency.
