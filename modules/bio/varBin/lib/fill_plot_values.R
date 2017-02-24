
chr8    3076521 3078128 .       1608    .       1608    22.992  1748    2
24.218  20.269  21.794  33.695
23.435  21.855  22.548  29.266
0       0       0       0

The following plost a pure calculation of modal copy number (used above), followed by that hmm result

d <- read.table(file('stdin'), header=FALSE, sep="\t")

abs <- d[9]
mabs <- median(abs)
cn <- mabs / abs * REF_PLOIDY


mabs/abs[rs]*2   "Strain Copy Number"
points(x,d[rs,10],col='red')

mnbc <- d[8]



The following plots a pure calculation of copy number change, followed by the hmm result
where nbc = normalized bin coverage
bs = bin size (after gap correction)
m = median, etc.
here 2 is being used as REF_PLOIDY
d[,c+16]-mnbc)/mnbc * cn "Copy Number Change"


