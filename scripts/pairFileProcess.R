###################################################################################################################################
######## PAIR FILE PROCESSING 
###################################################################################################################################
### Pair files are loaded and (Ringo) LOESS normalized


#./nimbleScanPair2rData.pl -input /mnt/CMF-resstore/2011-05/13/CMF_503942_S01.tif -designfile ../../arrayLayouts/100718_MM8_Economy_03_HX1.ndf -output ../all.pair.rarDir/CMF_503942_S01_Cy3_fr -verbose -dyeswap

#have to have cy3 and cy5 files in a directory named same as the file name
#eg '../all.pair.rarDir/CMF_503942_S01/CMF_503942_S01_Cy3_fr.pair

posfile       <- "../../arrayLayouts/100718_MM8_Economy_03_HX1.pos";
targetfile    <- "../all.pair.rarDir/CMF_504022_S01/CMF_504022_S01_targets.txt"
spottypefile  <- "../../arrayLayouts/100718_MM8_Economy_03_HX1_spottypes.txt"
RDatafile     <- "../all.pair.rarDir/CMF_504022_S01/CMF_5040222_S01.RData"
dyeswap       <- 1

library(Ringo)
RG <- readNimblegen(targetfile, spottypefile, path=NULL)

#adding this line in as is missing all the probes when readNimblegen runs controlStatus()
#may be ok for different pair files????
RG$genes$Status[RG$genes$GENE_EXPR_OPTION == "EXPERIMENTAL"] <- "Probe"
RG$genes$Status[RG$genes$GENE_EXPR_OPTION == "RANDOM"] <- "Random"


print("preprocess")
loess <- preprocess(RG, method="loess", returnMAList=TRUE)
head(loess)

pos <- read.delim(posfile);
head(pos)

#actual problem is that my probe info files don't match probe names in pair file. 

loess_dat <- data.frame(PROBE_ID=loess$genes$PROBE_ID, M=loess$M, A=loess$A)
colnames(loess_dat) <- c("PROBE_ID", "M", "A");
head(loess_dat)

if (dyeswap == 1) loess_dat$M <- -1 * loess_dat$M;

head(loess_dat)
head(pos)

merged <- merge(loess_dat, pos, by="PROBE_ID", sort=FALSE)

dat <- data.frame(seqname = merged$CHROMOSOME,                 start = merged$POSITION, 
                  end     = merged$POSITION+merged$LENGTH-1,  score = merged$M, 
                  A       = merged$A,                          probe = merged$PROBE_ID);

ord <- order(dat$seqname, dat$start);
dat <- dat[ord,];

dat$seqname  <- as.character(dat$seqname);
dat$probe    <- as.character(dat$probe);

save(dat, file=RDatafile);

q('no')
EOR

  ## run R with script and filename as argument
  print "Ringo/LOESS normalizing $in\n";
  print "|R --vanilla --args $posfile $targetfile $spottypefile $RDatafile $dyeswap >> $logfile 2>&1";
  open Rsession, "|R --vanilla --args $posfile $targetfile $spottypefile $RDatafile $dyeswap >> $logfile 2>&1" 
    or die "can't start R process: $!";
  print Rsession $Rcommand;
  close Rsession;
}
