
posfile       <- "../../arrayLayouts/100718_MM8_Economy_03_HX1.pos";
targetfile    <- "../all.pair.rarDir/CMF_504022_S01/CMF_504022_S01_targets.txt"
spottypefile  <- "../../arrayLayouts/100718_MM8_Economy_03_HX1_spottypes.txt"
RDatafile     <- "../all.pair.rarDir/CMF_504022_S01/CMF_5040222_S01.RData"
dyeswap       <- 1

library(Ringo)
RG <- readNimblegen(targetfile, spottypefile, path=NULL)
loess <- preprocess(RG, method="loess", returnMAList=TRUE)
