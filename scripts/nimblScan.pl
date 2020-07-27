#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use File::Basename;
use File::Path;



##################################################################################################################################
######## Perl wrapper to automate feature extraction of images of Nimblegen arrays.
###################################################################################################################################
#
# Based on a script originally written by Keith Ching.
# http://consultching.com/root/?p=8
# Adapted for zircon and features added/removed/fixed by Ludo Pagie.
# Further extended by Wouter Meuleman, in particular to include 
# Ringo/LOESS normalization and clean-up of code/procedure.
#
# The input is a (series of) multipage tif images from the Agilent
# scanner at the CMF. The images are supposed to have the red signal
# in the first image, green in the 2nd. Also, the images are expected
# to be oritented wrongly. Using using a tool called 'flipRotate' they
# are flipped and rotated, the images are then split. The two
# resulting images are named as the original with the file base
# extended with _Cy5/_Cy3.
#
# Images are subsequently analyzed with NimbleScan-2.5.
# Nimblescan makes an auto-alignment, with the 'local-alignment'
# optimization option turned off by default. 
# For every image it produces three data files:
#
# .pair file: raw data file, not including all nimblegen controls (use
# quantification data file to get all controls)
#
# .grd file: small file with info regarding alignment info
#
# _hyb.data: QC data file with all sorts of info regarding
# hybridization quality, also includes CHIP_ID and DESIGN_ID
# encodings, most precise 'QC score'
#
# Command line options (can be abbreviated):
# -images:        all images to be analyzed
# -designfile:    full path to designfile corresponding to the image files
# -output:        directory where to place all generated images and data files
# -localalign:    (optional) apply optimized local alignment of grid
# -dyeswap:       (optional) hybridisations are performed in dye-swap
# -verbose:       (optional) print some progess lines


#
# Modified by Carolyn 4/7/11
# Will now take manual alignments. 
#


###################################################################################################################################
######## INITIALIZATION
###################################################################################################################################
my $TIFFSPLIT = '/usr/bin/tiffsplit';
my $FLIPROTATE = '/usr/local/bin/flipRotate';
my $NIMBLESCAN = '/usr/local/bin/NimbleScan-2.5';
my ($ch1_col,$ch2_col) = ('Cy5','Cy3');
my %nimblescanOptions = ('outformat'  => [' -x probes',         # .pair (probesctrl doesn't give all controls anyway ..)
                                          '-x quantification',  # _raw.data (feature report, plus all nimblegen ctrls)
                                          '-x qcmetrics',       # _hyb.data; real qc data
#                                         'macro',          # .macro; sparse data file
#                                         'autoqc',    # _hyb.data (qcmetrics)
#                                         'uniformity', #.unif
#                                         'complexqc', #.cqc (????)
#                                         'dose', # crashes
#                                         'cqcprod' # doesn't work??
                                         ],
                                         'outdir'     => '-o ',
#                                        'skipalign'  => '-a',
                                         'localalign' => '-l',
                                         'imagefiles' => '-i "',
                                         'designfile' => '-n '
                        );
my $localalign = 0;
my $dyeswap = 0;
my $VERBOSE = 0;

###################################################################################################################################
######## COMMAND LINE PROCESSING
###################################################################################################################################

my $usage = "usage: nimblescanWrap.pl -images inputfiles -designfile designfile -output outputdir [-localalign] [-dyeswap] [-verbose]\n";

## Options:
my (@images, $designfile, $output) = ((),(),());
GetOptions('images=s{,}'     => \@images,
           'designfile=s'    => \$designfile,
           'output=s'        => \$output,
           'localalign'      => \$localalign,
           'dyeswap'         => \$dyeswap,
           'verbose'         => \$VERBOSE);
(@images and defined($designfile) and defined($output)) or die $usage;

## Check options
if(scalar(@images)<1) {
  die "'images' argument incorrect";
}
foreach my $im (@images) {
  die "cannot access image $im" if(! -e $im);
}
if(! -d $output) {
  $VERBOSE and print STDERR "==VERBOSE== Outputdir $output does not exist, creating\n";
  mkpath($output) or die "can't create directory $output";
}
$output .= '/' unless $output =~ /.+\/$/;
$nimblescanOptions{'outdir'} .= $output;

die 'no, or incorrect, designfile: $designfile' unless defined($designfile) and -e $designfile;
$nimblescanOptions{'designfile'} .= $designfile;

###################################################################################################################################
######## TIFF FILE PROCESSING (INCLUDING NIMBLESCAN)
###################################################################################################################################
### TIFFs are split up into two channels (red/green) and
### transposed, ie mirrored vertically and rotated 270 degrees. 
### Resulting images are named as follows: BASENAME_Cy5/Cy3.tif

foreach my $im (@images) {
  die "$im not a TIFF image (ie no '.tif' extension)"
    unless $im =~ m/(?:.*\/)*(.*)\.tif/;
  my $imagebasename = $1;  # base name used for all generated files

  # If the pair files already exist, go the next image.
  if(-f $output.$imagebasename . "_" . $ch1_col . "_fr.pair" and 
     -f $output.$imagebasename . "_" . $ch2_col . "_fr.pair") { 
    print "PAIR files exist, continuing.\n";
    next; 
  }

  # Split image
  print "Splitting $im\n";

  my $command = "$TIFFSPLIT $im $output$imagebasename";
  $VERBOSE and print "==VERBOSE== " , $command, "\n";
  my $logfile = $output . $imagebasename . ".log";
  `$command > $logfile 2>&1`;

  my $ch1Image = $output . $imagebasename . "_" . $ch1_col . ".tif";
  rename $output . $imagebasename . 'aaa.tif', $ch1Image;
  my $ch2Image = $output . $imagebasename . "_" . $ch2_col . ".tif";
  rename $output . $imagebasename . 'aab.tif', $ch2Image;

  unlink $output . $imagebasename . 'aac.tif';
  unlink $output . $imagebasename . 'aad.tif';

  # Flip and rotate images
  print "Rotating and flipping $im\n";

  my $ch1Image_orig = $ch1Image;
  $ch1Image =~ s/\.tif/_fr\.tif/g;
  $command = "$FLIPROTATE $ch1Image_orig $ch1Image";
  $VERBOSE and print "==VERBOSE== ", $command, "\n";
  `$command >> $logfile 2>&1`;
  unlink $ch1Image_orig;

  my $ch2Image_orig = $ch2Image;
  $ch2Image =~ s/\.tif/_fr\.tif/g;
  $command = "$FLIPROTATE $ch2Image_orig $ch2Image";
  $VERBOSE and print "==VERBOSE== ", $command, "\n";
  `$command >> $logfile 2>&1`;
  unlink $ch2Image_orig;


#skips this part, since need to do a manual alignment

#  # Perform grid alignment using NimbleScan 
#  $nimblescanOptions{'imagefiles'} = "-i $ch1Image";
#  $command = join(' ',($NIMBLESCAN, (@{$nimblescanOptions{'outformat'}}), 
#       @nimblescanOptions{'outdir','designfile','imagefiles'}));
#  $command .= ' ' . $nimblescanOptions{'localalign'} if $localalign;
#  print "Aligning grid for $ch1Image\n";
#  $VERBOSE and print "==VERBOSE== ", $command,"\n";
#  `$command >> $logfile 2>&1`;
#  #unlink $ch1Image;

#  $nimblescanOptions{'imagefiles'} = "-i $ch2Image";
#  $command = join(' ',($NIMBLESCAN, @{$nimblescanOptions{'outformat'}}, 
#       @nimblescanOptions{'outdir','designfile','imagefiles'}));
#  $command .= ' '.$nimblescanOptions{'localalign'} if $localalign;
#  print "Aligning grid for $ch2Image\n";
#  $VERBOSE and print "==VERBOSE== ", $command,"\n";
#  `$command >> $logfile 2>&1`;
  
 #unlink $ch2Image;

#$manually fillling this variable. 
$localalign = 1;

  ## check whether the auto-align option of NimbleScan made a mess
  if($localalign == 1) {
    print "Checking whether the local-alignment went OK\n";
    my $Rcommand = <<EOR;
args <- commandArgs(trailingOnly=TRUE)
if(is.na(args))
  stop("no filename argument given on command line")
fname <- as.character(args[1])
dat <- try(read.delim(fname,skip=1))
str(dat)
if(inherits(dat,"try-error")){
  stop(paste("can't read data file:",fname))
}
xdiff <- unlist(tapply(dat[['X_PIXEL']],dat[['X']],function(x)(diff(x))))
ydiff <- unlist(tapply(dat[['Y_PIXEL']],dat[['Y']],function(x)(diff(x))))
if(min(xdiff) < -3 | max(xdiff) > 3 | min(ydiff) < -3 | max(ydiff) > 3) {
  fname <- sub('_raw.data','.error',fname)
  write.table(x=xdiff,file=fname,quote=FALSE,sep='\t',row.names=FALSE)
  write.table(x=ydiff,file=fname,quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE,append=TRUE)
}
q('no')
EOR

    foreach my $col (($ch1_col, $ch2_col)) {
      my $newbasename = $output.$imagebasename."_".$col."fr";
      my $datafile = $newbasename . "_raw.data";

      ## run R with script and filename as argument
      open Rsession, "|R --vanilla --args $datafile >> $logfile 2>&1" or die "can't start R process: $!";
      print Rsession $Rcommand;
      close Rsession;

      ## check whether R discovered an error
      if(-e $newbasename.".error"){
        print STDERR "\nlocal align option of NimbleScan may have messed up your data in: \n$imagebasename\n";
        print STDERR "You should consider running without option '-localalign'\n\n";
      }
    }
  }
}


###################################################################################################################################
######## PAIR FILE PROCESSING 
###################################################################################################################################
### Pair files are loaded and (Ringo) LOESS normalized
foreach my $im (@images) {
#  $im =~ m/(.*\/*.*)\.tif/;
  $im =~ m/(?:.*\/)*(.*)\.tif/;
  my $imagebasename = $1;  # base name used for all generated files
  my $logfile = $output . $imagebasename . ".log";
  my $RDatafile = $output . $imagebasename . ".RData";

  print "Ringo/LOESS normalizing $im\n";

  $designfile =~ m/(.*\/*.*)\.ndf/;
  my $posfile = $1.".pos";
  my $spottypefile = $1."_spottypes.txt";
  my $targetfile = $output . $imagebasename . "_targets.txt";

  open(TARGETFILE, "> $targetfile");
  print(TARGETFILE "SlideNumber\tSpecies\tFileNameCy3\tFileNameCy5\tCy3\tCy5\n");
  print(TARGETFILE $imagebasename . "\tspecies\t", $output.$imagebasename . "_Cy3_fr.pair\t" . 
                   $output.$imagebasename . "_Cy5_fr.pair\t");
  if ($dyeswap) {
    print(TARGETFILE "Experiment\tControl\n");
  } else {
    print(TARGETFILE "Control\tExperiment\n");
  }
  close(TARGETFILE);

  my $Rcommand = <<EOR;
args <- commandArgs(trailingOnly=TRUE)
if(is.na(args[1]))
  stop("no filename argument given on command line")
if(is.na(args[2]))
  stop("no annotation filename argument given on command line")
if(is.na(args[3]))
  stop("no dyeswap-indication given on command line")

posfile       <- as.character(args[1]);
targetfile    <- as.character(args[2]);
spottypefile  <- as.character(args[3]);
RDatafile     <- as.character(args[4]);
dyeswap       <- args[5];

library(Ringo)
RG <- readNimblegen(targetfile, spottypefile, path=NULL)
loess <- preprocess(RG, method="loess", returnMAList=TRUE)

pos <- read.delim(posfile);

loess_dat <- data.frame(PROBE_ID=loess\$genes\$PROBE_ID, M=loess\$M, A=loess\$A)
colnames(loess_dat) <- c("PROBE_ID", "M", "A");

if (dyeswap == 1) loess_dat\$M <- -1 * loess_dat\$M;
merged <- merge(loess_dat, pos, by="PROBE_ID", sort=FALSE)

dat <- data.frame(seqname = merged\$CHROMOSOME,                 start = merged\$POSITION, 
                  end     = merged\$POSITION+merged\$LENGTH-1,  score = merged\$M, 
                  A       = merged\$A,                          probe = merged\$PROBE_ID);

ord <- order(dat\$seqname, dat\$start);
dat <- dat[ord,];

dat\$seqname  <- as.character(dat\$seqname);
dat\$probe    <- as.character(dat\$probe);

save(dat, file=RDatafile);

q('no')
EOR

  ## run R with script and filename as argument
  open Rsession, "|R --vanilla --args $posfile $targetfile $spottypefile $RDatafile $dyeswap >> $logfile 2>&1" 
    or die "can't start R process: $!";
  print Rsession $Rcommand;
  close Rsession;
}

###################################################################################################################################
######## CLEAN-UP
###################################################################################################################################
#system("gzip -f " . $output . "/*.pair");
#system("gzip -f " . $output . "/*_raw.data");


