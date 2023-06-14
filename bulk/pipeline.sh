GDC_CLIENT=./
COUNTSDIR=./
PLOTSDIR=./

mkdir $COUNTSDIR/cases $COUNTSDIR/controls
./$GDC_CLIENT download --manifest bulkrna_manifest_cases.2023-06-09.txt --dir $COUNTSDIR/cases --log-file $COUNTSDIR/cases/download.log
./$GDC_CLIENT download --manifest bulkrna_manifest_controls.2023-06-09.txt --dir $COUNTSDIR/controls --log-file $COUNTSDIR/controls/download.log

Rscript getCountsMatrix.R $COUNTSDIR cases
Rscript getCountsMatrix.R $COUNTSDIR controls

Rscript ./getNormalization.R $COUNTSDIR $PLOTSDIR
