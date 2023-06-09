GDC_CLIENT=../../Apps/gdc-client
OUTDIR=/datos/ot/diana/systems-inmuno/counts

mkdir $OUTDIR/cases $OUTDIR/controls
./$GDC_CLIENT download --manifest bulkrna_manifest_cases.2023-06-09.txt --dir $OUTDIR/cases --log-file $OUTDIR/cases/download.log
./$GDC_CLIENT download --manifest bulkrna_manifest_controls.2023-06-09.txt --dir $OUTDIR/controls --log-file $OUTDIR/controls/download.log

