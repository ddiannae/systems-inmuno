GDC_CLIENT=./
COUNTSDIR=./
PLOTSDIR=./
ARACNEDIR=./
ARACNESIF=./
ANNOT=./


mkdir $COUNTSDIR/cases $COUNTSDIR/controls
./$GDC_CLIENT download --manifest bulkrna_manifest_cases.2023-06-09.txt --dir $COUNTSDIR/cases --log-file $COUNTSDIR/cases/download.log
./$GDC_CLIENT download --manifest bulkrna_manifest_controls.2023-06-09.txt --dir $COUNTSDIR/controls --log-file $COUNTSDIR/controls/download.log

Rscript getCountsMatrix.R $COUNTSDIR cases $ANNOT
Rscript getCountsMatrix.R $COUNTSDIR controls $ANNOT

Rscript ./getNormalization.R $COUNTSDIR $PLOTSDIR
mkdir $ARACNEDIR
Rscript ./getAracneMatrices.R $COUNTSDIR $ARACNEDIR
singularity run --bind $ARACNEDIR $ARACNESIF aracne_matrix.py -e $ARACNEDIR/vst_controls.tsv -o $ARACNEDIR/vst_controls -m $ARACNEDIR/vst_controls.adj -c 78
singularity run  --bind $ARACNEDIR $ARACNESIF aracne_matrix.py -e $ARACNEDIR/vst_cases.tsv -o $ARACNEDIR/vst_cases -m $ARACNEDIR/vst_cases.adj -c 78
singularity run --bind $ARACNEDIR $ARACNESIF aracne_matrix.py -e $ARACNEDIR/normalized_controls.tsv -o $ARACNEDIR/normalized_controls -m $ARACNEDIR/normalized_controls.adj -c 78
singularity run --bind $ARACNEDIR $ARACNESIF aracne_matrix.py -e $ARACNEDIR/normalized_cases.tsv -o $ARACNEDIR/normalized_cases -m $ARACNEDIR/normalized_cases.adj -c 78


Rscript getNetwork.R $ARACNEDIR/vst_controls.adj 100000 $ANNOT $ARACNEDIR/vst_controls_network_100000
Rscript getNetwork.R $ARACNEDIR/vst_cases.adj 100000 $ANNOT $ARACNEDIR/vst_cases_network_100000
Rscript getNetwork.R $ARACNEDIR/normalized_cases.adj 100000 $ANNOT $ARACNEDIR/normalized_cases_network_100000
Rscript getNetwork.R $ARACNEDIR/normalized_controls.adj 100000 $ANNOT $ARACNEDIR/normalized_controls_network_100000

Rscript getCommunities.R $ARACNEDIR/normalized_controls_network_10000_vertices.tsv $ARACNEDIR/normalized_controls_network_10000_interactions.tsv $ARACNEDIR/normalized_controls_network_10000
Rscript getCommunities.R $ARACNEDIR/normalized_cases_network_10000_vertices.tsv $ARACNEDIR/normalized_cases_network_10000_interactions.tsv $ARACNEDIR/normalized_cases_network_10000
Rscript getCommunities.R $ARACNEDIR/normalized_controls_network_100000_vertices.tsv $ARACNEDIR/normalized_controls_network_100000_interactions.tsv $ARACNEDIR/normalized_controls_network_100000
Rscript getCommunities.R $ARACNEDIR/normalized_cases_network_100000_vertices.tsv $ARACNEDIR/normalized_cases_network_100000_interactions.tsv $ARACNEDIR/normalized_cases_network_100000

Rscript getEnrichments.R $ARACNEDIR/normalized_cases_network_10000_membership.tsv $ARACNEDIR/normalized_cases.tsv $ARACNEDIR/normalized_cases_10000_enrichments.tsv
Rscript getEnrichments.R $ARACNEDIR/normalized_controls_network_10000_membership.tsv $ARACNEDIR/normalized_controls.tsv $ARACNEDIR/normalized_controls_10000_enrichments.tsv
Rscript getEnrichments.R $ARACNEDIR/normalized_cases_network_100000_membership.tsv $ARACNEDIR/normalized_cases.tsv $ARACNEDIR/normalized_cases_100000_enrichments.tsv
Rscript getEnrichments.R $ARACNEDIR/normalized_controls_network_100000_membership.tsv $ARACNEDIR/normalized_controls.tsv $ARACNEDIR/normalized_controls_100000_enrichments.tsv
