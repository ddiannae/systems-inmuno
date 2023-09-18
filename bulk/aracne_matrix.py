from MultiAracne import Aracne
import sys, os, argparse

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp', '-e', help="tsv file of expression matrix", type= str)
    parser.add_argument('--out', '-o', help="output directory for single MI files", type=str)
    parser.add_argument('--matrix', '-m', help="output file for final mi matrix", type=str, default="matrix.adj")
    parser.add_argument('--cores', '-c', help="number of cores available", type= int, default=2)

    args=parser.parse_args()

    exp_matrix = args.exp
    outdir = args.out
    outmatrix = args.matrix
    procs = int(args.cores)

    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    ma = Aracne(exp_matrix)
    #ma.run(processes=procs, outdir=outdir, pval=1)
    ma.build_triu_missing_genes(outdir, outmatrix)
