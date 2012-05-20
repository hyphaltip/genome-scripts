module load mavid
module load stajichlab
module load cndsrc

omap2hmap genomes < pre.map > pre.h.map
makeBreakpointGraph pre.h.map treefile
mkdir bp_alignments
makeBreakpointAlignmentInput --out-dir=bp_alignments
mavidAlignDirs --init-dir=bp_alignments
findBreakpoints pre.h.map treefile edges bp_alignments > breakpoints
breakMap breakpoints < pre.h.map > better.h.map
hmap2omap genomes < better.h.map > better.map
phits2constraints -i ../input < pairwisehits > constraints

mkdir alignments
makeAlignmentInput --map=better.map . alignments
mavidAlignDirs --init-dir=alignments

mkdir fsa_alignments
makeAlignmentInput --map=better.map . fsa_alignments
