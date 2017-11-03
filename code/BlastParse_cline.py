import sys

def BlastParseExtra(infile, genome_fa, best_hit_criteria, Eval_thresh, get_frags_switch = 1, window_size = 2000, verb = 0):
    """
    Usage: BlastParseExtra.py  <blast_xml_output> <genome_fasta> <best_hit_criteria> <Eval_thresh> <get_frags_switch> <window_size> <verb>
    
    <blast_xml_output>  -  absolute path to blast xml output
    <genome_fasta>      -  absolute path to genome fasta file
    <best_hit_criteria> -  factor by which the best e-value must be lower than the second best (recommended: 1e-5)
    <Eval_thresh>       -  e-value threshold for unique alignments
    <get_frags_switch>  -  switch (1 = on, 0 = off) for getting fragments around the mapping
    <window_size>       -  size of the window (in bp) around the mapping coordinates used to extract the subject scaffold segment 
    <verb>              -  (switch, 0 or 1) to make it chatty

    This script first filters the mappings in the <blast_xml_output> for uniq hits with evalues better than
    <Eval_thresh> or for multi hits where the best hit is <best_hit_criteria> orders of magnitude better than 
    the second.
    
    It will then retrieve a segment of the scaffold from <genome_fasta> which is + and - the <window_size> around
    the mapping coordinates for each query. If the ends of the scaffold are not within this window, then the 
    length of the segment will be (length of mapped query sequence + 2 x <window_size>). However if an end of a
    scaffold is within this window, the segment will be trimmed to this length.
    
    """


    import sys
    from Bio.Blast import NCBIXML
    from Bio import SeqIO
    import gzip

    if infile.endswith("gz"):
	handle = gzip.open(infile, 'r')
    else:
	handle = open(infile, 'r')

    blast = NCBIXML.parse(handle)

    good_blast_outs = {}
    multi_counter = 0
    unique_counter = 0
    ## From Alan's script: Returns blast hits only when the best e-value is 5 orders of magnitude better than the second best.

    for record in blast :
        if len(record.alignments)==1:
            if record.alignments[0].hsps[0].expect <= Eval_thresh:
                unique_counter += 1
                good_blast_outs[record.query] = {}
                good_blast_outs[record.query]["Ref_hit_id"] = str(record.alignments[0].hit_def)
                good_blast_outs[record.query]["Evalue"] = float(record.alignments[0].hsps[0].expect)
                good_blast_outs[record.query]["Hit_start_coord"] = int(record.alignments[0].hsps[0].sbjct_start)
                good_blast_outs[record.query]["Hit_end_coord"] = int(record.alignments[0].hsps[0].sbjct_end)
                #print "Uniq\t%s\t%s\t%s\t%s\t%s" % (record.query, good_blast_outs[record.query]["Ref_hit_id"], good_blast_outs[record.query]["Evalue"], good_blast_outs[record.query]["Hit_start_coord"], good_blast_outs[record.query]["Hit_end_coord"])


        elif len(record.alignments)>1:
            if all([record.alignments[0].hsps[0].expect <= Eval_thresh, record.alignments[0].hsps[0].expect < best_hit_criteria * record.alignments[1].hsps[0].expect]):
                multi_counter += 1
                good_blast_outs[record.query] = {}
                good_blast_outs[record.query]["Ref_hit_id"] = str(record.alignments[0].hit_def)
                good_blast_outs[record.query]["Evalue"] = float(record.alignments[0].hsps[0].expect)
                good_blast_outs[record.query]["Hit_start_coord"] = int(record.alignments[0].hsps[0].sbjct_start)
                good_blast_outs[record.query]["Hit_end_coord"] = int(record.alignments[0].hsps[0].sbjct_end)

    if verb == "1":

        print "Number of multi-alingments kept:", multi_counter
        print "Number of unique alingments kept:", unique_counter


    Rtemp_summary_out = open("%s/blast_%s_summary.out" % (infile.rpartition("/")[0], window_size), 'w')
    Rtemp_summary_out.write("query\tsubject\n")

    if get_frags_switch == "1":
        
        Rtemp_summary_out = open("%s/blast_%s_summary.out" % (infile.rpartition("/")[0], window_size), 'w')
        Rtemp_summary_out.write("query\tsubject\n")

        print "Getting subject scaffold segments from %s . . . " % (genome_fa)

        Rtemp_chunks = open("%s/blast_%s_chunks.fa" % (infile.rpartition("/")[0], window_size), 'w')

        segment_counter = 0

        Rtemp = SeqIO.parse(genome_fa, "fasta")

        recorded_scaffolds = []

        for scaffold in Rtemp:

            for query in good_blast_outs:

                if scaffold.id == good_blast_outs[query]["Ref_hit_id"]: ## If the scaffold has a hit

                    Rtemp_summary_out.write("%s\t%s\n" % (query, scaffold.id))

                    if scaffold.id not in recorded_scaffolds:

                        recorded_scaffolds.append(scaffold.id)

                        if good_blast_outs[query]["Hit_start_coord"] - window_size <= 0 and good_blast_outs[query]["Hit_end_coord"] + window_size >= len(scaffold.seq): # if the beginning and if the end of the scaffold is within the rang of the window

                            SeqIO.write(scaffold, Rtemp_chunks, 'fasta') ## just print whole scaffold
                            segment_counter += 1


                        elif good_blast_outs[query]["Hit_start_coord"] - window_size <= 0 and good_blast_outs[query]["Hit_end_coord"] + window_size < len(scaffold.seq): ## or if the begninning is in range of the window but the end isn't

                            SeqIO.write(scaffold[:good_blast_outs[query]["Hit_end_coord"]+ window_size], Rtemp_chunks, 'fasta') ## print from beginning to upper end of window
                            segment_counter += 1


                        elif good_blast_outs[query]["Hit_start_coord"] - window_size > 0 and good_blast_outs[query]["Hit_end_coord"] + window_size >= len(scaffold.seq): ## or the end of the scaffold is in range but the beginning isnt

                            SeqIO.write(scaffold[good_blast_outs[query]["Hit_start_coord"]- window_size:], Rtemp_chunks, 'fasta') ## print from lower window limit to the end of the scaffold
                            segment_counter += 1

                        elif good_blast_outs[query]["Hit_end_coord"] + window_size < len(scaffold.seq) and good_blast_outs[query]["Hit_start_coord"] - window_size > 0: ## or if neither end of the scaffold is in range of the window
                            SeqIO.write(scaffold[good_blast_outs[query]["Hit_start_coord"]- window_size:good_blast_outs[query]["Hit_end_coord"] + window_size], Rtemp_chunks, 'fasta') ## print from lower window limit to the end of the scaffold
                            segment_counter += 1
        
        print "%s sequence scaffold segments are in %s/blast_%s_chunks.fa" % (segment_counter, infile.rpartition("/")[0], window_size)
        
        Rtemp_chunks.close()
    
    Rtemp_summary_out.close()
    
#print good_blast_outs


### Cline compatibility

if len(sys.argv) < 8:
    sys.exit("\nNot enough arguments\n%s" % BlastParseExtra.__doc__)
else:
    print "running"
    input_file = sys.argv[1]
    genome_fasta = sys.argv[2]
    best_hit_eval_diff = float(sys.argv[3])
    eval_threshold = float(sys.argv[4])
    get_frags = sys.argv[5]
    window = int(sys.argv[6])
    verbose = sys.argv[7]

BlastParseExtra(input_file, genome_fasta, best_hit_eval_diff, eval_threshold, get_frags, window, verbose)

