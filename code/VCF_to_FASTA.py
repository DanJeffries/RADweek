import sys

def VCF_to_FASTA(vcf_path, cat_path):
    """
    USAGE: python VCF_to_FASTA.py <full/vcf/path> <full/catalog/path>
    
    """
    
    import gzip
    
    outpath = vcf_path.rpartition(".")[0]

    myvcf = open(vcf_path, 'r')

    if cat_path.endswith("gz"):
        my_catalog = gzip.open(cat_path, 'r')
    else:
        my_catalog = open(cat_path, 'r')

    fasta_out = open("%s.fa" % outpath, 'w')

    to_keep = []
    
    print "\nGetting tag numbers"

    tag_count = 0
    for i in myvcf.readlines():
        if not i.startswith("#"):
            to_keep.append(i.split()[2])
            tag_count += 1
    print "\nThere are %s tags in your VCF" % tag_count

    print "\nGetting sequences. . . gimme a minute . . ."

    for line in my_catalog.readlines():
        if not line.startswith("#"):
            if line.split()[2] in to_keep:
                fasta_out.write(">%s\n%s\n" % (line.split()[2], line.split()[8]))
    
    fasta_out.close()
    
    print "\nAll done, tag sequences are in %s.fa\n" % vcf_path.rpartition(".")[0]

if len(sys.argv) < 3:
    sys.exit("Not enough arguments\n%s" % VCF_to_FASTA.__doc__)
else:
    myvcf_path = sys.argv[1]
    mycat_path = sys.argv[2]   

VCF_to_FASTA(myvcf_path, mycat_path) 
