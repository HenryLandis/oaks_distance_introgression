#! /usr/bin/awk -f 

BEGIN {
    FS = "\t"
    OFS = FS
}

#print the header info as-is
/^##/ {
    print
    next
}

#add sample named "Reference" to the list of samples
/^#CHROM/ {
    print $0"\tReference"
    next
}

#add homozygous reference allele to every locus.
{
    print $0"\t0/0:20,0:20:PASS:99:0,0,255"
}