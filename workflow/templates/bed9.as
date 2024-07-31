table bed9
"Bed9"
(
string chrom;      "Chromosome (or contig, scaffold, etc.)"
uint   chromStart; "Start position in chromosome"
uint   chromEnd;   "End position in chromosome"
string name;       "Name of item"
uint   score;      "Score from 0-1000"
char[1] strand;    "+ or -"
uint thickStart;   "Start of where display should be thick (start codon)"
uint thickEnd;     "End of where display should be thick (stop codon)"
uint color;        "Primary RGB color for the decoration"
)