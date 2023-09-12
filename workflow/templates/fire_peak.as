table fire_peaks "FIRE Peaks"
(

string  chrom;  "Reference sequence chromosome or scaffold"
uint    chromStart; "Alignment start"
uint    chromEnd;   "Alignment end"
string  name;   "Position of the matching duplication"
uint    score;  "fake score"
char[1] strand; "+ or - or . for unknown"
float   signalValue;    "Measurement of average enrichment for the region",
float   pValue; "Statistical significance of signal value (-log10). Set to -1 if not used."
float   qValue; "Statistical significance with multiple-test correction applied (FDR -log10) Set to -1 if not used."
uint    peak;   "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."
)
