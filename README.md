## MInimap2 Contig ClassifieR (MICCR)

MICCR is a taxonomically classifier for contigs that leverage minimap2 to map sequences to a reference genome database.
Based on the mapping results, MICCR annotated the mapped regions of a contig by analyzing LCA taxonomy from the best hits.

### REQUIREMENTS
* python >= 3.0
* pandas >= 0.23.0

### RESULTS

MICCR provides two tsv files `PREFIX.ctg.tsv` and `PREFIX.lca_ctg.tsv`.
The '.ctg.tsv' has all 11 columns below to show all of taxonomic annotation regions by regions.
The '.lca_ctg.tsv' includes first 10 columns below to display the LCA results of qualified regions in each contig.

| COLUMN | NAME         | DESCRIPTION                                                         |
|--------|--------------|---------------------------------------------------------------------|
| 1      | CONTIG       | Name of the contig                                                  |
| 2      | LENGTH       | Length of the contig                                                |
| 3      | START        | Start position of taxonomic annotated region (0-based)              |
| 4      | END          | End position of taxonomic annotated region                          |
| 5      | LCA_TAXID    | LCA taxonomy ID of all qualified alignments mapped to this region   |
| 6      | LCA_RANK     | LCA taxonomy rank of all qualified alignments mapped to this region |
| 7      | LCA_NAME     | LCA taxonomy name of all qualified alignments mapped to this region |
| 8      | HIT_COUNT    | Number of accounted alignments                                      |
| 9      | AGG_LENGTH   | Additional annotated length this mapped region can provide          |
| 10     | AVG_IDENTITY | Approximate mapping identity based on minimap2                      |
| 11     | AGG_REGION   | Start and end positions of additional annotated length              |

### USAGE
```
usage: ./miccr.py [-h] (-i [FASTA] | -f [PAF]) [-d [FASTA/MMI]] [-dp [PATH]]
                  [-x {asm5,asm10,map-pb,map-ont}] [-t <INT>] [-c] [-o [DIR]]
                  [-p <STR>] [-mp <FLOAT>] [-if <FLOAT>] [--silent] [-v]

MInimap2 Contig ClassifieR (MICCR) 0.0.1

optional arguments:
  -h, --help            show this help message and exit
  -i [FASTA], --input [FASTA]
                        Input one or multiple contig files in FASTA format.
                        Use space to separate multiple input files.
  -f [PAF], --paf [PAF]
                        Input a PAF alignment file.
  -d [FASTA/MMI], --database [FASTA/MMI]
                        Name/path of readmapper's index [default: None]
  -dp [PATH], --dbPath [PATH]
                        Path of databases. If dbPath isn't specified but a
                        path is provided in "--database" option, this path of
                        database will also be used in dbPath. Otherwise, the
                        program will search "database/" in program directory.
                        [default: database/]
  -x {asm5,asm10,map-pb,map-ont}, --platform {asm5,asm10,map-pb,map-ont}
                        You can specify one of the following platform:
                        "asm5"    : Long assembly to reference mapping (avg divergence < 5%);
                        "asm10"   : Long assembly to reference mapping (avg divergence < 10%);
                        "asm20"   : Long assembly to reference mapping (avg divergence < 20%);
                        "map-pb"  : PacBio/Oxford Nanopore read to reference mapping;
                        "map-ont" : Slightly more sensitive for Oxford Nanopore to reference mapping;
                        [default: 'asm10']
  -t <INT>, --numthreads <INT>
                        Number of cpus [default: 1]
  -c, --stdout          Output to STDOUT [default: False]
  -o [DIR], --outdir [DIR]
                        Output directory [default: .]
  -p <STR>, --prefix <STR>
                        Prefix of the output file [default:
                        <INPUT_FILE_PREFIX>]
  -mp <FLOAT>, --minLcaProp <FLOAT>
                        Classify contigs by finding LCA of mapped segments >
                        specified proportion of contig length
  -if <FLOAT>, --iqrfactor <FLOAT>
                        Specify a facter (f). Classify mapped segments with
                        agg_len > Q1+f*IQR, where Q1/3=first/third quartile of
                        agg_len and IQR=(Q3-Q1). [default: 1]
  --silent              Disable all messages.
  -v, --verbose         Provide verbose running messages.
```