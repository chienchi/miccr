## MInimap2 Contig ClassifieR (MICCR)

```
usage: ./miccr.py [-h] (-i [FASTA] | -f [PAF]) [-d [FASTA/MMI]] [-dp [PATH]]
                  [-x {asm5,asm10,map-pb,map-ont}] [-t <INT>] [--stdout]
                  [-o [DIR]] [-p <STR>] [-mlp <FLOAT>] [--silent] [-v]

MInimap2 Contig ClassifieR (MICCR) 0.0.1

optional arguments:
  -h, --help            show this help message and exit
  -i [FASTA], --input [FASTA]
                        Input one or multiple contig files in FASTA format.
                        Use space to separate multiple input files.
  -f [PAF], --paf [PAF]
                        Input a Minimap2 PAF file.
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
  --stdout              Output to STDOUT [default: False]
  -o [DIR], --outdir [DIR]
                        Output directory [default: .]
  -p <STR>, --prefix <STR>
                        Prefix of the output file [default:
                        <INPUT_FILE_PREFIX>]
  -mlp <FLOAT>, --minLcaProp <FLOAT>
                        LCA classified segments more than a proportion of
                        contig length [default: 0.1]
  --silent              Disable all messages.
  -v, --verbose         Provide verbose running messages.
```
