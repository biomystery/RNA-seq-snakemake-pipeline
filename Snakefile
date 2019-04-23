
import json
from os.path import join, basename, dirname

# globals ----------------------------

configfile: 'config.yml'
# Full path to an uncompressed FASTA file with all chromosome sequences.
CDNA = config['CDNA']

# Full path to a folder where intermediate output files will be created.
OUT_DIR = config['OUT_DIR']

FILES = json.load(open(config['SAMPLES_JSON']))


SAMPLES = sorted(FILES.keys())


# Functions -------------------------------------------------------------------
def rstrip(text, suffix):
        # Remove a suffix from a string.
        if not text.endswith(suffix):
                    return text
                    return text[:len(text)-len(suffix)]


# Rules ------------------------------

rule all:
	input:
		join(dirname(CDNA), 'kallisto', basename(CDNA).rstrip(".fa")),
		[OUT_DIR + "/" + x for x in expand('{sample}/abundance.tsv', sample = SAMPLES)],
		'abundance.tsv.gz'

rule kallisto_index:
	input:
		cdna = CDNA
	output:
		index = join(dirname(CDNA), 'kallisto', rstrip(basename(CDNA), '.fa'))
	log:
		 'logs/kallisto_index.log'
	shell:
		"""	
		kallisto index  --index={output.index} --make-unique {input.cdna}  >> {log} 2>&1
		"""

rule kallisto_quant:
	input:
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2'],
		index = rules.kallisto_index.output.index
	threads: 4
	output:
		join(OUT_DIR, '{sample}', 'abundance.tsv'),
		join(OUT_DIR, '{sample}', 'run_info.json')
	log:
	    'logs/{sample}_kallistos_quant.log'
	run:
            
	    shell('kallisto quant -t {threads}'+' -i {input.index} {input.r1} {input.r2}'+
                  ' --output-dir='+join(OUT_DIR,'{wildcards.sample}')+' >> {log} 2>&1')


rule collate_kallisto:
    input:
        expand(join(OUT_DIR, '{sample}', 'abundance.tsv'), sample= SAMPLES)
    output:
        'abundance.tsv.gz'
    benchmark:
        "benchmarks/collate.txt"
    run:
        import gzip

        b = lambda x: bytes(x, 'UTF8')

        # Create the output file.
        with gzip.open(output[0], 'wb') as out:

            # Print the header.
            header = open(input[0]).readline()
            out.write(b('sample\t' + header))

            for i in input:
                sample = basename(dirname(i))
                lines = open(i)
                # Skip the header in each file.
                lines.readline()
                for line in lines:
                        fields = line.strip().split('\t')
                        if float(fields[4]) > 0:
                                out.write(b(sample + '\t' + line))
