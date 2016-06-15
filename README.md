Convert VCF to mutect context for trinucleotide analysis[¶](#Convert-VCF-to-mutect-context-for-trinucleotide-analysis) {#Convert-VCF-to-mutect-context-for-trinucleotide-analysis}
======================================================================================================================

#### written by xiuque[¶](#written-by-xiuque)

Usage[¶](#Usage) {#Usage}
----------------

You will need to have bedtools in your path You will need to have Hongs'
Mutect\_context.py script in the same folder

VCFtoMutec.py -i infile -o outfile -g

supply full path is preferred

Dependencies[¶](#Dependencies) {#Dependencies}
------------------------------

In []:

    #!/usr/bin/env python

    ## to get the context sequences
    import pybedtools
    import sys
    import os
    from collections import defaultdict

    # Or whatever you call hong's script
    import pandas as pd
    import itertools
    import argparse

Create an object to hold information[¶](#Create-an-object-to-hold-information) {#Create-an-object-to-hold-information}
------------------------------------------------------------------------------

VariantFileObject is suppose to be a vcf file object that will convert
allow conversion of vcf file to other format and vice versa. Now it only
has vcf -\> bed -\> mutec/fasta

In []:

    class VariantFileObject(object):
        def __init__(self, file_name):
            self._filename = file_name
            self._variant = []
            self._format = ""
            self.__name__ = os.path.basename(file_name)
        
        

### Class method 1 : Read Vcf and store it as attributes in the object[¶](#Class-method-1-:-Read-Vcf-and-store-it-as-attributes-in-the-object) {#Class-method-1-:-Read-Vcf-and-store-it-as-attributes-in-the-object}

Follows the VCF4.1 format.
[link](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)

In []:

    ## read vcf and store the vcf in 2 self._variant and self._position.
        # give each snp an id. position is used to find the sequences
        # should change this function to be more generic, for instance if line starts with ('#') and not ##, make it the column line 
        #append all columns into self._variant/
        
        def readvcf(self):
            with open(self._filename) as f:
                snp_id = 0
                for line in f.readlines():
                    if not line.startswith("#"):
                        line = line.strip()
                        fields = line.split('\t')
                        snp_id = snp_id + 1
                        self._variant.append([fields[0],int(fields[1]),snp_id,fields[3],fields[4]])
                        
            self._variant = pd.DataFrame(self._variant, columns=['chr','pos','snp_id','ref','alt'])
            self._format = "vcf"
            print self._variant
        

#### The VCF is now store as in the VariantFileObject with the attribute : self.\_variant[¶](#The-VCF-is-now-store-as-in-the-VariantFileObject-with-the-attribute-:-self._variant) {#The-VCF-is-now-store-as-in-the-VariantFileObject-with-the-attribute-:-self._variant}

### Class Method 2 : Drop indels[¶](#Class-Method-2-:-Drop-indels) {#Class-Method-2-:-Drop-indels}

Drop indels based on the variant size

In []:

     ## function to drop indels for from either vcf or bed formation, might be sueful in the future 
        def drop_indels(self):
            if self._format == "vcf":
                ## remove anything that has more tahn 2 character in either the ref or the alt
                self._variant = self._variant[self._variant.ref.map(len) == 1]
                self._variant = self._variant[self._variant.alt.map(len) == 1]
            elif self._format == "bed":
                ## not used now, maybe useful in future
                name = self._variant['name']
                self._variant.join(name.apply(lambda x: pd.Series(x.split(':'))))
                self._variant = self._variant[self._variant.iloc[-1].map(len) == 1]
                self._variant = self._variant[self._variant.iloc[-2].map(len) == 1]

### Class Method 3 : Convert VCF to Bed[¶](#Class-Method-3-:-Convert-VCF-to-Bed) {#Class-Method-3-:-Convert-VCF-to-Bed}

The name field in a bed file will have its name in the following format
"snp\_id,ref,alt". will including a -2, +4 to the region in the vcf
file. VCF file is 1-based : having the 1st based at position 1. bed is 1
based for the starting, having the 1st based at position 1 and 0-based
and half open. I should seperate the conversion and the addition of
steps seperately

In []:

        # convert vcf to bed, wiht each of the fasta header containing information about the vcf
        def bedformat(self):
            if self._format != "bed":
                self._variant['start'] = self._variant.pos - 3
                self._variant['end'] = self._variant.pos + 3
                self._variant['name'] = self._variant.snp_id.map(str) + ';' + self._variant.ref + ';' + self._variant.alt
                self._variant = self._variant[['chr','start','end','name']]
                self._format = "bed"
        

### Class Method 4 : Write Bedfiles[¶](#Class-Method-4-:-Write-Bedfiles) {#Class-Method-4-:-Write-Bedfiles}

Writes the bedformat to a temporary file. Should replace this in future
so that i will not ned to write to the disk. But it is easier this way

In []:

    # write bed result to bed file. 
        def bedtofile(self):
            if self._format == "bed":
                bedfile = "%s_tmp.ed" % self._filename
                self.bedfile = bedfile
                self._variant.to_csv(bedfile, sep='\t',index=False, header=None)
                return bedfile

### Class Method 5 : Fasta to mutec[¶](#Class-Method-5-:-Fasta-to-mutec) {#Class-Method-5-:-Fasta-to-mutec}

In []:

        # convert fasta files from bed to muteect
        def to_mutect(self, bedtools):
            rows = []
            columns = ["snp_id","ref_allele","alt_allele","context"]
            self.out_file = "%s_mutect.csv" % self._filename
            fh_out = open(self.out_file, 'w+')
            fh_out.write('Something is brewing\n')
            os.unlink(self.bedfile)
            with open(Bedtools.seqfn) as f:
                for header, sequence in itertools.izip_longest(*[f]*2):
                    header = header.strip()
                    sequence = sequence.strip()
                    sequence = sequence[:3] + "x" + sequence[3:]
                    header = header.split(';')
                    header.append(sequence)
                    rows.append(header)

            df = pd.DataFrame(rows, columns=columns)
            df['judgement'] = 'KEEP'
            df.to_csv(fh_out, sep='\t')
            fh_out.close()

### Run Main Function[¶](#Run-Main-Function) {#Run-Main-Function}

In []:

            
            
    if __name__ == '__main__':
        
        parser = argparse.ArgumentParser()
        parser.add_argument("-i", dest="inFileName", required=True,
                            help="Input Vcf F")
        parser.add_argument("-o", dest="outFileName", required=True,
                            help="the output file")
        parser.add_argument("-v", dest="verbosity", type=int, default=0,
                            help="verbosity level 0-2 [default=0]")
        parser.add_argument("-g", dest="genome", required=True, 
                            help="Genome file")
        args = parser.parse_args()

In []:

        genome_fasta = args.genome
        ## read in vcf file and store it as vcf
        VF = VariantFileObject(args.inFileName)
        VF.readvcf()
        
        ## drop indels
        VF.drop_indels()
        
        ## convert to bed and add nucleotide 
        VF.bedformat()
        
        ## Use bedtools
        Bedtools = pybedtools.BedTool(VF.bedtofile())
        Bedtools = Bedtools.sequence(genome_fasta, name=True)
        VF.to_mutect(Bedtools)
















    __author__ = 'quek'
