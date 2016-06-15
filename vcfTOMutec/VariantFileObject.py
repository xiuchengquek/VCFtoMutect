
import argparse
import itertools

class Variant:
    def __init__(self, chr, pos, id, ref, alt):
        self.chr = chr
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        if len(ref) == 1 and len(alt) == 1:
            self.snp = True
        else :
            self.snp = False

    def toBedFormat(self, add_start=0, add_end=0):

        if self.snp == False:
            ## add logic to handle indels/next time not needed now
            pass


        start = self.pos + add_start
        end = self.pos + add_end
        return "%s\t%s\t%s\t%s" % (self.chr, start, int(end), self.id)


    def __str__(self):
       return str(self.__dict__)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


class VariantFileObject:
    def __init__(self, variantClass):
        self._variant = []
        self.variantClass = variantClass


    def load_variant_file(self, variantFile):
          with open(variantFile) as f:
            snp_id = 0
            for line in f.readlines():
                if not line.startswith("#"):
                    line = line.strip()
                    fields = line.split('\t')
                    snp_id = snp_id + 1
                    v = self.variantClass(chr = fields[0], pos = int(fields[1]), id=snp_id, ref=fields[3], alt=fields[4])
                    self._variant.append(v)

    def drop_indels(self):
        self._variant = [ x for x in self._variant if x.snp is True]

    def write_bed_file(self, bedfile, add_start=0, add_end=0):
        with open(bedfile, 'w') as f:
            for v in self._variant:
                f.write("%s\n" % v.toBedFormat(add_start,add_end))
        return bedfile

def to_mutect(bedtoolsObj, outfile):
    columns = ["snp_id","ref_allele","alt_allele","context"]
    fh_out = open(outfile, 'w+')
    fh_out.write("%s\n" % "\t".join(columns))

    with open(bedtoolsObj.seqfn) as f:
        for header, sequence in itertools.izip_longest(*[f]*2):
            header = header.strip()
            sequence = sequence.strip()
            sequence = sequence[:3] + "x" + sequence[3:]
            header = header.split(';')
            header.append(sequence)
            header.append('KEEP')
            fh_out.write("%s\n" % "\t".join(header))

    fh_out.close()


if __name__ == '__main__':
    import pybedtools

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="inFileName", required=True,
                        help="Input Vcf F")

    parser.add_argument("-b", dest="bedfile", required=True,
                        help="bedfile containing region")

    parser.add_argument("-o", dest="outFileName", required=True,
                        help="the output file")

    parser.add_argument("-v", dest="verbosity", type=int, default=0,
                        help="verbosity level 0-2 [default=0]")

    parser.add_argument("-g", dest="genome", required=True,
                        help="Genome file")


    args = parser.parse_args()


    VF = VariantFileObject(Variant)
    VF.load_variant_file(args.inFileName)
    VF.drop_indels()
    bedfile = VF.write_bed_file(args.bedfile, -3, 3)
    genome_fasta = args.genome
    Bedtools = pybedtools.BedTool(bedfile)
    Bedtools = Bedtools.sequence(genome_fasta, name=True)
    to_mutect(Bedtools, args.outFileName)







