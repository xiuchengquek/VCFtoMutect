#!/usr/bin/env python
###############################################################################
# Extract mutational signatures from MuTect calls.
#
# H.C. Lee    04/04/2014
###############################################################################
import argparse
import csv


class Signature(object):
    # a particular signature, one of six classes: C>A, C>G, C>T, T>A, T>C, T>G.
    # Each class has up to 16 sub-types based on the pre- and post-nucleotides
    def __init__(self, classId, subIdList):
        self._id = classId
        self._subClass = {}  # used to hold information about sub-classes

        for subId in subIdList:
            self._subClass[subId] = 0


    ## tally is a counter for the numer of sub-classes
    def tally(self, subId):
        self._subClass[subId] += 1


    ##count the number of subclasses. and convert them as a percentage of total mutation
    def toString(self, totalMuts):
        leStr = ""
        for subId in sorted(self._subClass):
            nMuts = self._subClass[subId]
            mutPerc = float(100 * nMuts) / totalMuts
            leStr += '\t'.join([self._id, subId, str(nMuts),
                                str(mutPerc) + '\n'])

        return leStr



def pre_post_combos(refAl):
    # generate a list of 16 sub-class signatures based on the given reference
    # allele
    subs = []
    nucleotides = ['A', 'C', 'G', 'T']

    for nuc in nucleotides:
        for i in xrange(4):
            subs.append(''.join([nuc, refAl, nucleotides[i]]))

    return subs

## do a re
def reverse_complement(nucleotide):
    assert nucleotide in 'ACGTNx', "Unknown nucleotide!"
    # reverse-complement
    return {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'x': 'x',
        'N': 'N'
    }[nucleotide]


def reverse_info(contextStr, refAl, altAl):
    # return information about the reference strand
    revContStr = ''.join([reverse_complement(c) for c in contextStr[::-1]])
    return (revContStr, reverse_complement(refAl), reverse_complement(altAl))




def parse_mutational_context(inFileName, outFileName):
    # if ref_allele is a Purine (A or G), reverse_info is called to simplify
    # the classes of mutational context.
    nMuts = 0  # total number of mutations
    mutSig = {}  # dict of mutational signature Id (e.g. C>A) to Signature

    with open(inFileName, 'rU') as inFile:
        inFile.readline()  # remove first line

        for row in csv.DictReader(inFile, dialect='excel-tab'):
            if row['judgement'] != 'KEEP':
                continue
            nMuts += 1
            context = row['context']
            refAl = row['ref_allele']
            altAl = row['alt_allele']

            if refAl == 'A' or refAl == 'G':
                (context, refAl, altAl) = reverse_info(context, refAl, altAl)

            sigId = ''.join([refAl, '>', altAl])
            subId = ''.join([context[2], refAl, context[4]])

            try:
                mutSig[sigId].tally(subId)
            except KeyError:
                mutSig[sigId] = Signature(sigId, pre_post_combos(refAl))
                mutSig[sigId].tally(subId)

    with open(outFileName, 'w') as outFile:
        for sigId in sorted(mutSig):
            outFile.write(mutSig[sigId].toString(nMuts))
            #outFile.write('\t'.join([context, refAl, altAl + '\n']))






if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", dest="inFileName", required=True,
                        help="the input MuTect call_stats.out file")
    parser.add_argument("-o", dest="outFileName", required=True,
                        help="the output file")
    parser.add_argument("-v", dest="verbosity", type=int, default=0,
                        help="verbosity level 0-2 [default=0]")

    args = parser.parse_args()

    parse_mutational_context(args['inFileName'] , args['outFileName'])


