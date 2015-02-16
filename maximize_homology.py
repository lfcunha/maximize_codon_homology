#!/usr/bin/env python

__author__ = 'Luis Cunha, PhD'
__DATE__ = "Feb 15th 2015"

'''
Description:
Script to Maximize the DNA sequence homology of protein homologues

it requires BioPython for reading inputs

Input:
    1) fasta file with multiple dna sequences
    2) .aln (clustalw) file with protein multiple seq alignment (with gaps)
    3) output prefix
    4) optional -v flag: prints result to console

This script maximizes sequence similarity of the 3rd nucleotide of each codon.
It optimizes not only same aminoacid codons, as well as codons coding different aa.

After optimziation, it validates that optimized dna sequence codes the same aa sequence as the 
original sequence, throwing an error otherwise

As is, it does not optimize the first and second positions of the codon.

Only Ser, Arg, and Leu have 6 codons, with alternative 1st/2nd positions. However,
both positions in 2 of the codons are different from the other four, requiring that both
nucleotides are taken into consideration. The frequency of an optimization being possible
is too low to justify it

Usage: maximize_homology [-h] [--dna FASTA] [--protein_alignment ALIGNMENT] [--output OUTPUT_PREFIX]
'''

import operator
import sys


class TranslateSequence:
    '''
	translate the given DNA sequence to protein sequence
    '''
    def __init__(self):
        self.codon_table = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA":"T","ACC":"T","ACG":"T","ACT":"T","AGA":"R","AGC":"S","AGG":"R","AGT":"S","ATA":"I","ATC":"I","ATG":"M","ATT":"I","CAA":"Q","CAC":"H","CAG":"Q","CAT":"H","CCA":"P","CCC":"P","CCG":"P","CCT":"P","CGA":"R","CGC":"R","CGG":"R","CGT":"R","CTA":"L","CTC":"L","CTG":"L","CTT":"L","GAA":"E","GAC":"D","GAG":"E","GAT":"D","GCA":"A","GCC":"A","GCG":"A","GCT":"A","GGA":"G","GGC":"G","GGG":"G","GGT":"G","GTA":"V","GTC":"V","GTG":"V","GTT":"V","TAA":"","TAC":"Y","TAG":"","TAT":"Y","TCA":"S","TCC":"S","TCG":"S","TCT":"S","TGA":"","TGC":"C","TGG":"W","TGT":"C","TTA":"L","TTC":"F","TTG":"L","TTT":"F"}

    def translate(self, sequence_):
        g=sequence_
        splitted_g = [g[0+i:3+i] for i in range(0, len(g), 3)]
        if len(splitted_g[-1]) != 3:
            splitted_g = splitted_g[:-1]
        trans = [self.codon_table[c] for c in splitted_g]
        return "".join(trans)



class optimizeHomology():
    def __init__(self, proteins, dnas, translator):
        self.trans = translator
        self.proteins = proteins
        self.dnas= dnas
        #assert that the length of the protein (without gaps) is the same as the protein resulting from the dna translation => make sure no deletions introduced in processing sequences (e.g. accidental deletion of a line)
        for sequenceName in self.proteins:
            assert len("".join([aa for aa in self.proteins[sequenceName] if aa != "-"])) == len(self.trans.translate(self.dnas[sequenceName])), "error matching protein sequence to translated dna sequence: " + sequenceName
        #ditionary for optimized sequence
        self.optimizedSequences = {protein: [] for protein in self.proteins}
        #auxiliary dictionaries
        #self.invertedCodonTable = {'': ['tag', 'taa', 'tga'], 'A': ['gca', 'gcc', 'gcg', 'gct'], 'C': ['tgt', 'tgc'], 'E': ['gag', 'gaa'], 'D': ['gat', 'gac'], 'G': ['ggt', 'ggg', 'gga', 'ggc'], 'F': ['ttt', 'ttc'], 'I': ['atc', 'ata', 'att'], 'H': ['cat', 'cac'], 'K': ['aag', 'aaa'], 'M': ['atg'], 'L': ['ctc', 'ctg', 'cta', 'ctt', 'tta', 'ttg'], 'N': ['aac', 'aat'], 'Q': ['caa', 'cag'], 'P': ['cct', 'ccg', 'ccc', 'cca'], 'S': ['agc', 'agt', 'tcg', 'tca', 'tcc', 'tct'], 'R': ['agg', 'aga', 'cga', 'cgc', 'cgg', 'cgt'], 'T': ['acc', 'aca', 'act', 'acg'], 'W': ['tgg'], 'V': ['gta', 'gtc', 'gtg', 'gtt'], 'Y': ['tat', 'tac']}
        self.aa_nt_space_per_position = {'': {0: ['t'], 1: ['a', 'g'], 2: ['g', 'a']}, 'A': {0: ['g'], 1: ['c'], 2: ['a', 'c', 'g', 't']}, 'C': {0: ['t'], 1: ['g'], 2: ['t', 'c']}, 'E': {0: ['g'], 1: ['a'], 2: ['g', 'a']}, 'D': {0: ['g'], 1: ['a'], 2: ['t', 'c']}, 'G': {0: ['g'], 1: ['g'], 2: ['t', 'g', 'a', 'c']}, 'F': {0: ['t'], 1: ['t'], 2: ['t', 'c']}, 'I': {0: ['a'], 1: ['t'], 2: ['c', 'a', 't']}, 'H': {0: ['c'], 1: ['a'], 2: ['t', 'c']}, 'K': {0: ['a'], 1: ['a'], 2: ['g', 'a']}, 'M': {0: ['a'], 1: ['t'], 2: ['g']}, 'L': {0: ['c', 't'], 1: ['t'], 2: ['c', 'g', 'a', 't']}, 'N': {0: ['a'], 1: ['a'], 2: ['c', 't']}, 'Q': {0: ['c'], 1: ['a'], 2: ['a', 'g']}, 'P': {0: ['c'], 1: ['c'], 2: ['t', 'g', 'c', 'a']}, 'S': {0: ['a', 't'], 1: ['g', 'c'], 2: ['c', 't', 'g', 'a']}, 'R': {0: ['a', 'c'], 1: ['g'], 2: ['g', 'a', 'c', 't']}, 'T': {0: ['a'], 1: ['c'], 2: ['c', 'a', 't', 'g']}, 'W': {0: ['t'], 1: ['g'], 2: ['g']}, 'V': {0: ['g'], 1: ['t'], 2: ['a', 'c', 'g', 't']}, 'Y': {0: ['t'], 1: ['a'], 2: ['t', 'c']}}

    def countNucleotideOccurrences(self, local_codon_, ntIndex):
        localcodon=local_codon_
        count = {"a": 0, "t": 0, "c": 0, "g": 0}
        for protein in self.proteins:
            if localcodon[protein][ntIndex] == "-":
                continue
            count[localcodon[protein][ntIndex]] += 1
        return count

    def optimize(self):
        # DNAindex for each sequence (gaps in the protein sequence do not advance the DNA sequence)
        DNAindex = {protein:0 for protein in self.proteins}
        # run through the length of the protein, one aa at a time
        # and optimize each codon for each sequence
        proteinNames = [protein for protein in self.proteins]
        for aminoacidIndex in range(len(proteins[proteinNames[0]])):
            localAA={}
            localCodon={}
            for sequenceName in self.proteins:
                aa = self.proteins[sequenceName][aminoacidIndex]
                localAA[sequenceName] = aa
                if aa != "-":
                    localCodon[sequenceName] = self.dnas[sequenceName][DNAindex[sequenceName]:DNAindex[sequenceName]+3].lower()
                    DNAindex[sequenceName] += 3
                else:
                    localCodon[sequenceName] = "---"

            # Potentially CAN OPTIMIZE EACH NT of the CODON
            # however, since only Ser, Arg, Leu have more than one alternative in the first position, and only one of them has more than one in the second
            # and allowing opimization of those nt require considering both the first and second position together
            # skip optimization of the first two nt and focus on the 3rd position
            for ntIndex in range(2, 3): #3rd codon only
                #FIRST GET COUNT OF EACH NT AT THAT POSITION
                count = self.countNucleotideOccurrences(localCodon, ntIndex)
                sorted_count = sorted(count.items(), key=operator.itemgetter(1), reverse=True)

                # FOR EACH  SEQUENCE, OPTIMIZE THE NT AT THE CURRENT INDEX
                # CODE HAS THE SIDE EFFECT OF MODIFYING THE localCodon VARIABLE
                # APPEND OPTIMIZED localCodon TO GROWING SEQUENCE IN optimizedSequences
                for sequenceName in proteins:
                    nt = localCodon[sequenceName][ntIndex]
                    if nt == "-": continue
                    if count[nt] == len(localCodon): continue
                    for j in range(len(sorted_count)):
                        nt = localCodon[sequenceName][ntIndex]
                        nt_ = sorted_count[j][0]
                        if nt == nt_: break
                        if nt_ in self.aa_nt_space_per_position[localAA[sequenceName]][ntIndex]:
                            if (count[nt_] > count[nt]):
                                strList = list(localCodon[sequenceName])
                                strList[ntIndex] = nt_
                                #Ser, leu and arg, which have 6 codons, divided in two groups according to the first 2 nt positions
                                #for one of the groups, there's only two nt options for the 3rd option, so we can't let the 3rd position have all four nt options for every casde
                                # we would have to consider what the first two nt are to decide on the nt possibilities on the 3rd position
                                # alternatively, I let it take any option, but then check if the aa coded is the same as original
                                # and reject if not
                                if self.trans.translate("".join(strList).upper()) == proteins[sequenceName][aminoacidIndex]:
                                    localCodon[sequenceName] = "".join(strList)
            for sequenceName in proteins:
                self.optimizedSequences[sequenceName].append(localCodon[sequenceName])

    @property
    def optimizedSequences(self):
        return self.__optimizedSequences

def usage():
    print "\n\tScript to maximize the homology of protein sequences\n"
    print "\tRequirements: \n\t\tfasta input file with dna sequences\n\t\t.aln file of t-coffee alignment (clustalw format)"
    print "\tUsage: maximize_homology [-h] [--dna FASTA] [--protein_alignment ALIGNMENT] [--output OUTPUT_PREFIX] \n\n"
    sys.exit(0)

if __name__ == "__main__":
    import argparse
    from Bio import SeqIO
    from Bio import AlignIO

    trans = TranslateSequence()

    #get arguments
    parser = argparse.ArgumentParser(description='Optimize sequence homology while retaining aminoacid sequence')
    parser.add_argument('--dna', "-d", help="fasta file with dna sequences")
    parser.add_argument('--protein_alignment', "-p", help=".aln clustalw sequence alignment")
    parser.add_argument('--output', "-o", help="output prefix")
    parser.add_argument("--verbose", "-v", help="print output to console", action='store_true')
    args = parser.parse_args()
    if args.dna==None or args.protein_alignment==None or args.output==None:
        usage()
    if args.dna!=None and args.protein_alignment!=None and args.output!=None:
        fastaInput = args.dna
        proteinAlignment = args.protein_alignment
        output = args.output

    #read input files with biopython help and build dna and proteinAlignment dictionaries
    def readInputs(dnaFile, proteinAlignmentFile):
        dnas = {}

        try:
            f =open(dnaFile, "rU")
        except:
            print "Can't open " + dnaFile + " file\n"
            sys.exit(0)
        else:
            with f as handle:
                for record in SeqIO.parse(handle, "fasta") :
                    dnas[record.id]=str(record.seq.upper())

        proteins = {}
        align = AlignIO.read(proteinAlignmentFile, "clustal")
        for i in range(len(align)):
            proteins[align[i].id] = align[i].seq
        return dnas, proteins

    (dnas, proteins) = readInputs(fastaInput, proteinAlignment)

    #optimize sequence homology
    optimizer = optimizeHomology(proteins, dnas, trans)
    optimizer.optimize()
    optimizedSequences=optimizer.optimizedSequences

    originalSequenceOutput = output + "_original.fasta"
    optimizedSequenceOutput = output + "_optimized.fasta"

    #empty file of exists
    file = open(originalSequenceOutput, "w")
    file.close()
    file = open(optimizedSequenceOutput, "w")
    file.close()

    try:
        f = open(originalSequenceOutput, "ab")
    except:
        print "Can't open " + originalSequenceOutput + " file for writing \n"
        sys.exit(0)
    else:
        try:
            g = open(optimizedSequenceOutput, "ab")
        except:
            print "Can't open " + optimizedSequenceOutput + " file for writing\n"
            sys.exit(0)
        else:
            with f as original:
                with g as optimized:
                    for sequenceName in proteins:
                        #assert that the protein obtained after optimization has the same sequence has the initial protein
                        assert trans.translate("".join([x for x in " " + "".join(optimizedSequences[sequenceName]).upper() if x!="-"]).strip()) == "".join([x for x in proteins[sequenceName] if x !="-"]), "error at sequenceName: " + sequenceName
                        header = ">" + sequenceName + "\n"
                        original.write(header)
                        optimized.write(header)
                        original.write(dnas[sequenceName] + "\n")
                        optimized.write("".join([x for x in " " + "".join(optimizedSequences[sequenceName]).upper() if x!="-"]).strip() + "\n")

    def printToConsole():
        print ">" + sequenceName +" " + "protein sequence translated from optimized dna sequence"
        print trans.translate("".join([x for x in " " + "".join(optimizedSequences[sequenceName]).upper() if x!="-"]).strip())
        print ">" + sequenceName + " " + "original protein sequence"
        print "".join([x for x in proteins[sequenceName] if x !="-"])
        print ">" + sequenceName +" optimized"    #" " + "optimized dna sequence"
        print "".join([x for x in " " + "".join(optimizedSequences[sequenceName]).upper() if x!="-"]).strip()
        print ">" + sequenceName +" original"     #" " + "original dna sequence"
        print dnas[sequenceName]

    if args.verbose:
    	printToConsole()



