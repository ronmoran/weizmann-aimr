from Bio.Alphabet import DNAAlphabet
from Bio.Data.CodonTable import CodonTable, TranslationError, ambiguous_dna_by_id
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, ExactPosition
from Bio.SeqRecord import SeqRecord

from src.Base import Base

__author__ = "ronmor"


def _cleave_after_stop_codon(nucleotide_seq, table):
    """
    :param nucleotide_seq: The DNA sequence to cleave
    :type nucleotide_seq: Seq
    :param table:
    :type table: CodonTable
    :return:
    :rtype: list[Seq]
    :raises: ValueError if not a nucleotide sequence
    """
    _check_if_nucleotide_sequence(nucleotide_seq)
    nucleotide_seq = nucleotide_seq.upper()
    stop_codons = table.stop_codons
    cleaved_at_stop = []
    codon = ""
    after_stop_index = 0
    for index, letter in enumerate(nucleotide_seq):
        codon += letter
        if codon in stop_codons:
            cleaved_at_stop.append(nucleotide_seq[after_stop_index:index + 1])
            after_stop_index = index + 1
        if len(codon) == 3:
            codon = ""
    return cleaved_at_stop


def _check_if_nucleotide_sequence(nucleotide_sequence):
    if not isinstance(nucleotide_sequence, Seq):
        raise ValueError("Expected a sequence, got  %s instead" % type(nucleotide_sequence))
    elif not isinstance(nucleotide_sequence.alphabet, DNAAlphabet):
        raise ValueError("Expected DNA alphabet, found %s instead" % type(nucleotide_sequence.alphabet))


class DownstreamAnalyzer(object):
    def __init__(self, downstream_sequence, coding_sequence_start_index, genbank_file, is_complementary=False):
        """
        :param downstream_sequence: The sequence downstream of the a gene which we want to analyze.
        :type downstream_sequence: Seq
        :param coding_sequence_start_index: The start index of the gene in the genbank-file features list.
        :type coding_sequence_start_index: int
        :param genbank_file:
        :type genbank_file: SeqRecord
        :param is_complementary:
        :type is_complementary: bool
        :raise: ValueError if the object doesn't hold a genbank documentation or if the genbank file isn't DNA
        """
        self.__downstream_seq = downstream_sequence
        self.__start_index = coding_sequence_start_index
        self. __is_complementing = is_complementary
        if genbank_file is not None:
            if not isinstance(genbank_file, SeqRecord):
                raise ValueError(
                    "genbank file type expected to be of type SeqRecord, found %s instead" % type(self.__genbank_file))
            _check_if_nucleotide_sequence(genbank_file.seq)
        self.__genbank_file = genbank_file

    @property
    def downstream_seq(self):
        return self.__downstream_seq

    @property
    def is_complementing(self):
        return self.__is_complementing

    def find_possible_proteins_in_downstream_sequence(self, table_id=11):
        """
        For the downstream sequence, get all possible CDSs and return their matching proteins.
        ORF is searched for the entire sequence and not for the complement strand (3 total).
        Use only for sequences known not to have introns!
        :param table_id: ID of translation table as appears on https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        Defaults to 11.
        :type table_id: int
        :return: ORF index (0-2) to list of proteins (type Seq)
        :rtype: dict[int, list[Seq]]
        """
        table = ambiguous_dna_by_id[table_id]
        orfs = self._find_possible_orfs_on_positive_strand(self.__downstream_seq)
        possible_proteins = {}
        for orf_num, orf in enumerate(orfs):
            cleaved_frame = _cleave_after_stop_codon(orf, table)
            for possible_cds in cleaved_frame:
                try:
                    possible_protein = possible_cds.translate(table=table, cds=True)
                except TranslationError:
                    continue
                else:
                    possible_proteins[orf_num] = possible_protein
        return possible_proteins

    @staticmethod
    def _find_possible_orfs_on_positive_strand(nucleotide_sequence):
        """
        Only find the orfs in the downstream sequence given, WITHOUT the complementing strand.
        :param nucleotide_sequence: the sequence to find the
        :rtype: [Seq]
        :return:
        """
        for frame in range(3):
            length = 3 * ((len(nucleotide_sequence) - frame) // 3)  # Multiple of three
            yield nucleotide_sequence[frame:frame + length]

    def _find_next_gene_index_in_genbank(self):
        """
        :return:
        :rtype: int
        """
        feature = None
        for feature in self.__genbank_file.features:
            if feature.location.start >= self.__start_index:
                break
        index = self.__genbank_file.features.index(feature) - 1 if feature else None
        if self.__is_complementing:
            return index
        elif index is not None:
            try:
                while self.__start_index < self.__genbank_file.features[index].location.start:
                    index += 1
            # End of list
            except IndexError:
                index -= 1
        return index

    def generate_downstream_cdss(self):
        """
        find all the features (usually CDSs) downstream of the gene, going towards the 3' end in both strands.
        :return:
        """
        next_feature_index = self._find_next_gene_index_in_genbank()
        if next_feature_index is None:
            raise ValueError("Could not find a feature downstream")
        downstream_feature_index = next_feature_index
        while abs(self.__genbank_file.features[downstream_feature_index].location.start -
                  self.__genbank_file.features[next_feature_index].location.start) <= len(self.__downstream_seq):
            if not (self.__genbank_file.features[downstream_feature_index].strand == Base.COMPLEMENT
                    and not self.__is_complementing):
                feature_to_yield = self.__genbank_file.features[downstream_feature_index]
                if hasattr(feature_to_yield, 'type') and feature_to_yield.type == 'CDS' and \
                        isinstance(feature_to_yield.location.start,  ExactPosition) and \
                        isinstance(feature_to_yield.location.end, ExactPosition):
                    yield feature_to_yield
            if self.__is_complementing:
                downstream_feature_index -= 1
            else:
                downstream_feature_index += 1
            if downstream_feature_index < 0 or downstream_feature_index > len(self.__genbank_file.features)-1:
                break

    def translate_feature(self, feature, table=11):
        """
        :param feature: a feature of a genome. has to be RNA or DNA
        :type feature: SeqFeature
        :param table: The table used for translation: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        :type table: int
        :rtype: Seq
        """
        is_cds = True if feature.type == 'CDS' else False
        return feature.extract(self.__genbank_file).seq.translate(table=table, cds=is_cds)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            raise ValueError("Expected instance of %s got %s instead" % (type(self), type(other)))
        return self.__downstream_seq == other.downstream_seq
