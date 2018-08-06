from Bio import SeqFeature, SeqRecord

from src.Base import Base, DiminishingDict

__author__ = "ronmor"


class ProteinOnDNA(Base):
    _accession_to_object = DiminishingDict()

    def __init__(self, accession_number, protein_location=None, protein_accession_no=None):
        super(ProteinOnDNA, self).__init__(accession_number)
        self.protein_location = protein_location
        self.__protein_accession_no = protein_accession_no
        self.__genbank_file = None

    @property
    def protein_accession_no(self):
        return self.__protein_accession_no

    @property
    def genbank_file(self):
        return self.__genbank_file

    def set_genbank_by_accession(self, start=None, stop=None):
        """
        :return:
        :rtype: _io.TextIOWrapper
        """
        if not self._accession_no:
            print("No accession number to search by. Quiting...")
            return
        nuccore = self._search_entrez_in_genbank_format('nuccore', seq_start=start, seq_stop=stop)
        if not nuccore:
            self._set_genbank_file(self._search_entrez_in_genbank_format('nucleotide', seq_start=start, seq_stop=stop))
        else:
            self._set_genbank_file(nuccore)

    def find_protein_location(self, protein_accession=None, force_search=False):
        if isinstance(self.protein_location, SeqFeature.FeatureLocation) and not force_search:
            return
        if protein_accession is None:
            protein_accession = self.__protein_accession_no
        if protein_accession is None:
            return
        if not isinstance(self.__genbank_file, SeqRecord.SeqRecord):
            self.set_genbank_by_accession()
        if not isinstance(self.__genbank_file, SeqRecord.SeqRecord):
            raise ValueError("Genbank file was not acquired. Cannot find protein location")
        feature = None
        for feature in self.__genbank_file.features:
            if 'protein_id' in feature.qualifiers:
                if protein_accession in feature.qualifiers['protein_id']:
                    break
        return feature.location

    def get_next_n_nucleotide_of_feature_location(self, n, location=None):
        """
        :param n:
        :param location:
        :return:
        """
        if location is None:
            location = self.protein_location
        if location:
            if location.strand == self.COMPLEMENT:
                start, stop = self._get_reverse_indices(n, location)
            else:
                start = location.end
                stop = start + n
            if self.__genbank_file is None:
                self.set_genbank_by_accession(start=start + 1, stop=stop)
                seq = self.__genbank_file.seq
            elif isinstance(self.__genbank_file, SeqRecord.SeqRecord):
                seq = self.__genbank_file.seq[start:stop]
            else:
                return
            return seq if location.strand != self.COMPLEMENT else seq.reverse_complement()

    @staticmethod
    def _get_reverse_indices(downstream_nucleotide_count, location):
        end_index = location.start - 1 if location.start > 0 else 0
        start_index = end_index - downstream_nucleotide_count
        start_index = start_index if start_index >= 0 else 0
        return start_index, end_index

    @classmethod
    def get_end_index_of_feature_location(cls, location):
        """
        :param location:
        :return:
        """
        if location.strand == cls.COMPLEMENT:
            start, stop = cls._get_reverse_indices(0, location)
        else:
            start = location.end
        return start

    def __hash__(self):
        return hash(self._accession_no)

    def _set_genbank_file(self, genbank_file):
        self.__genbank_file = genbank_file

    def set_entrez_file(self, efetched_file):
        self._set_genbank_file(efetched_file)

    @classmethod
    def fetchall(cls):
        cls.search_multiple_with_entrez_genbank_format(db='nucleotide')
