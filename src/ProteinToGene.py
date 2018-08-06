import re
import json
from Bio import SeqRecord, SeqFeature
from src.Base import Base, DiminishingDict

__author__ = "ronmor"


class ProteinToGene(Base):
    _accession_to_object = DiminishingDict()
    __location_format = re.compile("(?:(?P<complement>complement)\(|^)(?P<accession>[A-Z_0-9]*\.\d):" +
                                   "(?P<start>\d+)\.{2,}(?P<end>\d+)")

    def __init__(self, accession_number, db_type, similar_proteins_accessions=()):
        super(ProteinToGene, self).__init__(accession_number)
        self.__dont_fetch = False
        if db_type not in ('emb', 'gb', 'dbj'):
            self.__dont_fetch = True
        self.__db_type = db_type
        self.__genpept = None
        self._similar_proteins_accessions = similar_proteins_accessions
        self.__genome_accession = None

    def accession_no(self):
        return self._accession_no

    @property
    def genome_accession_no(self):
        return self.__genome_accession

    def _set_genome_accession(self, genome_accession):
        self.__genome_accession = genome_accession

    @property
    def similar_proteins_accessions(self):
        """
        :return:
        :rtype: list[str]
        """
        return self._similar_proteins_accessions

    def set_genpept_by_accession(self):
        """
        :return:
        :rtype: _io.TextIOWrapper
        """
        if self.__dont_fetch:
            print("The protein is found on a non supported DB (got '%s'). EMBL and genbank are supported. "
                  "Try similar proteins" % self.__db_type)
            print("similar proteins: %s" % ', '.join(self._similar_proteins_accessions))
        self._set_genpept(self._search_entrez_in_genbank_format(db='protein'))

    def _set_genpept(self, entrez_genpept_file):
        self.__genpept = entrez_genpept_file

    def set_genome_accession_no(self, protein_from_genpept=None):
        """
        :param protein_from_genpept:
        :type protein_from_genpept: SeqRecord
        """
        if not protein_from_genpept:
            protein_from_genpept = self.__genpept
        genome_accession = None
        if hasattr(protein_from_genpept, 'annotations'):
            genome_accession = self._accession_format.search(protein_from_genpept.annotations['db_source']).group()
        self.__genome_accession = genome_accession

    def __eq__(self, other_protein):
        """
        :param other_protein:
        :type other_protein: ProteinToGene
        :return:
        """
        if isinstance(other_protein, ProteinToGene) and self._accession_no and self.__genome_accession and \
            self._accession_no == other_protein.accession_no and \
                self.__genome_accession == other_protein.genome_accession_no:
            return True
        else:
            return False

    def __hash__(self):
        """
        :return:
        """
        return hash(self.__genome_accession) if self.__genome_accession else super(ProteinToGene, self).__hash__()

    def extract_location_on_gene(self):
        """
        :return:
        :rtype: SeqFeature.FeatureLocation
        """
        if self.__genpept is None:
            return
        cds_set = set(filter(lambda feature: feature.type == 'CDS', self.__genpept.features))
        if len(cds_set) != 1:
            return
        cds = cds_set.pop()
        coding_regions = cds.qualifiers.get('coded_by', [])
        if len(coding_regions) != 1:
            return
        coding_regions = coding_regions[0]
        res = self.__location_format.search(coding_regions)
        if res is None:
            return
        if self._accession_no is None:
            self._accession_no = res.group("accession")
        strand = self.COMPLEMENT if bool(res.group("complement")) else None
        start = int(res.group("start"))
        end = int(res.group("end"))
        return SeqFeature.FeatureLocation(start, end, strand)

    def set_entrez_file(self, entrez_file):
        self._set_genpept(entrez_file)

    @classmethod
    def fetchall(cls):
        cls.search_multiple_with_entrez_genbank_format(db='protein')

    def to_json(self):
        jsoned_dict = {"similar_protein_accessions": self._similar_proteins_accessions,
                       "accession_number": self._accession_no,
                       "db_type": self.__db_type,
                       "genome_accession": self.__genome_accession}
        return json.dumps(jsoned_dict)

    @classmethod
    def from_json(cls, json_file):
        similar_proteins = json_file["similar_protein_accessions"]
        db_type = json_file["db_type"]
        accession_number = json_file["accession_number"]
        obj = cls(accession_number, db_type, similar_proteins)
        obj._set_genome_accession(json_file["genome_accession"])
        return obj

    @classmethod
    def list_to_json(cls, objs_iter):
        """
        json all objects that are instances of this class
        :param objs_iter: iterable of objects of this class: set, tuple, list etc... 
        :return:
        """
        json_obj_list = []
        for obj in objs_iter:
            if isinstance(obj, cls):
                json_obj_list.append(obj.to_json())
        return json.dumps(json_obj_list)
