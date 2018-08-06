import re
from abc import ABCMeta, abstractmethod
from Bio import Entrez, SeqIO, SeqRecord
from tqdm import tqdm
__author__ = "ronmor"


class DiminishingDict(dict):
    def __iter__(self):
        self.previous_key = None
        self.dict_iterator = iter(self.keys())
        return self

    def __next__(self):
        finished = False
        try:
            key = next(self.dict_iterator)
        except StopIteration as e:
            finished = True
            raise e
        finally:
            if self.previous_key:
                self[self.previous_key] = None
            if finished:
                self.clear()
        self.previous_key = key
        return key


class Base(object):
    __metaclass__ = ABCMeta
    _accession_format = re.compile("([A-Z_0-9]*\.\d)")
    COMPLEMENT = -1
#    _accession_to_object = DiminishingDict()

    def __init__(self, accession_number):
        if accession_number is None:
            raise ValueError("Accession number must not be None")
        self._accession_no = accession_number
        self._add_object_to_accessions_dict(accession_number)

    @property
    @abstractmethod
    def _accession_to_object(self):
        """should map accession number to object of children classes. Class variable"""
        pass

    @classmethod
    def get_accession_to_object(cls):
        """
        Iterating over this dict will null its values, so iterate carefully!
        :rtype: dict [str, Base]
        """
        return cls._accession_to_object

    @property
    def accession_no(self):
        return self._accession_no

    def _search_entrez_in_genbank_format(self, db, email='ron.moran@weizmann.ac.il', seq_start=None, seq_stop=None):
        """
        :param db:
        :param email:
        :param seq_start:
        :param seq_stop:
        :return:
        :rtype: SeqRecord
        """
        Entrez.email = email
        handle = Entrez.efetch(db=db, id=self._accession_no, rettype='gb', retmode='text', seq_start=seq_start,
                               seq_stop=seq_stop)
        try:
            return SeqIO.read(handle, 'gb')
        except ValueError:
            print("More or less than one record was received, aborting Entrez search")

    @classmethod
    def search_multiple_with_entrez_genbank_format(cls, db, return_format='gb', email='ron.moran@weizmann.ac.il'):
        """
        :param db:
        :param return_format:
        :type return_format:
        :param email:
        :return:
        """
        Entrez.email = email
        if not isinstance(cls._accession_to_object, DiminishingDict):
            raise ValueError("You know...")
        handle = Entrez.efetch(db=db, id=list(cls._accession_to_object.keys()), rettype=return_format)
        fetched_files_iterator = SeqIO.parse(handle, return_format)
        for fetched_file in tqdm(fetched_files_iterator):
            file_id = fetched_file.id
            obj = cls._accession_to_object[file_id]
            if not isinstance(obj, cls):
                continue
            else:
                obj.set_entrez_file(fetched_file)

    def _add_object_to_accessions_dict(self, accession):
        self._accession_to_object[accession] = self

    @abstractmethod
    def set_entrez_file(self, efetched_file):
        pass

    @abstractmethod
    def fetchall(self):
        pass
