# coding=utf-8
"""Contains models related with GO Terms"""

# python standard imports
import gzip
import re
from urllib.request import urlretrieve

# django imports
from django.db import models, transaction

from pseudoenzymes.settings import GENE_ONTOLOGY_FILE


class TermQuerySet(models.QuerySet):
    """Some predefined querysets for the GoTerm model"""

    def catalytic(self):
        """return all catalytic go terms"""
        return self.children_of(Term.objects.filter(name="catalytic activity"))

    def children_of(self, parents):
        """recursively find all terms that are children of these parents

        parents are included in the resulting query
        """
        ids = set(parents.values_list("id", flat=True))
        no_previous_ids = 0
        while no_previous_ids != len(ids):
            no_previous_ids = len(ids)
            ids.update(Term.objects.filter(relation1__term2__id__in=ids)\
                       .values_list("id", flat=True))
        return self.filter(id__in=ids)


class Term(models.Model):
    """Gene Ontology Terms"""
    id = models.IntegerField(primary_key=True)
    name = models.TextField()
    definition = models.TextField()
    ASPECTS = {"molecular_function", "biological_process", "cellular_component"}
    aspect = models.TextField(choices=[(aspect, aspect) for aspect in ASPECTS])

    objects = TermQuerySet.as_manager()

    def __repr__(self):
        return self.code

    def __str__(self):
        return f"{self.code} - {self.name}"

    @property
    def code(self):
        """Returns a go term in this format GO:0001234"""
        return f"GO:{self.id:07}"

    @classmethod
    def download_ontology(cls):
        """Downloads the ontology file"""
        urlretrieve("http://purl.obolibrary.org/obo/go.obo", GENE_ONTOLOGY_FILE)

    @classmethod
    def create_from_ontology_file(cls):
        """Reads the ontology file and adds all GO terms to the database"""
        inside_quotes = re.compile(r'"[^"\\]*(?:\\.[^"\\]*)*"')
        objs = []
        with open(GENE_ONTOLOGY_FILE, 'r') as obo_file:
            info = {}
            for line in obo_file:
                if line == "\n" and "id" in info:
                    objs.append(cls(**info))
                    info = {}
                elif line.startswith("id: GO:"):
                    info["id"] = int(line.replace("id: GO:", ""))
                elif line.startswith("def:"):
                    info["definition"] = inside_quotes.findall(line)[0][1:-1]\
                        .replace("\\\"", "\"").replace("\\n", "\n")
                elif line.startswith("namespace"):
                    info["aspect"] = line.split()[1]
                elif line.startswith("name:"):
                    info["name"] = line.split(":", maxsplit=1)[1].strip()
        cls.objects.bulk_create(objs)
        # cls.objects.catalytic_filter().update(is_catalytic=True)


class Relation(models.Model):
    """Relationship between go terms"""
    term1 = models.ForeignKey(
            "Term",
            on_delete=models.CASCADE,
            related_name="relation1",
            db_index=True
            )
    term2 = models.ForeignKey(
            "Term",
            on_delete=models.CASCADE,
            related_name="relation2",
            db_index=True
            )
    relation = models.CharField(max_length=255, db_index=True)

    objects = models.Manager()

    def __repr__(self):
        return f"{self.term1} {self.relation} {self.term2}"

    def __str__(self):
        return self.__repr__()

    class Meta:
        unique_together = ['term1', 'relation', 'term2']

    @classmethod
    @transaction.atomic
    def create_from_ontology_file(cls):
        """Read the ontology file and add all the is_a GO relations to the database"""
        objs = []
        with open(GENE_ONTOLOGY_FILE, 'r') as obo_file:
            term1 = None
            for line in obo_file:
                if line.startswith("[Term]"):
                    term1 = "reading"
                elif line.startswith("id:") and term1 == "reading":
                    term1 = int(line.replace("id: GO:", ""))
                elif line.startswith("is_a") and term1:
                    term2 = int(line.replace("is_a: GO:", "").split("!")[0])
                    objs.append(cls(term1_id=term1,term2_id=term2, relation="is_a"))
                elif line == "\n":
                    term1 = None
        cls.objects.bulk_create(objs)


# class SequenceGoQuerySet(models.QuerySet):
    # """Some predefined querysets for the SequenceGo associations"""

    # def experimental(self):
        # """filters associations supported by experimental evidence"""
        # return self.filter(eco_term__in=common.models.EcoTerm.objects.experimental())
    
    # def enables(self, enables=True):
        # """filters associations with the enable qualifier"""
        # return self.filter(qualifier="enables" if enables else "NOT|enables")
    
    # def catalytic(self, enables=True):
        # """filters associations to catalytic go terms"""
        # return self.filter(go_term__in=common.models.GoTerm.objects.catalytic())


# class SequenceGo(models.Model):
    # """Through table to link go with sequence"""
    # go_term = models.ForeignKey("GoTerm", on_delete=models.CASCADE)
    # sequence = models.ForeignKey("Sequence", on_delete=models.CASCADE)
    # eco_term = models.ForeignKey("EcoTerm", on_delete=models.CASCADE)
    # qualifier = models.CharField(max_length=127, db_index=True)

    # objects = SequenceGoQuerySet.as_manager()

    # # noinspection PyMissingOrEmptyDocstring
    # class Meta:
        # unique_together = ["go_term", "sequence", "qualifier", "eco_term"]

    # @classmethod
    # def create_objects(cls):
        # """Populates this table by reading from the gpa file"""
        # info_tuples = cls._read_gpa_tuples()
        # print(len(info_tuples))
        # existing = set(cls.get_unique_together_values())
        # to_delete = existing - info_tuples
        # print("Should be deleting {} Sequence<->Go links".format(len(to_delete)))
        # for t in to_delete:
            # print(cls.objects.get(go_term_id=t[0], sequence_id=t[1], qualifier=t[2], eco_term_id=t[3]).delete())
        # to_create = info_tuples - existing
        # print("Creating {} Sequence<->Go links".format(len(to_create)))
        # cls.objects.bulk_create(
            # [cls(go_term_id=t[0], sequence_id=t[1], qualifier=t[2], eco_term_id=t[3]) for t in to_create],
            # batch_size=5000
        # )


    # @classmethod
    # def _read_gpa_tuples(cls):
        # """Reads information from gpa file and return tuples with information to add"""
        # from common.models import Sequence, GoTerm
        # seq_id = dict(Sequence.objects.reviewed().get_idx2id())
        # go_id = dict(GoTerm.objects.get_idx2id())
        # eco_id = dict(EcoTerm.objects.get_idx2id())
        # info = set()
        # with gzip.open(files.GO_GCRP_GPA, "rt") as gpa_file:
            # for line in gpa_file:
                # if line.startswith("UniProtKB"):
                    # words = line.split()
                    # uniprot_id = words[1]
                    # if uniprot_id in seq_id:
                        # qualifier = words[2]
                        # go_term = int(words[3].split(":")[1])
                        # eco_term = words[5].split(":")[1].strip()
                        # if go_term in go_id:
                            # info.add((go_id[go_term], seq_id[uniprot_id], qualifier, eco_id[eco_term]))
        # return info

