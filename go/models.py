# coding=utf-8
"""Contains models related with GO Terms"""

# python standard imports
import gzip
import re
from urllib.request import urlretrieve

# django imports
from django.db import models, transaction

from pseudoenzymes.settings import GENE_ONTOLOGY_FILE, GO_GPA_FILE
import uniprot.models as uniprot
import eco.models as eco


class TermQuerySet(models.QuerySet):
    """Some predefined querysets for the GoTerm model"""

    def catalytic(self):
        """return all catalytic go terms"""
        return self.children_of(Term.objects.filter(name="catalytic activity"))

    def functional(self):
        """go terms related with the function of the protein"""
        return self.children_of(Term.objects.filter(name="molecular_function"))

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

    def go_graph(self):
        """return graph representing the relationships among go terms"""
        import networkx as nx
        graph = nx.DiGraph()
        for relation in Relation.objects.filter(relation="is_a"):
            graph.add_edge(relation.term2_id, relation.term1_id)
        return graph

    def go_to_ancestors(self):
        """return dict of go terms pointing to ancestors"""
        import networkx as nx
        term_to_ancestors = {}
        graph = self.go_graph()
        for term in self:
            if term.id in graph:
                term_to_ancestors[term.id] = nx.ancestors(graph, term.id)
        return term_to_ancestors


class Term(models.Model):
    """Gene Ontology Terms"""
    id = models.IntegerField(primary_key=True)
    name = models.TextField()
    definition = models.TextField()
    ASPECTS = {"molecular_function", "biological_process", "cellular_component"}
    aspect = models.TextField(choices=[(aspect, aspect) for aspect in ASPECTS])
    uniprot_entries = models.ManyToManyField(
            "uniprot.Entry",
            through="TermUniProtEntry",
            related_name="go_terms"
            )

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
        existing = set(cls.objects.values_list("id", flat=True))
        objs = []
        with open(GENE_ONTOLOGY_FILE, 'r') as obo_file:
            info = {}
            for line in obo_file:
                if line == "\n" and "id" in info:
                    if info["id"] not in existing:
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
        print(f"Creating {len(objs)} Go Terms")
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
        # TODO remove old relations
        existing = set(cls.objects.values_list("term1", "relation", "term2"))
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
                    if (term1, "is_a", term2) not in existing:
                        objs.append(cls(term1_id=term1,term2_id=term2, relation="is_a"))
                elif line == "\n":
                    term1 = None
        cls.objects.bulk_create(objs)
        print(f"Creating {len(objs)} go term relations")


class TermUniProtEntryQuerySet(models.QuerySet):
    """Some predefined querysets for the GoTerm model"""

    def catalytic(self):
        """return all catalytic go uniprot associations"""
        return self.filter(term__in=Term.objects.catalytic(), qualifier="enables")

    def not_catalytic(self):
        """return associations with proof of non catalytic activity for a given reaction

        the protein might have other catalytic activities
        this does not include proteins that were not tested to be catalytic
        """
        return self.filter(term__in=Term.objects.catalytic(), qualifier="NOT|enables")

    def functional(self):
        """associations related with the function of the protein"""
        return self.filter(term__in=Term.objects.functional())

    def experimental(self):
        """return all experimentally supported associations"""
        return self.filter(eco_term__in=eco.Term.objects.experimental())

class TermUniProtEntry(models.Model):
    """Through table to link go with UniProt entries"""
    term = models.ForeignKey(
            "Term",
            on_delete=models.CASCADE,
            related_name="uniprot_associations")
    uniprot_entry = models.ForeignKey(
            "uniprot.Entry",
            on_delete=models.CASCADE,
            related_name="go_associations")
    eco_term = models.ForeignKey("eco.Term", on_delete=models.CASCADE)
    qualifier = models.CharField(max_length=127, db_index=True)

    objects = TermUniProtEntryQuerySet.as_manager()

    class Meta:
        unique_together = ["term", "uniprot_entry", "qualifier", "eco_term"]

    @classmethod
    def create_from_gpa_file(cls):
        """Populates this table by reading from the gpa file"""
        objs_data = cls.read_gpa_tuples()
        objs = [cls(term_id=d[0], uniprot_entry_id=d[1],
                qualifier=d[2], eco_term_id=d[3]) for d in objs_data]
        print(f"Creating {len(objs)} Uniprot<->Go links")
        cls.objects.bulk_create(objs, batch_size=100000)


    @classmethod
    def read_gpa_tuples(cls):
        """Read gpa file and return tuples with data to add"""
        acs = set(uniprot.Entry.objects.values_list("ac", flat=True))
        terms = set(Term.objects.values_list("id", flat=True))
        info = set()
        with gzip.open(GO_GPA_FILE, "rt") as gpa_file:
            for line in gpa_file:
                if line.startswith("UniProtKB"):
                    words = line.split()
                    uniprot_id = words[1]
                    if uniprot_id in acs:
                        qualifier = words[2]
                        term = int(words[3].split(":")[1])
                        if term not in terms:
                            print("(skipping) term not found", term, line)
                            # if just a couple of missing terms, is probably a new annotation
                            continue
                        eco_term = words[5].split(":")[1].strip()
                        info.add((term, uniprot_id, qualifier, eco_term))
        return info

