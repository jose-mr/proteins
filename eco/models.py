# python imports
import re
from urllib.request import urlretrieve
from collections import Counter

# django imports
from django.db import models

# pseudoenzymes imports
from pseudoenzymes.settings import ECO_ONTOLOGY_FILE

class TermQuerySet(models.QuerySet):
    """Some predefined querysets for the Eco Term model"""

    def experimental(self):
        """return all experimental eco terms"""
        return self.children_of(Term.objects.filter(name="experimental evidence"))

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
    """Evidence Ontology Terms"""
    # ids are not unique as integers. see 000156 for example
    id = models.CharField(primary_key=True, max_length=7)
    name = models.TextField(blank=False)
    definition = models.TextField()

    objects = TermQuerySet.as_manager()

    def __repr__(self):
        return self.code

    def __str__(self):
        return f"{self.id} - {self.name}"

    @property
    def code(self):
        """Returns a eco term in this format ECO:0001234"""
        return f"ECO:{self.id}"

    @classmethod
    def download_ontology(cls):
        """download eco ontology file"""
        urlretrieve("https://raw.githubusercontent.com/evidenceontology/"
                    "evidenceontology/master/eco.obo", ECO_ONTOLOGY_FILE)

    @classmethod
    def create_from_ontology_file(cls):
        """Reads the ontology file and adds all the ECO terms to the database"""
        inside_quotes = re.compile(r'"[^"\\]*(?:\\.[^"\\]*)*"')
        objs = []
        info = {}
        with open(ECO_ONTOLOGY_FILE, 'r') as obo_file:
            for line in obo_file:
                if line == "\n" and "id" in info:
                    objs.append(cls(**info))
                    info = {}
                elif line.startswith("id: ECO:"):
                    info["id"] = line.strip().replace("id: ECO:", "")
                elif line.startswith("def:"):
                    info["definition"] = inside_quotes.findall(line)[0][1:-1]\
                        .replace("\\\"", "\"").replace("\\n", "\n")
                elif line.startswith("name:"):
                    info["name"] = line.split(":", maxsplit=1)[1].strip()
        cls.objects.bulk_create(objs)
        print(f"Creating {len(objs)} ECO terms")


class Relation(models.Model):
    """Relations between eco terms"""
    term1 = models.ForeignKey(
            "Term",
            on_delete=models.CASCADE,
            related_name="relation1"
            )
    term2 = models.ForeignKey(
            "Term",
            on_delete=models.CASCADE,
            related_name="relation2"
            )
    relation = models.CharField(max_length=255, db_index=True)

    def __repr__(self):
        return f"{self.term1} {self.relation} {self.term2}"

    def __str__(self):
        return self.__repr__()

    class Meta:
        unique_together = ['term1', 'relation', 'term2']

    @classmethod
    def create_from_ontology_file(cls):
        """read the ontology file and add all the is_a ECO relations to the database"""
        objs = []
        with open(ECO_ONTOLOGY_FILE, 'r') as obo_file:
            term1 = None
            for line in obo_file:
                if line.startswith("[Term]"):
                    term1 = "reading"
                elif line.startswith("id: ECO") and term1 == "reading":
                    term1 = line.replace("id: ECO:", "").strip()
                elif line.startswith("is_a: ECO:") and term1 not in [None, "reading"]:
                    term2 = line.replace("is_a: ECO:", "").split("!")[0].split("{")[0].strip()
                    objs.append(cls(term1_id=term1, term2_id=term2, relation="is_a"))
                elif line == "\n":
                    term1 = None
        cls.objects.bulk_create(objs)
        print(f"Creating {len(objs)} ECO terms relations")


