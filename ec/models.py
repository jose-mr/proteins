# coding=utf-8
"""Includes the models related with the EC classification"""

# # python imports
from Bio import SeqIO
from collections import defaultdict
import gzip
from urllib.request import urlretrieve
from lxml import etree

# # django imports
from django.db import models

from pseudoenzymes.settings import (EC_DAT_FILE, EC_CLASSES_FILE, 
                                    SWISSPROT_DAT_FILE, EC_INTENZ_XML
                                    )
import uniprot.models as uniprot


# class ECQuerySet(models.QuerySet):
    # """Some predefined querysets for the EC model"""

    # def current(self):
        # """Filters only current ECs. Ignores preliminary, transferred and deleted"""
        # return self.filter(is_transferred=False, is_deleted=False)

    # def get_code2id(self):
        # """Returns a dictionary of EC code to database id pairs"""
        # return {ec.code: ec.id for ec in self}


class Entry(models.Model):
    """enzyme commission number"""
    number = models.TextField(primary_key=True)
    name = models.TextField()
    is_preliminary = models.BooleanField(default=False)
    is_deleted = models.BooleanField(default=False)
    transferred = models.TextField()
    # systematic_name = models.TextField()
    description = models.TextField()
    uniprot_entries = models.ManyToManyField(
            "uniprot.Entry",
            through="EntryUniProtEntry",
            related_name="ec_entries"
            )

    # objects = ECQuerySet.as_manager()

    def __str__(self):
        return str(self.number)

    @classmethod
    def download_classes_file(cls):
        """download sib classes file"""
        urlretrieve("https://ftp.expasy.org/databases/enzyme/enzclass.txt",
                    EC_CLASSES_FILE)

    @classmethod
    def download_ec_dat_file(cls):
        """download ec dat file with EC info"""
        urlretrieve("https://ftp.expasy.org/databases/enzyme/enzyme.dat",
                    EC_DAT_FILE)

    @classmethod
    def create_classes_classes_files(cls):
        """creates ec entries form and synonyms from the intenz xml file"""
        classes = []
        with open(EC_CLASSES_FILE, "r") as classes_file:
            for line in classes_file:
                if len(line) > 1 and line[1] == ".":
                    number = line[:10].replace(" ", "")
                    description = line[10:].strip()
                    classes.append(cls(number=number, description=description))
        print(f"Creating {len(classes)} classes, subclasses, and subsubclasses")
        cls.objects.bulk_create(classes)

    @classmethod
    def create_entries_from_dat_file(cls):
        """creates ec entries and links to uniprot from the ec dat file"""
        ec_info,  number_to_acs = cls._read_info_from_dat_file()
        objs = []
        for info in ec_info:
            if info["name"].startswith("Transferred entry:"):
                info["transferred"] = info["name"]
                info["name"] = ""
            if info["name"].startswith("Deleted entry."):
                info["is_deleted"] = True
                info["name"] = ""
            if "n" in info["number"]:
                info["is_preliminary"] = True
            objs.append(cls(**info))
        print(f"Creating {len(objs)} EC complete entries")
        cls.objects.bulk_create(objs)

    @classmethod
    def _read_info_from_dat_file(cls):
        """creates ec entries and links to uniprot from the ec dat file"""

        with open(EC_DAT_FILE, "r") as dat_file:
            objs = []
            info = defaultdict(str)
            all_info = []
            number_to_acs = defaultdict(set)
            for line in dat_file:
                if line.startswith("ID"):
                    info["number"] = line.split()[1] 
                elif line.startswith("DE"):
                    if info["name"]:
                        space = "" if info["name"].endswith("-") else " "
                        info["name"] = f"{info['name']}{space}{line[5:].strip()}"
                    info["name"] += line[5:].strip()
                elif line.startswith("CC"):
                    info["description"] += line[9:]
                elif line.startswith("DR"):
                    number_to_acs[info["number"]].update(
                        [ac_gene.split(",")[0].strip()
                         for ac_gene in line[5:-2].split(";")])
                elif line.startswith("//"):
                    all_info.append(info.copy())
                    info = defaultdict(str)
        return all_info, number_to_acs

    @classmethod
    def download_intenz_xml_file(cls):
        """download intenz xml file with EC info
        this does not seem to be updated anymore"""
        urlretrieve("ftp://ftp.ebi.ac.uk/pub/databases/intenz/xml/ASCII/intenz.xml",
                    EC_INTENZ_XML)

    @classmethod
    def create_synonyms_from_intenz_file(cls):
        """create EC synonyms from the intenz xml file

        this file is not updated anymore, so using it only for synonyms (which are
        not on the sibs file
        """
        synonyms = []
        for ec1 in etree.parse(str(EC_INTENZ_XML)).getroot():
            ec1_no = ec1.get("ec1")
            for ec2 in cls._find_children_with_tag(ec1, 'ec_subclass'):
                ec2_no = ec2.get("ec2")
                for ec3 in cls._find_children_with_tag(ec2, 'ec_sub-subclass'):
                    ec3_no = ec3.get("ec3")
                    for ec4 in cls._find_children_with_tag(ec3, 'enzyme'):
                        ec4_no = ec4.get("ec4")
                        if ec4.get("preliminary") == "true":
                            ec4_no = f"n{ec4_no}"
                        number = ".".join([ec1_no, ec2_no, ec3_no, ec4_no])
                        synonym_parent = cls._find_children_with_tag(ec4, 'synonyms')
                        if synonym_parent:
                            for syn in synonym_parent[0]:
                                synonyms.append(Synonym(entry_id=number, name=syn.text))
        Synonym.objects.bulk_create(synonyms, ignore_conflicts=True)

    @staticmethod
    def _find_children_with_tag(elem, tag):
        """return children that contains this tag in xml object"""
        return [child for child in elem if etree.QName(child).localname == tag]


class EntryUniProtEntry(models.Model):
    """Through table to link ec entries with UniProt entries"""
    entry = models.ForeignKey("Entry", on_delete=models.CASCADE)
    uniprot_entry = models.ForeignKey("uniprot.Entry", on_delete=models.CASCADE)
    # eco_term = models.ForeignKey("eco.Term", on_delete=models.CASCADE)

    # objects = SequenceGoQuerySet.as_manager()

    # # noinspection PyMissingOrEmptyDocstring
    class Meta:
        unique_together = ["entry", "uniprot_entry"]

    @classmethod
    def create_from_uniprot_dat_file(cls):
        """find and create all uniprot <-> ec pairs in the dat file

        this file is more complete than the ec dat file from sibs"""
        objs = []
        entries = set(Entry.objects.all().values_list("number", flat=True))
        uniprot_entries = set(uniprot.Entry.objects.all().values_list("ac", flat=True))
        with gzip.open(SWISSPROT_DAT_FILE, "rb") as dat_file:
            for record in SeqIO.parse(dat_file, "swiss"):
                for word in record.description.split():
                    if word.startswith("EC="):
                        number = word[3:].strip(";")
                        if number not in entries:
                            print(number)
                        if record.id not in uniprot_entries:
                            print(record.id)
                        objs.append(cls(entry_id=number, uniprot_entry_id=record.id))
        print(f"Creating {len(objs)} Uniprot<->EC associations")
        cls.objects.bulk_create(objs, ignore_conflicts=True)


class Synonym(models.Model):
    """synonyms of Enzyme names"""
    entry = models.ForeignKey(
            'Entry',
            on_delete=models.CASCADE,
            related_name="synonyms"
    )
    name = models.TextField()

    class Meta:
        unique_together = ('entry', 'name')
