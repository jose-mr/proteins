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

from pseudoenzymes.settings import EC_DAT_FILE, EC_CLASSES_FILE, SWISSPROT_DAT_FILE
import uniprot.models as uniprot


# class ECQuerySet(models.QuerySet):
    # """Some predefined querysets for the EC model"""

    # def current(self):
        # """Filters only current ECs. Ignores preliminary, transferred and deleted"""
        # return self.filter(is_transferred=False, is_deleted=False)

    # def get_code2id(self):
        # """Returns a dictionary of EC code to database id pairs"""
        # return {ec.code: ec.id for ec in self}

    # def codes(self, codes):
        # """Returns ec objects with these codes"""
        # return self.annotated_code().filter(Q(ec1234__in=codes) | Q(ec123__in=codes) | Q(ec12__in=codes))



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
        cls._create_uniprot_relations(number_to_acs)

    @classmethod
    def _create_uniprot_relations(cls, number_to_acs):
        relations = []
        for number, acs in number_to_acs.items():
            for ac in acs:
                relations.append(EntryUniProtEntry(entry_id=number, uniprot_entry_id=ac))
        print(f"Creating {len(relations)} uniprot associations")
        EntryUniProtEntry.objects.bulk_create(relations)


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

        
        return
        print("oi")
        objs = []
        synonyms = []
        for ec1 in etree.parse(str(EC_INTENZ_XML)).getroot():
            ec1_no = ec1.get("ec1")
            objs.append(cls._init_obj([ec1_no, 0, 0, 0], xml_obj=ec1))
            for ec2 in cls._find_children_with_tag(ec1, 'ec_subclass'):
                ec2_no = ec2.get("ec2")
                objs.append(cls._init_obj([ec1_no, ec2_no, 0, 0], xml_obj=ec2))
                for ec3 in cls._find_children_with_tag(ec2, 'ec_sub-subclass'):
                    ec3_no = ec3.get("ec3")
                    objs.append(cls._init_obj([ec1_no, ec2_no, ec3_no, 0], xml_obj=ec3))
                    for ec4 in cls._find_children_with_tag(ec3, 'enzyme'):
                        ec4_no = ec4.get("ec4")
                        objs.append(cls._init_obj([ec1_no, ec2_no, ec3_no, ec4_no],
                                                  xml_obj=ec4))
                        is_preliminary = ec4.get("preliminary") == "true"
                        number = cls.components_to_number([ec1_no, ec2_no, ec3_no, ec4_no],
                                                          is_preliminary)
                        synonym_parent = cls._find_children_with_tag(ec4, 'synonyms')
                        if synonym_parent:
                            for syn in synonym_parent[0]:
                                synonyms.append(Synonym(entry_id=number, name=syn.text))
        cls.objects.bulk_create(objs)
        Synonym.objects.bulk_create(synonyms, ignore_conflicts=True)

    @classmethod
    def download_intenz_xml_file(cls):
        """download intenz xml file with EC info
        this does not seem to be updated anymore"""
        urlretrieve("ftp://ftp.ebi.ac.uk/pub/databases/intenz/xml/ASCII/intenz.xml",
                    EC_INTENZ_XML)

    @classmethod
    def components_to_number(cls, ec_components, is_preliminary):
        return "{0}.{1}.{2}.{4}{3}".format(
                                           *[c if c else '-' for c in ec_components],
                                           "n" if is_preliminary else "")

    @classmethod
    def create_from_intenz_file(cls):
        """creates EC objects and synonyms from the intenz xml file"""
        objs = []
        synonyms = []
        for ec1 in etree.parse(str(EC_INTENZ_XML)).getroot():
            ec1_no = ec1.get("ec1")
            objs.append(cls._init_obj([ec1_no, 0, 0, 0], xml_obj=ec1))
            for ec2 in cls._find_children_with_tag(ec1, 'ec_subclass'):
                ec2_no = ec2.get("ec2")
                objs.append(cls._init_obj([ec1_no, ec2_no, 0, 0], xml_obj=ec2))
                for ec3 in cls._find_children_with_tag(ec2, 'ec_sub-subclass'):
                    ec3_no = ec3.get("ec3")
                    objs.append(cls._init_obj([ec1_no, ec2_no, ec3_no, 0], xml_obj=ec3))
                    for ec4 in cls._find_children_with_tag(ec3, 'enzyme'):
                        ec4_no = ec4.get("ec4")
                        objs.append(cls._init_obj([ec1_no, ec2_no, ec3_no, ec4_no],
                                                  xml_obj=ec4))
                        is_preliminary = ec4.get("preliminary") == "true"
                        number = cls.components_to_number([ec1_no, ec2_no, ec3_no, ec4_no],
                                                          is_preliminary)
                        synonym_parent = cls._find_children_with_tag(ec4, 'synonyms')
                        if synonym_parent:
                            for syn in synonym_parent[0]:
                                synonyms.append(Synonym(entry_id=number, name=syn.text))
        cls.objects.bulk_create(objs)
        Synonym.objects.bulk_create(synonyms, ignore_conflicts=True)

    @classmethod
    def _init_obj(cls, ecs, xml_obj=None):
        """create EC entry from the ec codes and the ec xml object"""
        name_key = 'accepted_name' if ecs[3] else 'name'
        description = cls._find_child_text_with_tag(xml_obj, 'description')
        is_preliminary = xml_obj.get("preliminary") == "true"
        comments = cls._find_children_with_tag(xml_obj, 'comments')
        comments = "\n".join([c.text for c in comments[0]]) if comments else ""
        return cls(
                name=cls._find_child_text_with_tag(xml_obj, name_key),
                number = cls.components_to_number(ecs, is_preliminary),
                is_deleted = bool(cls._find_children_with_tag(xml_obj, 'deleted')),
                is_transferred = bool(cls._find_children_with_tag(xml_obj, 'transferred')),
                is_preliminary=is_preliminary,
                description=description,
                comments=comments,
                systematic_name = cls._find_child_text_with_tag(xml_obj, 'systematic_name'),
        )

    @staticmethod
    def _find_children_with_tag(elem, tag):
        """return children that contains this tag in xml object"""
        return [child for child in elem if etree.QName(child).localname == tag]

    @classmethod
    def _find_child_text_with_tag(cls, elem, tag):
        """return the text of the child with this tag of a xml object"""
        results = cls._find_children_with_tag(elem, tag)
        if results:
            return results[0].text
        else:
            return ''


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
    def create_from_dat_file(cls):
        """find and craete all uniprot <-> ec pairs in the dat file"""
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
        created = cls.objects.bulk_create(objs, ignore_conflicts=True)
        print(len(created))
        return


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
