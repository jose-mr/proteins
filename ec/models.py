# coding=utf-8
"""Includes the models related with the EC classification"""

# # python imports
from urllib.request import urlretrieve
from lxml import etree

# # django imports
from django.db import models

from pseudoenzymes.settings import EC_INTENZ_XML


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
    is_preliminary = models.BooleanField()
    is_deleted = models.BooleanField()
    is_transferred = models.BooleanField()
    systematic_name = models.TextField()
    comments = models.TextField()
    description = models.TextField()

    # objects = ECQuerySet.as_manager()

    @classmethod
    def download_intenz_xml_file(cls):
        """downloads intenz xml file with EC info"""
        urlretrieve("ftp://ftp.ebi.ac.uk/pub/databases/intenz/xml/ASCII/intenz.xml",
                    EC_INTENZ_XML)

    def __str__(self):
        return str(self.number)

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
