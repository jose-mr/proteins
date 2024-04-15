from django.shortcuts import render
from django.views.generic.list import ListView
from django_filters.views import FilterView
from django_tables2.views import SingleTableMixin
import django_tables2 as tables

import cath.models

class SuperfamiliesTable(tables.Table):
    class Meta:
        model = cath.models.Superfamily
        fields = ("number", "name", "uniprot_entries_count",
                  "uniprot_entries_ec_count", "uniprot_entries_no_ec_count")
        per_page = 50

import django_filters
class SuperfamilyFilter(django_filters.FilterSet):
    class Meta:
        model = cath.models.Superfamily
        fields = {'name': ["contains", ],
                  'number': ["startswith",]}

    @property
    def qs(self):
        queryset = super().qs.superfamilies().annotate_uniprot_count()\
                             .annotate_uniprot_enzyme_ec_count()
        return queryset


    # queryset = cath.models.Superfamily.objects.superfamilies()


class SuperfamiliesView(SingleTableMixin, FilterView):
    """class based view for the rules page"""

    # template_name = 'common/rules/rules.html'
    table_class = SuperfamiliesTable
    model = cath.models.Superfamily
    filterset_class = SuperfamilyFilter

    # def get_queryset(self):
        # queryset = self.queryset\
            # .annotate(step_count=Count("step_rules", filter=Q(step_rules__complete_in_step=True), distinct=True))\
            # .annotate(partial_step_count=Count("step_rules", filter=Q(step_rules__complete_in_step=False), distinct=True))\
            # .annotate(entry_count=Count("entry_rules__entry_id", distinct=True))\
            # .annotate(exclusive_arrow=Max(Cast('entry_rules__exclusive_in_entry', IntegerField())))\
            # .prefetch_related("entry_rules")\
            # .order_by("-step_count", "-partial_step_count")
        # if (entry_macie_id := self.request.GET.get("entry")) is not None:
            # queryset = queryset.filter(arrows__step__mechanism__reaction__entry__macie_id=entry_macie_id)
        # if (version := self.request.GET.get("version")) is not None:
            # queryset = queryset.filter(version=version)
        # else:
            # queryset = queryset.filter(version="0.2")
        # # queryset = [rule for rule in queryset if "#16" in rule.reaction_smarts]
        # return queryset

