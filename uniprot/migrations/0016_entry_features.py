# Generated by Django 5.0.4 on 2025-06-26 10:26

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('uniprot', '0015_remove_entry_taxid'),
    ]

    operations = [
        migrations.AddField(
            model_name='entry',
            name='features',
            field=models.JSONField(default=[]),
            preserve_default=False,
        ),
    ]
