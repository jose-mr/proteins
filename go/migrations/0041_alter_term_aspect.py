# Generated by Django 5.0.4 on 2025-03-17 11:53

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('go', '0040_alter_term_aspect'),
    ]

    operations = [
        migrations.AlterField(
            model_name='term',
            name='aspect',
            field=models.TextField(choices=[('molecular_function', 'molecular_function'), ('biological_process', 'biological_process'), ('cellular_component', 'cellular_component')]),
        ),
    ]
