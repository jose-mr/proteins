# Generated by Django 5.0.4 on 2024-04-09 15:14

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('ec', '0002_entry_is_deleted_entry_is_preliminary_and_more'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='entry',
            name='id',
        ),
        migrations.AlterField(
            model_name='entry',
            name='number',
            field=models.TextField(primary_key=True, serialize=False),
        ),
    ]
