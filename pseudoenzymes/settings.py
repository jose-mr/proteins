"""
Django settings for pseudoenzymes project.

Generated by 'django-admin startproject' using Django 5.0.4.

For more information on this file, see
https://docs.djangoproject.com/en/5.0/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/5.0/ref/settings/
"""

from pathlib import Path
import os

# Build paths inside the project like this: BASE_DIR / 'subdir'.
BASE_DIR = Path(__file__).resolve().parent.parent


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/5.0/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'django-insecure-(+13g()(lr1-g^o7i(3e)@&=4=bg$r-k)-2@rzexh(woykoiya'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

ALLOWED_HOSTS = []


# Application definition

INSTALLED_APPS = [
    # 'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    # 'django.contrib.sessions',
    # 'django.contrib.messages',
    # 'django.contrib.staticfiles',
    'uniprot.apps.UniprotConfig',
    'go.apps.GoConfig',
    'eco.apps.EcoConfig',
    'ec.apps.EcConfig',
    'cath.apps.CathConfig',
    'stats.apps.StatsConfig',
    'wpdb.apps.WpdbConfig',
    'django_extensions',
    'django_filters',
    'django_tables2'
]

MIDDLEWARE = [
    # 'django.middleware.security.SecurityMiddleware',
    # 'django.contrib.sessions.middleware.SessionMiddleware',
    # 'django.middleware.common.CommonMiddleware',
    # 'django.middleware.csrf.CsrfViewMiddleware',
    # 'django.contrib.auth.middleware.AuthenticationMiddleware',
    # 'django.contrib.messages.middleware.MessageMiddleware',
    # 'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'pseudoenzymes.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',

            ],
        },
    },
]

WSGI_APPLICATION = 'pseudoenzymes.wsgi.application'


# Database
# https://docs.djangoproject.com/en/5.0/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'NAME': os.environ.get("PROTEINS_DB_NAME"),
        'USER': os.environ.get("PROTEINS_ADMIN"),
        'PASSWORD': os.environ.get("PROTEINS_ADMIN_PASSWORD"),
        'HOST': 'exade7',
        'PORT': '',
    }
}



# Password validation
# https://docs.djangoproject.com/en/5.0/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/5.0/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_TZ = True



CACHES = {
    "default": {
        "BACKEND": "django.core.cache.backends.memcached.PyMemcacheCache",
        "LOCATION": "127.0.0.1:11211",
    }
}

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/5.0/howto/static-files/

STATIC_URL = 'static/'

# Default primary key field type
# https://docs.djangoproject.com/en/5.0/ref/settings/#default-auto-field

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'

# CUSTOM settings
DATA_FOLDER = BASE_DIR / "data"
SWISSPROT_DAT_FILE = DATA_FOLDER / "uniprot/uniprot_sprot.dat.gz"
SWISSPROT_ACS_FILE = DATA_FOLDER / "uniprot/acs.txt"

GO_DATA_FOLDER = DATA_FOLDER / "go"
GO_DATA_FOLDER.mkdir(parents=True, exist_ok=True)
GENE_ONTOLOGY_FILE = GO_DATA_FOLDER / "go.obo"

GO_GPA_FILE = GO_DATA_FOLDER / "goa_uniprot_all.gpa.gz"

ECO_DATA_FOLDER = DATA_FOLDER / "eco"
ECO_DATA_FOLDER.mkdir(parents=True, exist_ok=True)
ECO_ONTOLOGY_FILE = ECO_DATA_FOLDER / "eco.obo"

EC_DATA_FOLDER = DATA_FOLDER / "ec"
EC_DATA_FOLDER.mkdir(parents=True, exist_ok=True)
EC_DAT_FILE = EC_DATA_FOLDER / "enzyme.dat"
EC_CLASSES_FILE = EC_DATA_FOLDER / "enzclass.txt"
EC_INTENZ_XML = ECO_DATA_FOLDER / "intenz.xml"

CATH_DATA_FOLDER = DATA_FOLDER / "cath"
CATH_NAMES_FILE = CATH_DATA_FOLDER / "cath-b-newest-names.gz"

INTERPRO_DATA_FOLDER = DATA_FOLDER / "interpro"
INTERPRO_DAT_FILE = INTERPRO_DATA_FOLDER / "protein2ipr.dat.gz"
INTERPRO_ONLY_G3_SP_DAT_FILE = INTERPRO_DATA_FOLDER / "g3d_swissprot_only.txt"

PDB_DATA_FOLDER = DATA_FOLDER / "pdb"
PDB_DATA_FOLDER.mkdir(parents=True, exist_ok=True)
PDB_ENTRIES_IDX = PDB_DATA_FOLDER / "entries.idx"
PDB_UNIPROT_SIFTS = PDB_DATA_FOLDER / "uniprot_pdb.tsv.gz"
PDB_UNIPROT_DAT_FILE = PDB_DATA_FOLDER / "uniprot_pdb.dat.gz"

# ANALYSIS OUTPUT
OUT_FOLDER = BASE_DIR / "out"

# MSA and Trees
MSA_FOLDER = OUT_FOLDER / "msa"
MSA_FOLDER.mkdir(parents=True, exist_ok=True)

MSA_BY_DOMAIN_ARCHITECTURE = MSA_FOLDER / "domain_architecture"
MSA_BY_DOMAIN_ARCHITECTURE.mkdir(parents=True, exist_ok=True)

MSA_BY_DOMAIN_STRIP = MSA_FOLDER / "domain_strip"
MSA_BY_DOMAIN_STRIP.mkdir(parents=True, exist_ok=True)

TREE_OUT = OUT_FOLDER / "trees"
TREE_OUT.mkdir(parents=True, exist_ok=True)

NOTUNG_FOLDER = OUT_FOLDER / "notung"
NOTUNG_FOLDER.mkdir(parents=True, exist_ok=True)

# software
PROGRAMS_FOLDER = BASE_DIR / "programs"
NOTUNG_JAR = PROGRAMS_FOLDER / "notung.jar"
