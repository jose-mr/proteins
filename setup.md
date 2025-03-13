# Setup database (only do this if database not created already)
the secrets file contains the environment variables required to connect to the database

create a database and user using these commands in PSQL.
replace the XXX with the variables in the secrets file
replace utc timezone for your timezone

```
CREATE DATABASE XXXX;
CREATE USER XXX_ADMIN WITH ENCRYPTED PASSWORD 'XXXX_ADMIN_PASSWORD';
ALTER ROLE XXX_ADMIN SET client_encoding TO 'utf8';
ALTER ROLE XXX_ADMIN SET default_transaction_isolation TO 'read committed';
ALTER ROLE XXX_ADMIN SET timezone TO 'UTC';
GRANT ALL PRIVILEGES ON DATABASE XXXX TO XXX_ADMIN;
ALTER DATABASE XXXX OWNER TO XXX_ADMIN;


CREATE USER XXX_GUEST WITH ENCRYPTED PASSWORD 'XXXX_GUEST_PASSWORD';

```

# install miniconda and create a new environment
```
conda create --name proteins python=3.12
```

install dependencies
```
pip install -r requirements.txt
```

