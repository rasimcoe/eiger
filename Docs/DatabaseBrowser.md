

# Introduction, Purpose and Scope

The *eiger* SQL database should primarily accessed in your code
through a python API that incorporates the psycopg2 library, as documented [here](https://github.com/rasimcoe/eiger/blob/main/Docs/EIGERDataHandling_r02.pdf), but it is extremely handy
to have a way of running queries interactively in a GUI, and having
the results returned as browsable lists.

This can be done using the free software
[pgadmin4](https://www.pgadmin.org/download/). This is the industry
standard tool used to access postgresql databases, and is well
suported in the community. Installation instructions are included with
the install package, this document explains how to configure pgadmin4
for use with the *eiger* AWS database server, and common ways to run
queries.

For interactive queries the password structure can be somewhat
cumbersome, because we have set up the database using IAM credentials
fom AWS, and these use machine generated passwords that are changed
every 10 minutes.

# Refreshing your AWS IAM credentials 

In the process of configuring the eiger database package, you should have already installed the awstools package. This contains a set of command line tools that you can use to generate new password sets with AWS. The most convenient way to do this is to set up an alias in your .bash_profile or equivalent configuration file, I call this alias "eigerpass" on my system.  Note that this requires your EIGER_DBSERVER environment variable to be set.

> alias eigerpass="aws rds generate-db-auth-token --hostname $EIGERDB_SERVER --port 5432 --region us-east-2 --username $USER | pbcopy" 

In this case, calling "eigerpass" from the command line generates a new password and automatically copies it to your clipboard (like clicking Command-C).

# Configuring **pgadmin4** for use with eiger's AWS database

We assume that your have downloaded and installed the package,
e.g. using the .dmg file for apple workstations.




