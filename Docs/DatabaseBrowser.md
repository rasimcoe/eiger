

# 

The *eiger* SQL database should primarily accessed in your code
through a python API that incorporates the psycopg2 library, as documented [here](https://github.com/rasimcoe/eiger/blob/main/Docs/EIGERDataHandling_r02.pdf), but it is extremely handy
to have a way of running queries interactively in a GUI, and having
the results returned as browsable lists.

This can be done using the free software
[pgadmin4](https://www.pgadmin.org/download/). This is the industry
standard tool used to access postgresql databases, and is well
suported in the community.

