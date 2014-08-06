all: databases blast

databases:
	doit -n 4 prep_databases

blast:
	doit -n 4 blast

clean:
	doit clean
