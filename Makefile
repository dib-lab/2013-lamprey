all: databases blast

wdir:
	mkdir -p work
	ln -fs ../metadata.ini work/
	ln -fs ../metadata.spec.ini work/

databases: wdir
	doit --dir work -n 4 prep_databases

blast: wdir
	doit --dir work -n 4 blast

clean: wdir
	doit clean --dir work
