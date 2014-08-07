CPUS=4

all: databases blast

wdir:
	mkdir -p work
	ln -fs ../metadata.ini work/
	ln -fs ../metadata.spec.ini work/

databases: wdir
	doit --dir work -n $(CPUS) prep_databases

blast: wdir
	doit --dir work -n $(CPUS) blast

clean: wdir
	doit clean --dir work
