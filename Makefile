CPUS=4

all: databases blast

wdir:
	mkdir -p work
	ln -fs ../metadata.ini work/
	ln -fs ../metadata.spec.ini work/

get_databases: wdir
	doit --dir work -n $(CPUS) get_databases


gunzip_databases: wdir
	doit --dir work -n $(CPUS) gunzip_databases

fix_database_names: wdir
	doit --dir work -n $(CPUS) fix_database_names

make_blastdbs: wdir
	doit --dir work -n $(CPUS) make_blastdbs


blast: wdir
	doit --dir work -n $(CPUS) blast

clean: wdir
	doit clean --dir work
