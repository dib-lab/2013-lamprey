CPU=4

all: databases blast

databases:
	doit -n $(CPU) prep_databases

blast:
	doit -n $(CPU) blast

clean:
	doit clean
