import cruzdb, sys

db = sys.argv[1]
db_file = "sqlite:////tmp/%s.db" % (db,)
cruzdb.Genome(db).mirror(["refGene", "knownGene"], db_file)
