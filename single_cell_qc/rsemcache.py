import pandas
from collections.abc import Mapping

class RSEMCache(Mapping):
    def __init__(self, cache_name='rsem-genes.h5'):
        self._cache_name = cache_name
        self._store = None
        self._experiments = None
        self._libraries = None

        self.open(self._cache_name)

    def open(self, cache_name=None):
        if cache_name is not None:
            self._cache_name = cache_name
        self._store = pandas.HDFStore(self._cache_name, 'r')

    def close(self):
        self._store.close()

    @property
    def experiments(self):
        """Return experiments
        """
        if self._experiments is None:
            self.update_metadata()

        return self._experiments

    @property
    def libraries(self):
        if self._libraries is None:
            self.update_metadata()
        return self._libraries

    def update_metadata(self):
        self._experiments = {}
        self._libraries = {}
        for i, row in self._store['metadata'].iterrows():
            lib_id = str(row.library_id)
            name = row.experiment_name
            self._libraries[lib_id] = name
            self._experiments.setdefault(name, []).append(lib_id)

    def __getitem__(self, key):
        return self._store[key]

    def __setitem__(self, key, value):
        raise NotImplemented("The Cache is read only")

    def __delitem__(self, key):
        raise NotImplemented("The cache is read only")

    def __iter__(self):
        for name in self._store.keys():
            if name != '/metadata':
                yield name

    def __len__(self):
        # hide metadata table
        return len(self._store) - 1

    def get_gene_expression(self,
                            replicates,
                            gene_ids=None,
                            quantification='FPKM'):
        results = {}
        for lib_id in replicates:
            lib = self._store['/genes/library_{}'.format(lib_id)]
            lib.index = lib['gene_id']
            results[str(lib_id)] = lib['FPKM']

        df =  pandas.DataFrame(results)
        if gene_ids:
            df = df.loc[list(gene_ids)]
        return df                

def test():
    cache = RSEMCache()
    names = list(cache.keys())
    assert '/metadata' not in names
    print('names:', len(names))
    table = cache[names[0]]
    print('{} shape: {}'.format(names[0], table.shape))
    print('experiments:', len(cache.experiments))
    print('libraries:', len(cache.libraries))
    cache.close()

if __name__ == '__main__':
    test()
