import pandas


class RSEMCache:
    def __init__(self, cache_name='rsem-genes.h5'):
        self.cache_name = cache_name

    def open(self, cache_name=None):
        if cache_name is not None:
            self.cache_name = cache_name
        self.store = pandas.HDFStore(self.cache_name, 'r')

    def close(self):
        self.store.close()

    def experiments():
        """Return experiments
        """

