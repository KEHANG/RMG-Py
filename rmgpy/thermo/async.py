import sys
import traceback


import logging as logging


class WorkerWrapper(object):
    __name__ = 'WorkerWrapper'

    def __init__(self, myfn):
        self.myfn = myfn

    def __call__(self, *args, **kwargs):
        try:
            return self.myfn(*args, **kwargs)
        except:
            type, value, tb = sys.exc_info()
            lines = traceback.format_exception(type, value, tb)
            print ''.join(lines)
            raise

class AbstractEngine(object):
    """Template for asynchronous computations that need to be stored."""

    def __init__(self):
        self.db = {}#dictionary of identifier -> values
    

    def submit_future(self, identifier):
        """
        Method needs to be overridden in subclass.
        Submits a future with the given identifier by attaching 
        the evaluator method that should be called on the remote worker.

        """
        pass


    def isValidID(self, identifier):
        """
        Method needs to be overridden in subclass.
        """
        return True

    def parse(self, identifier):
        """
        Method needs to be overridden in subclass.
        """
        return identifier


    def submit(self, obj):
        """
        Checks whether the provided object
        is already present in the database.


        If present, the associated value is returned right away.

        If not present, it is added to the database,
        and value is computed.

        """
        logging.debug("Submitting %s...", obj)

        if not self.isValidID(obj):
            logging.error("%s is not a valid identifier. Exiting...", obj)
            raise Exception

        key = self.parse(obj)

        if key in self.db:
            logging.debug("Value found for identifier %s. Not computing any data.", obj)
        else:
            logging.debug("Identifier %s not found in database. Computing data...", obj)
            data = self.submit_future(obj)
            self.call_back(data)
    

    def status(self):
        """Reports the state of the engine."""
        logging.debug("These are the keys of the db: %s", self.db.keys())

    def update_db(self, key, data):
        """Adds data to database, if necessary."""
        logging.debug("Task for the key %sis finished. Updating the database...", key)
        
        if key in self.db:
            logging.debug("Key %s already in database. Not adding it again...", key)
        else:
            self.db[key] = data

        self.status()

    def search_db(self, key):
        """Searches the database and returns the found data, if possible."""
        try:
            data = self.db[key]
            logging.debug("Identifier %s found in the database!", key)
            return data
        except KeyError, e:
            logging.debug("Could not find the key %s in the database...", key)
            return None

    def search(self, key):
        """Searches in the database"""
        
        logging.debug("Searching db...")
        data = self.search_db(key)

        if data is not None:
            return data

        return None

    def get_data(self, obj):
        """
        Retrieves the data object associated
        with the identifier.

        """
    
        key = self.parse(obj)
        
        data = self.search(key)

        if data is not None:
            self.update_db(key, data)
            return data
        
        logging.debug("Submitting %s, since it was not found yet...", key)
        self.submit(obj)
        data = self.search(key)
        
        if data is not None:
            self.update_db(key, data)
            return data
        else:
            logging.error("Data for key %s is None...", key)
            raise Exception

        return None