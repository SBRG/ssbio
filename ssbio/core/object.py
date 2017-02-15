import os
import os.path as op
import ssbio.utils
import pickle
from copy import deepcopy
import pandas as pd
import ssbio.utils
import logging
log = logging.getLogger(__name__)


class Object(object):
    """Cobra core object with additional methods to update and get attributes"""

    def __init__(self, id=None, description=None):
        self.id = id
        self.description = description

    def __str__(self):
        return str(self.id)

    def __repr__(self):
        return "<%s %s at 0x%x>" % (self.__class__.__name__, self.id, id(self))

    def update(self, newdata, overwrite=False, only_keys=None):
        """Add/update any attributes from a dictionary.

        Args:
            newdata: Dictionary of attributes
            overwrite: If existing attributes should be overwritten if provided in newdata
            only_keys: List of keys to update

        Examples:
            >>> myobj = Object(id='hi', description='blankname')
            >>> myobj.update({'description':'withname'}, overwrite=False)
            >>> myobj.get_dict() == {'id': 'hi', 'description':'blankname'}
            True

            >>> myobj = Object(id='hi', description='blankname')
            >>> myobj.update({'description':'withname'}, overwrite=True)
            >>> myobj.get_dict() == {'id': 'hi', 'description':'withname'}
            True

            >>> myobj = Object(id='hi', description='blankname')
            >>> myobj.update({'description':'withname', 'randomkey':'randomval'}, overwrite=True)
            >>> myobj.get_dict() == {'id': 'hi', 'description': 'withname', 'randomkey': 'randomval'}
            True

            >>> myobj = Object(id='hi', description='blankname')
            >>> myobj.update({'description':'withname', 'randomkey':'randomval'}, overwrite=True, only_keys='description')
            >>> myobj.get_dict() == {'id': 'hi', 'description': 'withname'}
            True

            >>> myobj = Object(id='hi', description='blankname')
            >>> myobj.update({'description':'withname', 'randomkey':'randomval'}, overwrite=True, only_keys='randomkey')
            >>> myobj.get_dict() == {'id': 'hi', 'description': 'blankname', 'randomkey': 'randomval'}
            True

            >>> myobj = Object(id='hi', description='blankname')
            >>> myobj.update({'description':'withname', 'randomkey':'randomval'}, overwrite=False)
            >>> myobj.get_dict() == {'id': 'hi', 'description': 'blankname', 'randomkey': 'randomval'}
            True

            >>> myobj = Object(id='hi', description='blankname')
            >>> myobj.update({'description':'withname', 'randomkey':'randomval'}, overwrite=False, only_keys='randomkey')
            >>> myobj.get_dict() == {'id': 'hi', 'description': 'blankname', 'randomkey': 'randomval'}
            True

        """
        # Filter for list of keys in only_keys
        if only_keys:
            only_keys = ssbio.utils.force_list(only_keys)
            newdata = {k:v for k,v in newdata.items() if k in only_keys}

        # Update attributes
        for key, value in newdata.items():
            # Overwrite flag overwrites all attributes
            if overwrite:
                setattr(self, key, value)
            else:
                # Otherwise check if attribute is None and set it if so
                if hasattr(self, key):
                    if not getattr(self, key):
                        setattr(self, key, value)
                    else:
                        continue
                # Or just add a new attribute
                else:
                    setattr(self, key, value)

    def get_dict(self, only_keys=None, exclude_attributes=None, df_format=False):
        """Get a dictionary of this object's attributes. Optional format for storage in a Pandas DataFrame.

        Args:
            keys (str, list): Attributes that should be returned. If not provided, all are returned.
            exclude_attributes (str, list): Attributes that should be excluded.
            df_format (bool): If dictionary values should be formatted for a dataframe
                (everything possible is transformed into strings, int, or float -
                if something can't be transformed it is excluded)

        Returns:
            dict: Dictionary of attributes

        """

        # Choose attributes to return, return everything in the object if a list is not specified
        if not only_keys:
            keys = list(self.__dict__.keys())
        else:
            keys = ssbio.utils.force_list(only_keys)

        # Remove keys you don't want returned
        if exclude_attributes:
            exclude_attributes = ssbio.utils.force_list(exclude_attributes)
            for x in exclude_attributes:
                if x in keys:
                    keys.remove(x)

        # Copy attributes into a new dictionary
        df_dict = {}
        for k, orig_v in self.__dict__.items():
            if k in keys:
                v = deepcopy(orig_v)
                if df_format:
                    if v and not isinstance(v, str) and not isinstance(v, int) and not isinstance(v, float) and not isinstance(v, bool):
                        try:
                            df_dict[k] = ssbio.utils.force_string(deepcopy(v))
                        except TypeError:
                            log.warning('{}: excluding attribute from dict, cannot transform into string'.format(k))
                    else:
                        df_dict[k] = deepcopy(v)
                else:
                    df_dict[k] = deepcopy(v)
            # else:
            #     log.debug('{}: not copying attribute'.format(k))
        return df_dict

    def save_dataframes(self, outdir, prefix='df_'):
        """Save all attributes that start with "df" into a specified directory.

        Args:
            outdir (str): Path to output directory
            prefix (str): Prefix that dataframe attributes start with

        """
        # Get list of attributes that start with "df_"
        dfs = list(filter(lambda x: x.startswith(prefix), dir(self)))

        for df in dfs:
            outpath = ssbio.utils.outfile_maker(inname=df, outext='.csv', outdir=outdir)
            my_df = getattr(self, df)
            if not isinstance(my_df, pd.DataFrame):
                raise TypeError('{}: object is not a Pandas DataFrame'.format(df))

            my_df.to_csv(outpath)
            log.debug('{}: saved dataframe'.format(outpath))

        log.debug('Saved {} dataframes at {}'.format(len(dfs), outdir))

    def save_pickle(self, outname, protocol=2, outext='.pckl', outdir=None):
        """Save the object as a pickle file

        Args:
            outname (str): Basename of file
            protocol (int): Pickle protocol to use. Default is 2 to remain compatible with Python 2
            outext (str): Extension of file
            outdir (str): Path to output directory

        Returns:
            str: Path to pickle file

        """
        if not outdir:
            outdir = os.getcwd()

        outfile = ssbio.utils.outfile_maker(inname=outname, outext=outext, outdir=outdir)

        with open(outfile, 'wb') as f:
            pickle.dump(self, f, protocol=protocol)

        return outfile