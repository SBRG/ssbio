import unittest

from ssbio.core.object import Object

class TestObject(unittest.TestCase):
    """Unit tests for Object"""

    @classmethod
    def setUpClass(self):
        self.ob = Object(id='idtester', description='nametester')

    def test_attr(self):
        self.assertTrue(hasattr(self.ob, 'id'))
        self.assertTrue(hasattr(self.ob, 'description'))

    def test_update(self):
        new_stuff = {'newkey':'newvalue', 'dontadd':'dontadd', 'description':'newdescription'}
        self.ob.update(newdata=new_stuff, overwrite=False, only_keys=['newkey','description'])
        self.assertTrue(hasattr(self.ob, 'newkey'))
        self.assertEqual('nametester', self.ob.description)

        self.ob.update(newdata=new_stuff, overwrite=True, only_keys=['description'])
        self.assertEqual('newdescription', self.ob.description)

        self.ob.update(newdata=new_stuff, overwrite=True)
        self.assertEqual('newdescription', self.ob.description)
        self.assertEqual('newvalue', self.ob.newkey)
        self.assertEqual('dontadd', self.ob.dontadd)

    def test_get_dict(self):
        gotdict = self.ob.get_dict(only_attributes='id')
        self.assertEqual(gotdict, {'id':'idtester'})