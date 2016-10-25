import unittest

import ssbio.utils

class TestUtils(unittest.TestCase):
    """Unit tests for utils
    """

    def test_split_folder_and_path(self):
        test_path = '/this/is/a/hypothetical/path/to_a_file.myfile'
        result = ('/this/is/a/hypothetical/path', 'to_a_file', '.myfile')

        self.assertEqual(result, ssbio.utils.split_folder_and_path(test_path))