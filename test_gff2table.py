#!/usr/bin/env python3

import pandas
import unittest
import gff2table

class AttributeParserTest(unittest.TestCase):
    def test_tokenizer_space_sep(self):
        attribute_parser = gff2table.AttributesParser()

        arg = 'a "b"; c 3'
        expected = ['a', ' ', '"b"', ';', 'c', ' ', '3']
        self.assertEqual(list(attribute_parser.tokenize(arg)), expected)

        arg = 'alpha 3; beta "NULL"; delta "asdf;fdsa";'
        expected = ['alpha', ' ', '3', ';',
                    'beta', ' ', '"NULL"', ';',
                    'delta', ' ', '"asdf;fdsa"', ';']
        self.assertEqual(list(attribute_parser.tokenize(arg)), expected)

        arg = 'alpha 3;beta "NULL";   delta "asdf;asdf"  '
        expected = ['alpha', ' ', '3', ';',
                    'beta', ' ', '"NULL"', ';',
                    'delta', ' ', '"asdf;asdf"']
        self.assertEqual(list(attribute_parser.tokenize(arg)), expected)

        arg = 'alpha 3;beta "NULL";   delta "asdf;xzyz"; '
        expected = ['alpha', ' ', '3', ';',
                    'beta', ' ', '"NULL"', ';',
                    'delta', ' ', '"asdf;xzyz"', ';']
        self.assertEqual(list(attribute_parser.tokenize(arg)), expected)

    def test_tokenizer_equals_sep(self):
        attribute_parser = gff2table.AttributesParser(sep='=')

        arg = 'a="b"; c=3'
        expected = ['a', '=', '"b"', ';', 'c', '=', '3']
        self.assertEqual(list(attribute_parser.tokenize(arg)), expected)

        arg = 'alpha = 3; beta = "NULL"; delta = "asdf;fdsa";'
        expected = ['alpha', '=', '3', ';',
                    'beta', '=', '"NULL"', ';',
                    'delta', '=', '"asdf;fdsa"', ';']
        self.assertEqual(list(attribute_parser.tokenize(arg)), expected)

        arg = 'alpha=3;beta="NULL";   delta="asdf;asdf"  '
        expected = ['alpha', '=', '3', ';',
                    'beta', '=', '"NULL"', ';',
                    'delta', '=', '"asdf;asdf"']
        self.assertEqual(list(attribute_parser.tokenize(arg)), expected)

        arg = 'alpha=3;beta="NULL";   delta="asdf = xzyz;"; '
        expected = ['alpha', '=', '3', ';',
                    'beta', '=', '"NULL"', ';',
                    'delta', '=', '"asdf = xzyz;"', ';']
        self.assertEqual(list(attribute_parser.tokenize(arg)), expected)
        
    def test_call(self):
        attribute_parser = gff2table.AttributesParser()
        arg = 'alpha 3;beta "NULL";   delta "asdf"; '
        self.assertEqual(attribute_parser(arg), 3)
        self.assertEqual(attribute_parser.index, 1)
        self.assertEqual(attribute_parser.max_string['delta'], 4)
        
        arg = 'alpha 3;beta "Foo";   delta "asdf;asdf"; gamma 7'
        self.assertEqual(attribute_parser(arg), 4)
        self.assertEqual(attribute_parser.index, 2)
        self.assertEqual(attribute_parser.max_string['delta'], 9)

        arg = 'alpha 3; delta "asdf;asdf;fdsa"'
        self.assertEqual(attribute_parser(arg), 2)
        self.assertEqual(attribute_parser.index, 3)
        self.assertEqual(attribute_parser.max_string['delta'], 14)

        df = pandas.DataFrame(attribute_parser.terms)
        self.assertEqual(df.shape, (3,4))
        self.assertTrue(pandas.isnull(df['beta'][0]))
        self.assertTrue(pandas.isnull(df['beta'][2]))
        
        
if __name__ == '__main__':
    unittest.main()
