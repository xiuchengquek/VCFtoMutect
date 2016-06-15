
import unittest
import sys

if sys.version[0] == 3:
    import unittest.mock as mock
    builtin_open = 'builtins.open'


else:
    builtin_open = '__builtin__.open'
    import mock


from VariantFileObject import VariantFileObject, Variant



class TestVariantObjectt(unittest.TestCase):
    def test_variant_is_snp(self):
        v = Variant(1,10,1,'A','C')
        self.assertEqual(v.chr, 1)
        self.assertEqual(v.pos, 10)
        self.assertEqual(v.id, 1)
        self.assertEqual(v.alt, 'C')
        self.assertEqual(v.ref, 'A')
        self.assertEqual(v.snp, True)

    def test_variant_is_not_snp(self):
        v = Variant(1,10,1,'A','CT')
        self.assertEqual(v.chr, 1)
        self.assertEqual(v.pos, 10)
        self.assertEqual(v.id, 1)
        self.assertEqual(v.alt, 'CT')
        self.assertEqual(v.ref, 'A')
        self.assertEqual(v.snp, False)

    def test_variant_is_not_snp_deletion(self):
        v = Variant(1,10,1,'CT','C')
        self.assertEqual(v.chr, 1)
        self.assertEqual(v.pos, 10)
        self.assertEqual(v.id, 1)
        self.assertEqual(v.alt, 'C')
        self.assertEqual(v.ref, 'CT')
        self.assertEqual(v.snp, False)

    def test_variant_equality(self):
        v1 = Variant(1,10,1,'CT','C')
        v2 = Variant(1,10,1,'CT','C')
        self.assertTrue(v1 == v2)




    def test_to_bed_format_default(self):
        v = Variant(1,10,1,'CT','C')
        self.assertEqual(v.toBedFormat(), '1\t10\t10\t1')


class TestVariantFileObject(unittest.TestCase):
    def setUp(self):
        self.vcf_content = [
            '#',
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER',
            '1\t10\t-\tA\tC\t10\tPASS',
            '2\t10\t-\tA\tC\t10\tPASS',
            '3\t10\t-\tAT\tC\t10\tPASS'

        ]

        self.test_variant_objects  = [
            Variant('1',10 , 1 , 'A', 'C' ),
            Variant('2',10 , 2 , 'A', 'C' ),
            Variant('3',10 , 3 , 'AT', 'C' )
        ]



    def test_load_variant_file(self):
        v = VariantFileObject(Variant)
        m = mock.mock_open(read_data='\n'.join(self.vcf_content))
        m.return_value.__iter__ = lambda self: self
        m.return_value.__next__ = lambda self: self.readline()
        with mock.patch(builtin_open,m ):
            v.load_variant_file('test.vcf')
            self.assertEqual(len(v._variant) , 3)
            self.assertEqual(v._variant[0], self.test_variant_objects[0])
            self.assertEqual(v._variant[1], self.test_variant_objects[1])
            self.assertEqual(v._variant[2], self.test_variant_objects[2])

    def test_drop_indels(self):
        v = VariantFileObject(Variant)
        m = mock.mock_open(read_data='\n'.join(self.vcf_content))
        m.return_value.__iter__ = lambda self: self
        m.return_value.__next__ = lambda self: self.readline()
        with mock.patch(builtin_open,m ):
            v.load_variant_file('test.vcf')
            self.assertEqual(len(v._variant) , 3)
            v.drop_indels()
            self.assertEqual(len(v._variant) , 2)

    def test_write_bed_file(self):
        v = VariantFileObject(Variant)
        m = mock.mock_open(read_data='\n'.join(self.vcf_content))
        m.return_value.__iter__ = lambda self: self
        m.return_value.__next__ = lambda self: self.readline()

        with mock.patch(builtin_open,m ):
            v.load_variant_file('test.vcf')
            v.drop_indels()
            fh = mock.mock_open()
            with mock.patch(builtin_open, fh):
                v.write_bed_file('output.bed')
                calls = mock.call
                calls_list = [calls('output.bed' ,'w')]
                fh.assert_has_calls(calls_list, any_order=True)
                handle = fh()
                calls_list = [
                    calls("chr\tstart\tend\tid\n"),
                    calls("1\t10\t10\t1\n",),
                    calls("2\t10\t10\t2\n",)
                ]
                handle.write.assert_has_calls(calls_list)



                v.write_bed_file('output_with_penalty.bed', -3, 3)
                calls = mock.call
                calls_list = [calls('output_with_penalty.bed' ,'w')]
                fh.assert_has_calls(calls_list, any_order=True)
                handle = fh()
                calls_list = [
                    calls("chr\tstart\tend\tid\n"),
                    calls("1\t7\t13\t1\n",),
                    calls("2\t7\t13\t2\n",)
                ]
                handle.write.assert_has_calls(calls_list)




if __name__ == '__main__':
    unittest.main()
