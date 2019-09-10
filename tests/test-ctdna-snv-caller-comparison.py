import pytest
import os


class TestFilesExist:
    """ 
    Check correct output files are created.
    This is kind of redundant as Nextflow won't complete if files aren't
    made, but an extra check can't hurt.
    """

    out_path = 'data/test_data/output'
    sample_name = 'test_data'

    # processed BAM files
    def test_bam_rg(self):
        assert os.path.exists(
            f'{self.out_path}/processed_bams/{self.sample_name}_rg.bam'
        ) == True

    def test_bam_index(self):
        assert os.path.exists(
            f'{self.out_path}/processed_bams/{self.sample_name}_rg.bam.bai'
        ) == True

    # VCFs from each variant caller
    def test_vcf_mutect(self):
        assert os.path.exists(
            f'{self.out_path}/vcfs/{self.sample_name}_mutect.vcf'
        ) == True

    def test_vcf_other(self):
        assert os.path.exists(
            f'{self.out_path}/vcfs/{self.sample_name}_other.vcf'
        ) == True

    def test_combined_vcfs(self):
        assert os.path.exists(
            f'{self.out_path}/combined/{self.sample_name}_all.txt'
        ) == True



# TODO - test files are valid format
# TODO - test variants are called?
# TODO - these will require packages to be installed - use docker container?
