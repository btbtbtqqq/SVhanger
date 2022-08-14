from collections import Counter
import gzip
import argparse
import pandas as pd
import os as os

#Display 9 columes
pd.set_option('display.max_columns', 9)
#Display 50 rows
pd.set_option('display.max_rows', 10000)
VCF_COL_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']

class VCF_STAT():

    def __init__(self,vcf_file_input):
        """
        Open an optionally gzipped VCF file and generate an OrderedDict for each line.
        """
        self.filename = vcf_file_input

    def dataframe(self, large=True):
        """
        Open an optionally gzipped VCF file and return a pandas.DataFrame with
            each INFO field included as a column in the dataframe.
            Note: Using large=False with large VCF files. It will be painfully slow.
            :param filename:    An optionally gzipped VCF file.
            :param large:       Use this with large VCF files to skip the ## lines and
                        leave the INFO fields unseparated as a single column.

            example:
                vcf = VCF_control(file_name) : load VCF class
                pd_sv = vcf.dataframe(large=True) : save vcf into dataframe
        """

        if large:
            # Set the proper argument if the file is compressed.
            # 'zip','gzip','bz2','zstd'
            comp = 'gzip' if self.filename.endswith('.gz') else None
            # Count how many comment lines should be skipped.

            # Return a simple DataFrame without splitting the INFO column
            # allow to read the vcf in gz format
            self.vcf_pd = pd.read_table(self.filename, compression=comp,
                                 skip_blank_lines=True,
                                 comment='#',
                                 na_filter=True,
                                 names=VCF_COL_HEADER, usecols=range(10), index_col=None)
            return self.vcf_pd

        # Each column is a list stored as a value in this dict. The keys for this
        # dict are the VCF column names and the keys in the INFO column.

        # result =
        # (VCF_COL_HEADER) [
        # 'CHROM':[chr1,chr2,....],
        # 'POS':[2,3,4,5],
        # 'ID':[2,3,4,5],
        # 'REF':[2,3,4,5],
        # 'ALT':[2,3,4,5],
        # 'QUAL':[2,3,4,5],
        # 'FILTER':[2,3,4,5],
        # 'INFO':[2,3,4,5],
        # 'SAMPLE':[2,3,4,5]]


    def show_stat_GT(self,VCF_dataframe):
        """
        caculate the GT types and numeber
        """
        GT_list = []

        vcf_GT = VCF_dataframe
        for index, row in vcf_GT.iterrows():
            GT_NAME = str(row['INFO'])
            GT_NUM = str(row['SAMPLE'])

            GT_NUM = GT_NUM.split(':')[0]
            if GT_NAME.find('=') == -1:
                GT_NAME_snp = str(row['FORMAT']).split(':')[0]
                GT_list.append(GT_NAME_snp + ' ' + GT_NUM)
            else:
                GT_NAME = GT_NAME.split(';')[1].split('=')[1]
                GT_list.append(GT_NAME + ' ' + GT_NUM)

        count_GT = Counter(GT_list)

        for i in sorted(count_GT.items()):
            print(str(i[0])+' '+str(i[1]))


    def show_svnum(self,vcf_file):
        num_rows = int(len(vcf_file))
        return num_rows

    """
    [1]Filter VCF module
    """

    def filter_pass(self, vcf_pd):
        filter_pass = vcf_pd.query("FILTER.str.startswith('PASS').values")
        return filter_pass

    def filter_RE(self, vcf_pd):
        filter_RE = vcf_pd.query("~INFO.str.contains('RE\=1').values")
        return filter_RE

    def filter_BND(self, vcf_pd):
        filter_BND = vcf_pd.query("~INFO.str.contains('BND').values")
        return filter_BND

    def filter_haploid(self, vcf_pd):
        filter_hp = vcf_pd.query("~SAMPLE.str.contains('\.\/\.').values")
        #filter_hp = vcf_pd.query("~SAMPLE.str.contains('0/1').values")
        return filter_hp

    def filter_keep_haploid_only(self, vcf_pd):
        filter_vcf = vcf_pd.query("~SAMPLE.str.contains('\.\/\.').values")
        filter_hp = filter_vcf.query("~SAMPLE.str.contains('1|1').values")
        #filter_hp = vcf_pd.query("~SAMPLE.str.contains('0/1').values")
        return filter_hp

    """
    [2]Display information module
    """

    def show_maxlen(self,  INFO_values):
        NEWLIST = []
        df_svlen_list = INFO_values.tolist()
        for i in range(len(df_svlen_list)):
            NEWLIST.append(abs(int(df_svlen_list[i].split(';')[2].split('=')[1].strip())))
        len_max = max(NEWLIST)
        return len_max

    def show_minlen(self, INFO_values):

        NEWLIST=[]
        df_svlen_list = INFO_values.tolist()
        for i in range(len(df_svlen_list)):
            NEWLIST.append(abs(int(df_svlen_list[i].split(';')[2].split('=')[1].strip())))
        len_min = min(NEWLIST)
        return len_min

    def show_unsign(self, vcf_pd):
        len = 0
        for i in vcf_pd.SAMPLE:
            if i.split(':')[0] == '0/1':
                len += 1
            else:
                continue
        return len

    def show_BND(self, vcf_pd):
        len = 0
        for i in vcf_pd.ID:
            if i.find('BND') != -1:
                len += 1
            else:
                continue
        return len

    '''
    [3]HEAD information  module by generator 
    '''
    def head_generator(self):
        """head ## save in a generator for reading file quickly
            """
        vcf_file_open = gzip.open if self.filename.endswith('.gz') else open

        with vcf_file_open(self.filename) as vcf_file:
            for line in vcf_file:
                if line.startswith('##'):
                    yield line.strip()
                else:
                    continue
    def show_head(self):
        """print head to check
            """
        for k,v in enumerate(list(self.head_generator())):
            print(v)

    def head(self):
        """save head into list
            """
        head = list(self.head_generator())
        return (head)

    '''
    [4] Save SV length distribution by CSV module
    '''
    def len_destribution(self,vcf_pd):
        len_list = []
        chr_list = []
        for i in vcf_pd.INFO:
            len_list.append(abs(int(i.split(';')[2].split('=')[1])))
        for j in vcf_pd.CHROM:
            chr_list.append(j.strip())

        len_destr_df = pd.DataFrame({'CHR': chr_list,'SVLEN': len_list})
        print('##### Save SV length distribution #####')
        return len_destr_df

    def output_csv(self, VCF_dataframe):
        output_filename = os.path.splitext(self.filename)[0] + '_' + 'len_destribution' + '.csv'
        with open(output_filename, 'w') as vcf:
            VCF_dataframe.to_csv(output_filename, sep="\t", mode='a', index=False, header=True)
        print('##### Output SV length distribution CSV#####')
        vcf.close()



def main():
    """
    Display Statistical genotypes information and Length distribution of structural variants.
    (1)Statistical genotypes information fo structural variants from VCF file
    (2)Generate (phased_SV_length_destribution) in CSV
    (3)Display SV statistical result in Terminal
    """
    parser = argparse.ArgumentParser(description="1.Statistical genotypes information fo structural variants from VCF file" \
                                                 "2.Generate (phased_SV_length_destribution);" )

    parser.add_argument('-v', '--vcf-file', dest='sv_vcf', required=True, help='Input VCF file (structural variants)')
    args = parser.parse_args()
    print('#   Input VCF File  #')
    print(args.sv_vcf)

    #Os.path.isfile () : checks whether an object (with an absolute path) is a file
    #Os.path.isdir () : checks whether an object (with an absolute path) is a directory
    if not os.path.isfile(args.sv_vcf):
        print("ERROR: VCF file does not exist!")
        exit(-1)


    #file_name = 'your path'
    file_name = args.sv_vcf
    vcf_file = VCF_STAT(file_name)
    pd_sv = vcf_file.dataframe()

    """
    (1) show the PASS SV only
    """
    filtered_pass = vcf_file.filter_pass(pd_sv)
    print('[INFO] The PASS SV:', ' ', vcf_file.show_svnum(filtered_pass))
    vcf_file.show_stat_GT(filtered_pass)
    print('--------------------------------------------------------')


    """
    (2) show the filtered SV (Delect 'BND'/ Unphasing SV)
    """
    filtered_BND = vcf_file.filter_BND(filtered_pass)
    filtered_hp = vcf_file.filter_haploid(filtered_BND)
    print('[INFO] The filtered SV:', ' ',vcf_file.show_svnum(filtered_hp))
    vcf_file.show_stat_GT(filtered_hp)

    #show more information
    print('[INFO][SV count] The total deleted number of \'BND\' ：', vcf_file.show_BND(filtered_pass))
    print('[INFO][SV count] The total deleted number of \'Unphasing\' ：', vcf_file.show_unsign(filtered_BND))
    print('--------------------------------------------------------')
    num_1 = vcf_file.show_svnum(filtered_hp)
    num_2 = vcf_file.show_svnum(filtered_pass)
    print('[INFO][SV count] The total number of SV  (after filtering) ：', num_1)
    print('[INFO][SV count] The total number of deleted SV (after filtering) ：', (num_2 - num_1))

    """
    (3) show the filtered SV (Delect 'BND'/ Unphasing SV)
    """
    print('--------------------------------------------------------')
    print('[INFO] The longest SV ：', vcf_file.show_maxlen(filtered_hp['INFO'].values))
    print('[INFO] The shortest SV ：', vcf_file.show_minlen(filtered_hp['INFO'].values))
    """
    (4) Generate the length distribution of SV
    """
    # Output the CSV file containing SV length distribution
    len_df = vcf_file.len_destribution(filtered_hp)
    vcf_file.output_csv(len_df)


if __name__ == "__main__":
    main()
