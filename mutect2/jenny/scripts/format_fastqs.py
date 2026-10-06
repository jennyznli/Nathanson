
import yaml
import glob
import re
import os
import argparse
from collections import defaultdict
from yaml.representer import Representer
yaml.add_representer(defaultdict, Representer.represent_dict)

def subset_R2_fastqs(FASTQS):
    r2=[]
    for file in FASTQS:
        if re.search(r'.*_R2\.fastq\.gz',os.path.basename(file)):
            r2.append(file)
    return r2

def make_info(r2):
    sample_data=defaultdict(dict)
    for fq2 in r2:
        basename=os.path.basename(fq2).removesuffix('_R2.fastq.gz')
        parts=basename.split('_')
        sample_name=parts[0]
        LB=parts[1]
        run,lane,index=parts[-3],parts[-2],parts[-1]
        fq1=fq2.replace('_R2.fastq.gz','_R1.fastq.gz')
        files=[fq1,fq2]
        sample_data[sample_name][basename]={'LB':LB,'run':run,'lane':lane,'index':index,'files':files}
    return sample_data

def flat_fastq_config(sample_data):
    """Collapse nested run entries to sample -> files (alignment.smk / project-utils shape)."""
    flat={}
    for sample,runs in sample_data.items():
        files=[]
        for entry in runs.values():
            files.extend(entry['files'])
        flat[sample]={'files':sorted(files)}
    return flat

def write_output(flat_config):
    os.makedirs('config',exist_ok=True)
    os.makedirs('input',exist_ok=True)
    with open('config/fastq.yml','w') as outfile:
        yaml.dump(flat_config,outfile,default_flow_style=False)
    with open('input/sample.list','w') as file:
        for i in sorted(flat_config.keys()):
            file.write(f'{i}\n')
    print('Created: config/fastq.yml, input/sample.list')

def main(argv=None):
    p=argparse.ArgumentParser()
    p.add_argument('--dir',nargs='+',default=['FASTQ'],help='FASTQ dir(s); pass every batch dir to build one combined config')
    args=p.parse_args()

    all_fastqs=[]
    for d in args.dir:
        found=glob.glob(f'{d}/*.fastq.gz')
        print(f"{len(found)} found in {d}.")
        all_fastqs.extend(found)
    r2=subset_R2_fastqs(all_fastqs)
    sample_data=make_info(r2)

    validate=defaultdict(set)
    for sample_key,sample_vals in sample_data.items():
        for runid_key,runid_vals in sample_vals.items():
            validate[sample_key].add(runid_vals['LB'])

    for sample_key,LB in validate.items():
        if len(LB)!=1:
            print(f'WARNING! Multiple LB values for {sample_key}: {",".join(LB)}')

    write_output(flat_fastq_config(sample_data))

if __name__=='__main__':
    main()
