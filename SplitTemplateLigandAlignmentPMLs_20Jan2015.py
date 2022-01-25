#! /usr/bin/python

PML_file = 'TemplateLigandAlignmentPML_20Jan2015.txt'
input_blocks = file(PML_file).read().split("END_BLOCK\n")[1:-1]

for input_block in input_blocks:
    dir_name = input_block.split("\n")[0].split()[-1]
    print dir_name
    file(dir_name + '/' + dir_name + '.pml', 'w').write("\n".join([x for x in input_block.split("\n") if not '_COMMENT_' in x]))
